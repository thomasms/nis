"""
Fission pulse decay heat: the fission products of one fission, then decay (scripts/nisbench/pulse.py).

These check the fission yields, the decay chains, the decay energies and the solver together, with no
experiment needed:
    - the decay heat over all time is the delayed beta and gamma energy of ENDF MF=1 MT=458
    - the time steps of the solver give the same integral as a direct solve
    - two decay libraries give about the same decay heat with the same yields
    - a short irradiation through the reaction matrix gives the same heat per fission as the pulse

Needs data/jeff311.nfy (scripts/get_fission_yields.py) and the decay data from git lfs, else the tests skip.
"""

import math
import unittest

import numpy as np

import pynis as nis

from .testdata import ENDFBVII_DECAY, JEFF311_DECAY, JEFF311_NFY, JEFF4_REACTION, has_data

from nisbench.pulse import DEFAULT_TIMES, delayed_energy, pulse_decay_heat, pulse_inventory, released_energy

# ENDF/B-VII.1 MF=1 MT=458, the energy independent terms: the delayed beta (EB) and gamma (EGD) energy per
# fission in MeV, and their uncertainties
MT458 = {
    "U233": (5.16, 0.06, 5.01, 0.06),
    "U235": (6.50, 0.05, 6.33, 0.05),
    "Pu239": (5.31, 0.06, 5.17, 0.06),
    "Pu241": (6.58, 0.09, 6.40, 0.09),
}

# How close the integral must be to MT458. With the JEFF-3.1.1 decay data it is within 6 % (the beta energy is
# 3-6 % higher, the gamma energy 2-6 % lower), with the ENDF/B-VII.1 decay data within 4.3 %.
MT458_REL_TOL = {"jeff311": 0.08, "endfbvii1": 0.05}


def _decay(path):
    return nis.core.decay_mat_from_file(str(path)), nis.core.decay_energy_data_from_file(str(path))


@unittest.skipUnless(has_data(JEFF311_NFY) and has_data(JEFF311_DECAY), "needs data/jeff311.nfy and data/jeff311.dec")
class FissionPulseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.yields = nis.core.fission_yields_from_file(str(JEFF311_NFY))
        cls.libraries = {"jeff311": _decay(JEFF311_DECAY)}
        if has_data(ENDFBVII_DECAY):
            cls.libraries["endfbvii1"] = _decay(ENDFBVII_DECAY)

    def test_inventory_is_one_fission(self) -> None:
        for parent in MT458:
            with self.subTest(parent=parent):
                inventory = pulse_inventory(self.yields, parent)
                self.assertAlmostEqual(sum(inventory.values()), 2.0, delta=0.01)
                self.assertNotIn(parent, inventory)
        with self.assertRaises(ValueError):
            pulse_inventory(self.yields, "Fe56")

    def test_delayed_energy_is_mt458(self) -> None:
        for library, (decay_mat, energy) in self.libraries.items():
            for parent, (eb, _, egd, _) in MT458.items():
                with self.subTest(library=library, parent=parent):
                    beta, gamma = delayed_energy(decay_mat, energy, pulse_inventory(self.yields, parent))
                    self.assertAlmostEqual(beta / eb, 1.0, delta=MT458_REL_TOL[library])
                    self.assertAlmostEqual(gamma / egd, 1.0, delta=MT458_REL_TOL[library])

    def test_double_beta_emitters_are_left_out(self) -> None:
        # JEFF-3.1.1 gives Zr96, Mo100, Xe136 and others half-lives of ~1e27 s: with them the beta energy of U235
        # is ~0.9 MeV more than MT458
        decay_mat, energy = self.libraries["jeff311"]
        inventory = pulse_inventory(self.yields, "U235")
        beta, _ = delayed_energy(decay_mat, energy, inventory)
        beta_all, _ = delayed_energy(decay_mat, energy, inventory, max_half_life_s=math.inf)
        self.assertGreater(beta_all - beta, 0.5)

    def test_time_steps_give_the_direct_integral(self) -> None:
        # the integral of P(t) over 1e-6 s to 1e15 s, from the solver steps, is the direct solve
        decay_mat, energy = self.libraries["jeff311"]
        inventory = pulse_inventory(self.yields, "U235")
        times = np.logspace(-6, 15, 211)
        heats = pulse_decay_heat(decay_mat, energy, inventory, times)
        log_t = np.log(times)
        beta = np.trapezoid([h.time_s * h.beta for h in heats], log_t)
        gamma = np.trapezoid([h.time_s * h.gamma for h in heats], log_t)
        exact_beta, exact_gamma = delayed_energy(decay_mat, energy, inventory)
        self.assertAlmostEqual(beta / exact_beta, 1.0, delta=1e-3)
        self.assertAlmostEqual(gamma / exact_gamma, 1.0, delta=1e-3)

    def test_released_energy(self) -> None:
        # the energy released to t goes up to the delayed energy, and is the integral of the heat
        decay_mat, energy = self.libraries["jeff311"]
        inventory = pulse_inventory(self.yields, "U235")
        total = delayed_energy(decay_mat, energy, inventory)
        previous = (0.0, 0.0)
        # the half-lives go to 1e15 s (left out above), so all is released by 1e17 s
        for t in (1.0, 100.0, 1e4, 1e6, 1e9, 1e12, 1e17):
            released = released_energy(decay_mat, energy, inventory, t)
            self.assertTrue(previous[0] < released[0] <= total[0] * (1 + 1e-9))
            self.assertTrue(previous[1] < released[1] <= total[1] * (1 + 1e-9))
            previous = released
        self.assertAlmostEqual(previous[0] / total[0], 1.0, delta=1e-3)
        times = np.logspace(0, 4, 161)
        heats = pulse_decay_heat(decay_mat, energy, inventory, times)
        integral = np.trapezoid([h.time_s * h.beta for h in heats], np.log(times))
        exact = released_energy(decay_mat, energy, inventory, 1e4)[0] - released_energy(decay_mat, energy, inventory, 1.0)[0]
        self.assertAlmostEqual(integral / exact, 1.0, delta=1e-3)

    def test_fission_source_engine_is_exact(self) -> None:
        # an irradiation case of the benchmark tools (a fission source at 1 fission/s for T, then cooling) gives
        # the energy the pulse releases from t to t + T
        from nisbench.engines import _PynisData, NuclearData, _run_pynis
        from nisbench.format import parse_case
        T, cooling = 1000.0, [10.0, 300.0, 1e4]
        dts = [T] + [b - a for a, b in zip([0.0] + cooling, cooling)]
        case = parse_case({
            "format": "nis-benchmark-case", "version": 1, "benchmark": "test", "name": "source", "sample": {"atoms": {}},
            "steps": {"dt_s": dts, "flux_n_m2_s": [0.0] * len(dts), "fission_rate_per_s": [1.0] + [0.0] * len(cooling)},
            "spectrum": None, "fission_source": {"nuclide": "U235", "neutron_energy_ev": 0.0253},
            "measurements": [{"quantity": "beta_heat", "unit": "MeV/(fission/s)", "steps": [1, 2, 3],
                              "values": [1.0, 1.0, 1.0], "uncertainties": [0.0, 0.0, 0.0]}]})
        self.assertEqual(case.fissions, T)
        results = _run_pynis(case, _PynisData(NuclearData(decay=JEFF311_DECAY, fission_yields=JEFF311_NFY)))
        decay_mat, energy = self.libraries["jeff311"]
        inventory = pulse_inventory(self.yields, "U235")
        for t, result in zip(cooling, results[1:]):
            with self.subTest(time=t):
                b0 = released_energy(decay_mat, energy, inventory, t)[0]
                b1 = released_energy(decay_mat, energy, inventory, t + T)[0]
                # heat in W at 1 fission/s is the energy released from t to t+T (MeV) per s
                self.assertAlmostEqual(result["beta_heat"] / 1.602176634e-13 / (b1 - b0), 1.0, delta=1e-6)

    def test_heat_at_the_start(self) -> None:
        # at t = 0 the heat is sum(yield * lambda * E), the solver must give it at very short times
        decay_mat, energy = self.libraries["jeff311"]
        inventory = pulse_inventory(self.yields, "U235")
        heat = nis.core.compute_decay_heat(energy, nis.core.vector_from_inventory(inventory))
        start = (heat.light_particle + heat.electromagnetic) / 1.602176634e-13
        first = pulse_decay_heat(decay_mat, energy, inventory, [0.0, 1e-6])
        self.assertAlmostEqual(first[0].beta + first[0].gamma, start, delta=1e-12 * start)
        self.assertAlmostEqual((first[1].beta + first[1].gamma) / start, 1.0, delta=1e-5)

    def test_decay_heat_goes_down(self) -> None:
        # after a pulse there is no new production, so P(t) goes down (beta + gamma, from 1 s)
        decay_mat, energy = self.libraries["jeff311"]
        for parent in MT458:
            with self.subTest(parent=parent):
                heats = pulse_decay_heat(decay_mat, energy, pulse_inventory(self.yields, parent), DEFAULT_TIMES)
                values = [h.beta + h.gamma for h in heats if h.time_s >= 1.0]
                self.assertTrue(all(b < a for a, b in zip(values, values[1:])))

    def test_decay_libraries_agree(self) -> None:
        if "endfbvii1" not in self.libraries:
            self.skipTest("needs data/endfbvii1.dec")
        # with the same yields, JEFF-3.1.1 and ENDF/B-VII.1 decay data give t*P(t) within 6 %
        # from 0.1 s to 1e7 s (the largest difference, ~5 %, is at 0.1 s)
        times = [0.1, 0.3, 1, 3, 10, 30, 100, 300, 1e3, 3e3, 1e4, 3e4, 1e5, 1e6, 1e7]
        for parent in ("U235", "Pu239"):
            inventory = pulse_inventory(self.yields, parent)
            jeff = pulse_decay_heat(*self.libraries["jeff311"], inventory, times)
            endf = pulse_decay_heat(*self.libraries["endfbvii1"], inventory, times)
            for a, b in zip(jeff, endf):
                with self.subTest(parent=parent, time=a.time_s):
                    self.assertAlmostEqual((a.beta + a.gamma) / (b.beta + b.gamma), 1.0, delta=0.06)


@unittest.skipUnless(all(has_data(f) for f in (JEFF311_NFY, JEFF311_DECAY, JEFF4_REACTION)),
                     "needs data/jeff311.nfy, data/jeff311.dec and data/jeff4.reac")
class ShortIrradiationTests(unittest.TestCase):
    def test_short_irradiation_is_a_pulse(self) -> None:
        # U235 in a thermal flux for 1 ms, then cooling: the heat per fission is the heat of the pulse.
        # The fissions are the fission product atoms over the products per fission (decay keeps the atoms).
        from nisbench.engines import reaction_cache
        decay_mat, energy = _decay(JEFF311_DECAY)
        yields = nis.core.fission_yields_from_file(str(JEFF311_NFY))
        reaction_data = nis.core.reaction_data_from_file(str(reaction_cache(JEFF4_REACTION)))
        spectrum = nis.core.Spectrum([0.02, 0.03], [1e14])
        irradiation = decay_mat + nis.core.reaction_mat(reaction_data, spectrum, yields)

        model = nis.core.BasicModel()
        model.nuclides = {"U235": 1e20}
        model.dt = 1e-3
        nuclides = nis.core.compute_inventory(model, irradiation, nis.core.solver_cram48).nuclides

        pulse = pulse_inventory(yields, "U235")

        def fission_products(inventory):
            return sum(v for k, v in inventory.items() if 20 <= nis.zai.get_z(nis.zai.get_from_name(k)) <= 72)

        fissions = fission_products(nuclides) / fission_products(pulse)
        self.assertGreater(fissions, 0)

        # the U235 sample itself also decays, so look only to 1000 s, where its heat is < 0.3 %
        times = [1.0, 10.0, 100.0, 1000.0]
        expected = pulse_decay_heat(decay_mat, energy, pulse, times)
        previous = 0.0
        for t, e in zip(times, expected):
            model = nis.core.BasicModel()
            model.nuclides = nuclides
            model.dt = t - previous
            nuclides = nis.core.compute_inventory(model, decay_mat, nis.core.solver_cram48).nuclides
            previous = t
            heat = nis.core.compute_decay_heat(energy, nis.core.vector_from_inventory(nuclides))
            per_fission = (heat.light_particle + heat.electromagnetic) / 1.602176634e-13 / fissions
            with self.subTest(time=t):
                self.assertAlmostEqual(per_fission / (e.beta + e.gamma), 1.0, delta=1e-3)


__all__ = ["FissionPulseTests", "ShortIrradiationTests"]
