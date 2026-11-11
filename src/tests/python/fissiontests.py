"""
Fission yield tests with the JEFF-3.1.1 yields (data/jeff311.nfy and data/jeff311.sfy, see the README),
the JEFF-3.1.1 decay data and the JEFF-4.0 reaction data. The tests skip if a file is not there.
"""

import math
import unittest

import pynis as nis

from .testdata import JEFF311_DECAY, JEFF311_NFY, JEFF311_SFY, JEFF4_REACTION, has_data

from nisbench.engines import reaction_cache

U235 = 922350
CF252 = 982520
SECONDS_PER_YEAR = 365.25 * 86400


def _sum_atoms(nuclides: dict, names) -> float:
    return sum(nuclides.get(n, 0.0) for n in names)


@unittest.skipUnless(has_data(JEFF311_NFY), "needs data/jeff311.nfy")
class FissionYieldDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.yields = nis.core.fission_yields_from_file(str(JEFF311_NFY))

    def test_parents(self) -> None:
        self.assertEqual(len(self.yields), 19)
        for zai in (922350, 922380, 942390, 942410, 952421):
            self.assertIn(zai, self.yields)
        # thermal, 400 keV and 14 MeV
        self.assertEqual(len(self.yields[U235].energies), 3)

    def test_two_fragments_per_fission(self) -> None:
        # the independent yields add up to 2 products per fission
        for zai in self.yields.keys():
            for products in self.yields[zai].yields:
                self.assertTrue(1.95 < sum(products.values()) < 2.05, f"{zai}: {sum(products.values())}")


@unittest.skipUnless(all(has_data(f) for f in (JEFF311_NFY, JEFF311_DECAY, JEFF4_REACTION)),
                     "needs data/jeff311.nfy, data/jeff311.dec and data/jeff4.reac")
class ThermalFissionTests(unittest.TestCase):
    """1e22 atoms of U235 in a thermal flux for 1 h, then 1 day of cooling."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.decay_mat = nis.core.decay_mat_from_file(str(JEFF311_DECAY))
        cls.reaction_data = nis.core.reaction_data_from_file(str(reaction_cache(JEFF4_REACTION)))
        cls.yields = nis.core.fission_yields_from_file(str(JEFF311_NFY))
        cls.spectrum = nis.core.Spectrum([0.02, 0.03], [1e14])

    def run_u235(self, yields):
        mat = self.decay_mat + nis.core.reaction_mat(self.reaction_data, self.spectrum, yields)
        model = nis.core.BasicModel()
        model.nuclides = {"U235": 1e22}
        model.dt = 3600.0
        irradiated = nis.core.compute_inventory(model, mat, nis.core.solver_cram48).nuclides
        model.nuclides = irradiated
        model.dt = 86400.0
        cooled = nis.core.compute_inventory(model, self.decay_mat, nis.core.solver_cram48).nuclides
        return irradiated, cooled

    def test_long_lived_products_give_the_cumulative_yield(self) -> None:
        irradiated, cooled = self.run_u235(self.yields)
        # the fissions: U235 burnt minus the captures (U236)
        fissions = 1e22 - irradiated["U235"] - irradiated["U236"]
        # thermal cumulative yields of U235 (JEFF-3.1.1): Cs137 6.19 %, Sr90 5.78 %, Zr95 6.50 %
        for name, cumulative in (("Cs137", 0.0619), ("Sr90", 0.0578)):
            self.assertTrue(math.isclose(cooled[name] / fissions, cumulative, rel_tol=0.05),
                            f"{name}: {cooled[name] / fissions:.4f} per fission")
        # the mass 95 chain is in Zr95 or Nb95 after 1 day
        self.assertTrue(math.isclose(_sum_atoms(cooled, ["Zr95", "Nb95", "Nb95m"]) / fissions, 0.065, rel_tol=0.05))

    def test_fission_ratio_to_capture(self) -> None:
        irradiated, _ = self.run_u235(self.yields)
        fissions = 1e22 - irradiated["U235"] - irradiated["U236"]
        # thermal U235: 585 b fission, 99 b capture
        self.assertTrue(math.isclose(fissions / irradiated["U236"], 585 / 99, rel_tol=0.05))

    def test_no_yields_goes_to_the_sink(self) -> None:
        irradiated, _ = self.run_u235(None)
        # less than 1 atom (CRAM can leave rounding noise)
        self.assertLess(irradiated.get("Cs137", 0.0), 1.0)
        self.assertGreater(irradiated["U236"], 0.0)


@unittest.skipUnless(has_data(JEFF311_SFY) and has_data(JEFF311_DECAY), "needs data/jeff311.sfy and data/jeff311.dec")
class SpontaneousFissionTests(unittest.TestCase):
    def test_cf252(self) -> None:
        sf_yields = nis.core.fission_yields_from_file(str(JEFF311_SFY))
        self.assertIn(CF252, sf_yields)
        mat = nis.core.decay_mat_from_file(str(JEFF311_DECAY), sf_yields)
        model = nis.core.BasicModel()
        model.nuclides = {"Cf252": 1e20}
        model.dt = 0.1 * SECONDS_PER_YEAR
        cf = nis.core.compute_inventory(model, mat, nis.core.solver_cram48).nuclides

        decays = 1e20 - cf["Cf252"]
        # 3.09 % of Cf252 decays are spontaneous fission, Cs137 has a cumulative yield of about 6 %
        sf = decays * 0.0309
        self.assertTrue(0.04 < cf["Cs137"] / sf < 0.08, cf["Cs137"] / sf)

        # without the yields, no fission products
        no_sf = nis.core.compute_inventory(model, nis.core.decay_mat_from_file(str(JEFF311_DECAY)),
                                           nis.core.solver_cram48).nuclides
        self.assertLess(no_sf.get("Cs137", 0.0), 1.0)


__all__ = ["FissionYieldDataTests", "ThermalFissionTests", "SpontaneousFissionTests"]
