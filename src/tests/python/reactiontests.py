"""
Benchmarks for the neutron reaction matrix with real JEFF-4.0 data.

These are simple cases with known answers:
    - one group cross sections against reference values
    - spectrum handling (linearity and bin splitting)
    - the structure of the full matrix
    - burn-up and activation chains against the analytic Bateman solution

Reference cross sections are the 2200 m/s (0.0253 eV) values from the
Atlas of Neutron Resonances (Mughabghab) and 14.1 MeV values from evaluated
data (IRDFF-II and JEFF/ENDF). The data is group-wise, so the value in a group
is an average and we use a tolerance on the reference.
"""

from pathlib import Path
import math
import unittest

import numpy as np
import scipy.sparse

import pynis as nis


_data_dir = Path(__file__).parent.parent.parent.parent / "data"
_jeff4_reaction_file = _data_dir / "jeff4.reac"
_jeff311_decay_file = _data_dir / "jeff311.dec"

# LFS pointer files are small text files, the real file is ~1 GB
_has_reaction_data = (
    _jeff4_reaction_file.exists() and _jeff4_reaction_file.stat().st_size > 1e6
)

BARN = 1e-24
THERMAL_ENERGY = 0.0253
FUSION_ENERGY = 14.1e6
FLUX = 1e14

_ZAI_INDEX = {zai: i for i, zai in enumerate(nis.data.ALL_ZAIS())}


def index(name: str) -> int:
    return _ZAI_INDEX[nis.zai.get_from_name(name)]


def _states(ground_zai: int) -> list:
    return [ground_zai + i for i in range(10) if ground_zai + i in _ZAI_INDEX]


def one_group_spectrum(energy: float, flux: float = FLUX) -> nis.core.Spectrum:
    """A narrow group around the energy, so it falls in one data group."""
    return nis.core.Spectrum([energy * 0.9999, energy * 1.0001], [flux])


@unittest.skipUnless(_has_reaction_data, "needs data/jeff4.reac from git lfs")
class ReactionBenchmarkBase(unittest.TestCase):
    # read the tape once, it is big
    reaction_data = None

    @classmethod
    def setUpClass(cls) -> None:
        if ReactionBenchmarkBase.reaction_data is None:
            ReactionBenchmarkBase.reaction_data = nis.core.reaction_data_from_file(
                str(_jeff4_reaction_file)
            )

    def matrix(self, spectrum: nis.core.Spectrum) -> scipy.sparse.csc_matrix:
        # the matrix is sparse, a dense one is ~130 MB
        return nis.core.reaction_mat(self.reaction_data, spectrum)

    def check_xs(
        self,
        mat: scipy.sparse.csc_matrix,
        parent: str,
        daughter: str,
        expected_barns: float,
        rel_tol: float,
        flux: float = FLUX,
    ) -> None:
        # add the isomers of the daughter, the reference is for all states
        xs = sum(
            mat[index(parent), _ZAI_INDEX[zai]]
            for zai in _states(nis.zai.get_from_name(daughter))
        ) / flux / BARN
        self.assertTrue(
            math.isclose(xs, expected_barns, rel_tol=rel_tol),
            f"{parent} -> {daughter}: {xs:.5g} b, expected {expected_barns:.5g} b",
        )


class ThermalCrossSectionTests(ReactionBenchmarkBase):
    def test_thermal_reactions(self) -> None:
        mat = self.matrix(one_group_spectrum(THERMAL_ENERGY))
        cases = [
            # parent, daughter, reference in barns
            ("Co59", "Co60", 37.18),  # (n,g)
            ("Au197", "Au198", 98.65),  # (n,g)
            ("Mn55", "Mn56", 13.36),  # (n,g)
            ("Na23", "Na24", 0.530),  # (n,g)
            ("U238", "U239", 2.680),  # (n,g)
            ("H1", "H2", 0.3326),  # (n,g)
            ("B10", "Li7", 3842.0),  # (n,a)
            ("Li6", "He4", 940.0),  # (n,t) gives He4 + H3
            ("Li6", "H3", 940.0),
            ("He3", "H3", 5333.0),  # (n,p) gives H3 + H1
        ]
        for parent, daughter, expected in cases:
            with self.subTest(parent=parent, daughter=daughter):
                self.check_xs(mat, parent, daughter, expected, rel_tol=0.03)

    def test_thermal_fission_to_sink(self) -> None:
        mat = self.matrix(one_group_spectrum(THERMAL_ENERGY))
        u235 = index("U235")
        xs = mat[u235, nis.data.SINK_INDEX] / FLUX / BARN
        self.assertTrue(math.isclose(xs, 585.1, rel_tol=0.03), f"{xs} b")

    def test_loss_is_sum_of_channels(self) -> None:
        # for (n,g) only targets the loss rate is the rate to the daughter states
        mat = self.matrix(one_group_spectrum(THERMAL_ENERGY))
        co59 = index("Co59")
        co60, co60m = index("Co60"), index("Co60m")
        self.assertTrue(
            math.isclose(-mat[co59, co59], mat[co59, co60] + mat[co59, co60m], rel_tol=1e-6)
        )

    def test_isomer_production(self) -> None:
        # Ni62 (n,p) makes both Co62 and Co62m, from MF=10
        mat = self.matrix(one_group_spectrum(FUSION_ENERGY))
        ni62 = index("Ni62")
        self.assertGreater(mat[ni62, index("Co62")], 0.0)
        self.assertGreater(mat[ni62, index("Co62m")], 0.0)


class FusionCrossSectionTests(ReactionBenchmarkBase):
    def test_14mev_reactions(self) -> None:
        mat = self.matrix(one_group_spectrum(FUSION_ENERGY))
        cases = [
            ("Fe56", "Mn56", 0.110),  # (n,p)
            ("Al27", "Na24", 0.122),  # (n,a)
            ("Ni58", "Co58", 0.355),  # (n,p)
            ("Zr90", "Zr89", 0.62),  # (n,2n)
            ("Nb93", "Nb92", 1.24),  # (n,2n)
        ]
        for parent, daughter, expected in cases:
            with self.subTest(parent=parent, daughter=daughter):
                self.check_xs(mat, parent, daughter, expected, rel_tol=0.10)

    def test_14mev_fission_to_sink(self) -> None:
        mat = self.matrix(one_group_spectrum(FUSION_ENERGY))
        u238 = index("U238")
        xs = mat[u238, nis.data.SINK_INDEX] / FLUX / BARN
        self.assertTrue(math.isclose(xs, 1.15, rel_tol=0.10), f"{xs} b")

    def test_threshold_reactions_off_at_thermal(self) -> None:
        mat = self.matrix(one_group_spectrum(THERMAL_ENERGY))
        self.assertEqual(mat[index("Fe56"), index("Mn56")], 0.0)
        self.assertEqual(mat[index("Nb93"), index("Nb92")], 0.0)


def assert_sparse_allclose(a: scipy.sparse.csc_matrix, b: scipy.sparse.csc_matrix, rtol: float) -> None:
    # same non-zero entries, and the same values within rtol
    a, b = a.tocsc(), b.tocsc()
    a.sort_indices()
    b.sort_indices()
    np.testing.assert_array_equal(a.indptr, b.indptr)
    np.testing.assert_array_equal(a.indices, b.indices)
    np.testing.assert_allclose(a.data, b.data, rtol=rtol, atol=0)


class SpectrumTests(ReactionBenchmarkBase):
    def test_linear_in_flux(self) -> None:
        spectrum = nis.core.Spectrum([1.0, 1e3, 1e6, 2e7], [1e12, 2e12, 3e12])
        double = nis.core.Spectrum(spectrum.x, [2 * v for v in spectrum.y])
        assert_sparse_allclose(2 * self.matrix(spectrum), self.matrix(double), rtol=1e-12)

    def test_split_bins(self) -> None:
        # the flux is flat in lethargy within a bin, so splitting a bin
        # into two with equal lethargy widths gives the same rates
        one_bin = nis.core.Spectrum([1e3, 1e5], [2e12])
        two_bins = nis.core.Spectrum([1e3, 1e4, 1e5], [1e12, 1e12])
        assert_sparse_allclose(self.matrix(one_bin), self.matrix(two_bins), rtol=1e-10)

    def test_no_flux_gives_zero_matrix(self) -> None:
        # flux outside the data range (above 1 GeV)
        mat = self.matrix(nis.core.Spectrum([2e9, 3e9], [1e14]))
        self.assertEqual(mat.count_nonzero(), 0)


class MatrixStructureTests(ReactionBenchmarkBase):
    def test_fusion_spectrum_structure(self) -> None:
        # a rough fusion spectrum, thermal to 14 MeV peak
        spectrum = nis.core.Spectrum(
            [1e-5, 1.0, 1e3, 1e5, 1e6, 1e7, 1.3e7, 1.5e7, 2e7],
            [1e12, 2e12, 3e12, 3e12, 4e12, 1e12, 5e13, 1e10],
        )
        mat = self.matrix(spectrum).toarray()
        self.assertTrue(np.all(np.isfinite(mat)))

        diag = np.diag(mat).copy()
        off = mat - np.diag(diag)
        self.assertTrue(np.all(diag <= 0.0))
        self.assertTrue(np.all(off >= 0.0))
        # the sink row has no reactions
        self.assertFalse(np.any(mat[nis.data.SINK_INDEX, :]))
        # each reaction gives at least one product, more if a light
        # particle is emitted, so the outgoing rate is never below the loss
        out = off.sum(axis=1)
        self.assertTrue(np.all(out >= -diag * (1 - 1e-12)))
        # most nuclides in the library must have reactions
        self.assertGreater(np.count_nonzero(diag), 500)


class BurnupTests(ReactionBenchmarkBase):
    def inventory(self, mat: scipy.sparse.csc_matrix, nuclides: dict, dt: float) -> dict:
        model = nis.core.BasicModel()
        model.nuclides = nuclides
        model.dt = dt
        result = nis.core.compute_inventory(model, mat, nis.core.solver_cram48)
        return result.nuclides

    def test_au197_capture_chain(self) -> None:
        # Au197 -(n,g)-> Au198 -(n,g)-> Au199, no decay
        # compare with the analytic two step Bateman solution
        # part of the (n,g) goes to Au198m, so use the fraction to the ground state
        flux = 1e15
        mat = self.matrix(one_group_spectrum(THERMAL_ENERGY, flux))
        au197, au198 = index("Au197"), index("Au198")
        a = -mat[au197, au197]
        b = -mat[au198, au198]
        to_ground = mat[au197, au198]

        n0 = 1e20
        dt = 1.0 / a  # ~63% of Au197 burnt
        nuclides = self.inventory(mat, {"Au197": n0}, dt)

        expected_197 = n0 * math.exp(-a * dt)
        expected_198 = n0 * to_ground / (b - a) * (math.exp(-a * dt) - math.exp(-b * dt))
        self.assertTrue(math.isclose(nuclides["Au197"], expected_197, rel_tol=1e-8))
        self.assertTrue(math.isclose(nuclides["Au198"], expected_198, rel_tol=1e-8))

    def test_atoms_conserved_for_capture(self) -> None:
        # Co59 at thermal is (n,g) only, so the number of atoms stays the same
        flux = 1e15
        mat = self.matrix(one_group_spectrum(THERMAL_ENERGY, flux))
        n0 = 1e20
        nuclides = self.inventory(mat, {"Co59": n0}, 1e8)
        self.assertTrue(math.isclose(sum(nuclides.values()), n0, rel_tol=1e-6))

    def test_co60_activation(self) -> None:
        # classic activation benchmark: Co59 in a thermal flux, with Co60
        # decay added by hand (half life 5.2711 y) as no decay data is used here.
        # Low burn-up (~6e-4) so the activity is A = N0 sigma phi (1 - exp(-lambda t))
        flux = 1e11
        # LIL format, so we can add the decay entries
        mat = self.matrix(one_group_spectrum(THERMAL_ENERGY, flux)).tolil()
        co59, co60, ni60 = index("Co59"), index("Co60"), index("Ni60")
        lam = math.log(2) / (5.2711 * 365.25 * 24 * 3600)
        mat[co60, co60] -= lam
        mat[co60, ni60] += lam

        n0 = 1e24
        sigma_phi = -mat[co59, co59]
        dt = 5.2711 * 365.25 * 24 * 3600  # one half life
        nuclides = self.inventory(mat, {"Co59": n0}, dt)

        # full Bateman with burn-up of both nuclides
        a = sigma_phi
        b = -mat[co60, co60]
        expected_co60 = n0 * a / (b - a) * (math.exp(-a * dt) - math.exp(-b * dt))
        self.assertTrue(math.isclose(nuclides["Co60"], expected_co60, rel_tol=1e-8))

        # the simple textbook formula is within the burn-up effect
        activity = lam * nuclides["Co60"]
        saturation = n0 * sigma_phi * (1 - math.exp(-lam * dt))
        self.assertTrue(math.isclose(activity, saturation, rel_tol=1e-2))
        self.assertTrue(math.isclose(activity / (n0 * sigma_phi), 0.5, rel_tol=1e-2))


class ReachableTests(ReactionBenchmarkBase):
    """Solving only the reachable nuclides must give the same as the full matrix."""

    def check(self, mat, inventory: dict, dt: float) -> None:
        n0 = nis.core.vector_from_inventory(inventory)
        full = np.asarray(nis.core.solver_cram48(mat, n0, dt))
        reduced = np.asarray(nis.core.solve_reachable(mat, n0, dt, nis.core.solver_cram48))
        total = n0.sum()

        reachable = np.zeros(len(n0), dtype=bool)
        reachable[nis.core.reachable_indices(mat, n0)] = True
        # fewer nuclides, else this does not help
        self.assertLess(reachable.sum(), len(n0))
        # the full solve only has rounding noise where it can't reach
        self.assertTrue(np.all(reduced[~reachable] == 0.0))
        self.assertLess(np.abs(full[~reachable]).max(), 1e-12 * total)
        # the same everywhere, to rounding
        self.assertLess(np.abs(reduced - full).max(), 1e-12 * total)
        big = np.abs(full) > 1e-9 * total
        np.testing.assert_allclose(reduced[big], full[big], rtol=1e-8)

    def test_irradiation_and_decay(self) -> None:
        decay = nis.core.decay_mat_from_file(str(_jeff311_decay_file))
        spectrum = nis.core.Spectrum(
            [1e-5, 1.0, 1e3, 1e5, 1e6, 1e7, 1.3e7, 1.5e7, 2e7],
            [1e12, 2e12, 3e12, 3e12, 4e12, 1e12, 5e13, 1e10],
        )
        irradiation = decay + self.matrix(spectrum)
        inventory = {"Fe56": 1e24, "Cr52": 2e23, "Ni58": 1e23, "Co59": 1e20, "U235": 1e18}
        for dt in (60.0, 3600.0 * 24, 3600.0 * 24 * 365):
            with self.subTest(dt=dt):
                self.check(irradiation, inventory, dt)
                self.check(decay, inventory, dt)

        # the inventory after irradiation, then decay, as the FNS benchmark does
        model = nis.core.BasicModel()
        model.nuclides = inventory
        model.dt = 3600.0
        shutdown = nis.core.compute_inventory(model, irradiation, nis.core.solver_cram48).nuclides
        self.check(decay, shutdown, 3600.0)


__all__ = [
    "ThermalCrossSectionTests",
    "FusionCrossSectionTests",
    "SpectrumTests",
    "MatrixStructureTests",
    "BurnupTests",
    "ReachableTests",
]
