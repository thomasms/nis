"""
Checks of the decay matrix and the decay energy data made from each decay library. They find errors in
the data and in how NIS reads it (a wrong daughter, a lost branch, a wrong isomer).

Needs the decay data from git lfs (data/jeff311.dec and data/endfbvii1.dec), else the tests skip.
"""

import math
import unittest

import numpy as np
import scipy.sparse as sp

import pynis as nis

from .testdata import ENDFBVII_DECAY, JEFF311_DECAY, has_data

YEAR = 365.25 * 86400

# half-lives (s) that the libraries must give to 0.5 % (from the evaluated data, e.g. NUBASE/ENSDF)
KNOWN_HALF_LIVES = {
    "H3": 12.32 * YEAR,
    "C14": 5700 * YEAR,
    "Na24": 14.997 * 3600,
    "Co60": 5.2714 * YEAR,
    "Co60m": 10.467 * 60,
    "Ge75": 82.78 * 60,
    "Ge75m": 47.7,
    "As75m": 17.62e-3,
    "Sr90": 28.79 * YEAR,
    "Tc99m": 6.0072 * 3600,
    "I131": 8.0252 * 86400,
    "Cs137": 30.08 * YEAR,
    "U238": 4.468e9 * YEAR,
}

# isomers and the fraction that must go to the ground state (isomeric transition)
ISOMERIC_TRANSITIONS = {
    "Co60m": ("Co60", 0.99),
    "Ge75m": ("Ge75", 0.99),
    "As75m": ("As75", 0.99),
    "Tc99m": ("Tc99", 0.99),
    "Am242m": ("Am242", 0.99),
}

# the most nucleons a decay can lose to particles that NIS does not track as a nuclide (neutrons)
MAX_LOST_NUCLEONS_PER_DECAY = 4.0


class DecayDataChecks:
    """The checks for one decay library. A subclass sets DATA, and the nuclides with known errors in the data."""

    DATA = None
    # nuclides whose branching ratios add up to more than 1, or to less than 0.99, in the evaluation
    BRANCHING_SUM_NOT_1 = set()
    matrix = None

    @classmethod
    def setUpClass(cls) -> None:
        cls.matrix = nis.core.decay_mat_from_file(str(cls.DATA)).tocsr()
        cls.energy = nis.core.decay_energy_data_from_file(str(cls.DATA))
        cls.zais = np.array(nis.data.ALL_ZAIS())
        cls.sink = nis.data.SINK_INDEX
        cls.lam = -cls.matrix.diagonal()
        cls.off = (cls.matrix - sp.diags(cls.matrix.diagonal())).tocsr()
        cls.off.eliminate_zeros()
        cls.unstable = np.nonzero(cls.lam > 0)[0]

    def name(self, index: int) -> str:
        return "sink" if index == self.sink else nis.zai.get_name_from(int(self.zais[index]))

    def index(self, name: str) -> int:
        return int(np.argmax(nis.core.vector_from_inventory({name: 1.0})))

    def branches(self, name: str) -> dict:
        i = self.index(name)
        row = self.off[i]
        return {self.name(j): v / self.lam[i] for j, v in zip(row.indices, row.data)}

    def test_size(self) -> None:
        self.assertEqual(self.matrix.shape, (nis.data.ZAI_SIZE, nis.data.ZAI_SIZE))
        # the libraries have ~3600 radionuclides
        self.assertGreater(len(self.unstable), 3000)

    def test_signs(self) -> None:
        # a nuclide only loses atoms by its own decay, and only gives atoms to its daughters
        self.assertTrue(np.all(self.lam >= 0), "a decay constant below 0")
        self.assertTrue(np.all(self.off.data > 0), "a production rate below 0")
        self.assertTrue(np.all(np.isfinite(self.matrix.data)))

    def test_stable_nuclides_and_the_sink_do_not_decay(self) -> None:
        self.assertEqual(self.matrix[self.sink].nnz, 0)
        for name in ("H1", "He4", "C12", "O16", "Fe56", "As75", "Pb208"):
            with self.subTest(nuclide=name):
                i = self.index(name)
                self.assertEqual(self.lam[i], 0.0)
                self.assertEqual(self.off[i].nnz, 0)

    def test_branching_ratios(self) -> None:
        light = {self.index(n) for n in ("H1", "H2", "H3", "He3", "He4")} | {self.sink}
        for i in self.unstable:
            row = self.off[i]
            with self.subTest(nuclide=self.name(i)):
                # one branch to a heavy daughter is at most all of the decays
                heavy = [v for j, v in zip(row.indices, row.data) if j not in light]
                self.assertLessEqual(max(heavy, default=0.0), self.lam[i] * (1 + 1e-9))
                # each decay makes at least one atom (the branching ratios add up to 1) and a few at most.
                # Some evaluations add up to a little less than 1 (In129m in JEFF-3.1.1: 0.997).
                atoms_per_decay = row.data.sum() / self.lam[i]
                self.assertGreaterEqual(atoms_per_decay, 0.995)
                self.assertLessEqual(atoms_per_decay, 10.0)

    def test_mass_number_is_not_made(self) -> None:
        # a decay cannot make nucleons. It can lose a few to particles that are not tracked (neutrons).
        mass = np.append(self.zais // 10 % 1000, 0)
        for i in self.unstable:
            row = self.off[i]
            if self.sink in row.indices or self.name(i) in self.BRANCHING_SUM_NOT_1:
                # spontaneous fission or a daughter that is not tracked: the mass is not known
                continue
            with self.subTest(nuclide=self.name(i)):
                lost = mass[i] - (mass[row.indices] * row.data).sum() / self.lam[i]
                # some branching ratios add up to a little more than 1 (up to ~1e-5 in JEFF-3.1.1)
                self.assertGreaterEqual(lost, -1e-4 * mass[i])
                self.assertLessEqual(lost, MAX_LOST_NUCLEONS_PER_DECAY)

    def test_most_decays_conserve_mass_number(self) -> None:
        # beta, alpha, IT and EC decays conserve the mass number exactly (He4 is tracked)
        mass = np.append(self.zais // 10 % 1000, 0)
        conserved = 0
        for i in self.unstable:
            row = self.off[i]
            lost = mass[i] - (mass[row.indices] * row.data).sum() / self.lam[i]
            if abs(lost) <= 1e-6 * mass[i]:
                conserved += 1
        self.assertGreater(conserved / len(self.unstable), 0.85)

    def test_known_half_lives(self) -> None:
        for name, half_life in KNOWN_HALF_LIVES.items():
            with self.subTest(nuclide=name):
                lam = self.lam[self.index(name)]
                self.assertGreater(lam, 0.0)
                self.assertAlmostEqual(math.log(2) / lam / half_life, 1.0, delta=0.005)

    def test_isomeric_transitions(self) -> None:
        for isomer, (ground, fraction) in ISOMERIC_TRANSITIONS.items():
            with self.subTest(isomer=isomer):
                self.assertGreaterEqual(self.branches(isomer).get(ground, 0.0), fraction)

    def test_decay_constants_are_the_same_in_the_energy_data(self) -> None:
        np.testing.assert_allclose(np.array(self.energy.decay_constants), self.lam, rtol=1e-12, atol=0)

    def test_decay_energies(self) -> None:
        stable = self.lam == 0
        total = np.zeros_like(self.lam)
        for part in ("light_particle", "electromagnetic", "heavy_particle"):
            values = np.array(getattr(self.energy, part))
            with self.subTest(part=part):
                self.assertTrue(np.all(np.isfinite(values)))
                self.assertTrue(np.all(values >= 0))
                self.assertTrue(np.all(values[stable] == 0))
                # no decay gives more than 1 GeV (in eV)
                self.assertTrue(np.all(values <= 1e9))
            total += values
        # most radionuclides have a decay energy
        self.assertGreater(np.count_nonzero(total[self.unstable]) / len(self.unstable), 0.9)

    def test_alpha_decay_makes_he4(self) -> None:
        # U238 -> Th234 + He4
        branches = self.branches("U238")
        self.assertAlmostEqual(branches["Th234"], 1.0, delta=1e-5)
        self.assertAlmostEqual(branches["He4"], 1.0, delta=1e-5)

    def test_proton_emission(self) -> None:
        # Li5 -> He4 + p (JEFF-3.1.1 has it as proton emission, ENDF/B-VII.1 as alpha decay to H1).
        # The proton is an H1 atom, and the daughter has one nucleon less.
        branches = self.branches("Li5")
        self.assertAlmostEqual(branches.get("He4", 0.0), 1.0, delta=1e-6)
        self.assertAlmostEqual(branches.get("H1", 0.0), 1.0, delta=1e-6)
        self.assertNotIn("He5", branches)


@unittest.skipUnless(has_data(JEFF311_DECAY), "needs data/jeff311.dec from git lfs")
class Jeff311DecayDataTests(DecayDataChecks, unittest.TestCase):
    DATA = JEFF311_DECAY
    # errors in the JEFF-3.1.1 evaluations: the sums of the branching ratios are 1.01 (Ge81m, Er152),
    # 1.002 (Au183), 1.0002 (Hg186), 0.98 (Rn221) and 1.0015 (Pa228)
    BRANCHING_SUM_NOT_1 = {"Ge81m", "Er152", "Au183", "Hg186", "Rn221", "Pa228"}

    def test_known_data_errors(self) -> None:
        # the nuclides in BRANCHING_SUM_NOT_1 still have the error, if not, take them out of the list
        mass = np.append(self.zais // 10 % 1000, 0)
        for name in self.BRANCHING_SUM_NOT_1:
            with self.subTest(nuclide=name):
                i = self.index(name)
                row = self.off[i]
                lost = mass[i] - (mass[row.indices] * row.data).sum() / self.lam[i]
                self.assertTrue(lost < -1e-4 * mass[i] or lost > MAX_LOST_NUCLEONS_PER_DECAY)

    def test_ge75m_feeds_as75m(self) -> None:
        # the FNS As case: in JEFF-3.1.1 the As75m comes only from the beta branch of Ge75m
        self.assertAlmostEqual(self.branches("Ge75m")["As75m"], 3e-4, delta=1e-6)
        self.assertNotIn("As75m", self.branches("Ge75"))


@unittest.skipUnless(has_data(ENDFBVII_DECAY), "needs data/endfbvii1.dec from git lfs")
class EndfbVII1DecayDataTests(DecayDataChecks, unittest.TestCase):
    DATA = ENDFBVII_DECAY


@unittest.skipUnless(has_data(JEFF311_DECAY) and has_data(ENDFBVII_DECAY),
                     "needs data/jeff311.dec and data/endfbvii1.dec from git lfs")
class DecayLibrariesAgreeTests(unittest.TestCase):
    def test_half_lives_of_long_lived_nuclides(self) -> None:
        # the half-lives of nuclides that live longer than 1 hour are well measured, so the two libraries
        # must give about the same values
        jeff = -nis.core.decay_mat_from_file(str(JEFF311_DECAY)).diagonal()
        endf = -nis.core.decay_mat_from_file(str(ENDFBVII_DECAY)).diagonal()
        long_lived = (jeff > 0) & (endf > 0) & (jeff < math.log(2) / 3600) & (endf < math.log(2) / 3600)
        ratio = jeff[long_lived] / endf[long_lived]
        self.assertGreater(long_lived.sum(), 500)
        self.assertGreater(np.mean(np.abs(ratio - 1) < 0.05), 0.95)


__all__ = ["Jeff311DecayDataTests", "EndfbVII1DecayDataTests", "DecayLibrariesAgreeTests"]
