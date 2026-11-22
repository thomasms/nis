"""
Tests for the nuclear data API: decay data, cross sections and fission yields made and changed in python.

These need no data files. The data is made by hand, then the matrices are made from it.
"""

import math
import unittest

import numpy as np

import pynis as nis

core = nis.core

CO60 = 270600
NI60 = 280600
PO210 = 842100
PB206 = 822060
HE4 = 20040
FE56 = 260560
FE57 = 260570
MN56 = 250560
U235 = 922350
BA140 = 561400
MO100 = 421000

BARN_TO_CM2 = 1e-24
CO60_HALF_LIFE = 1.663415e8


def index(zai):
    return nis.data.ALL_ZAIS().index(zai)


def co60_decay_data():
    decay_data = core.DecayData()
    decay_data[CO60] = core.DecayNuclide(half_life=CO60_HALF_LIFE, modes=[core.DecayMode(rtyp=1.0, branching_ratio=1.0)],
                                         light_particle=9.577e4, electromagnetic=2.5036e6)
    decay_data[NI60] = core.DecayNuclide()
    return decay_data


class DecayDataTest(unittest.TestCase):

    def test_stable_by_default(self):
        nuclide = core.DecayNuclide()
        self.assertTrue(nuclide.is_stable)
        self.assertEqual(nuclide.decay_constant, 0.0)
        self.assertTrue(math.isinf(nuclide.half_life))
        self.assertEqual(len(nuclide.modes), 0)

    def test_decay_constant_sets_half_life(self):
        nuclide = core.DecayNuclide(half_life=10.0)
        self.assertAlmostEqual(nuclide.decay_constant, math.log(2)/10.0, delta=1e-15)
        nuclide.decay_constant = 0.5
        self.assertAlmostEqual(nuclide.half_life, math.log(2)/0.5, delta=1e-12)
        nuclide.decay_constant = 0.0
        self.assertTrue(nuclide.is_stable)

    def test_dict_interface(self):
        decay_data = co60_decay_data()
        self.assertEqual(len(decay_data), 2)
        self.assertIn(CO60, decay_data)
        self.assertEqual(sorted(decay_data.keys()), [CO60, NI60])
        del decay_data[NI60]
        self.assertNotIn(NI60, decay_data)

    def test_decay_matrix(self):
        mat = core.decay_mat(co60_decay_data()).toarray()
        lam = math.log(2)/CO60_HALF_LIFE
        self.assertEqual(mat.shape, (nis.data.ZAI_SIZE, nis.data.ZAI_SIZE))
        self.assertAlmostEqual(mat[index(CO60), index(CO60)], -lam, delta=lam*1e-14)
        self.assertAlmostEqual(mat[index(CO60), index(NI60)], lam, delta=lam*1e-14)
        self.assertEqual(np.count_nonzero(mat), 2)

    def test_change_in_place(self):
        """decay_data[zai] is the nuclide in the data, not a copy, so changes to it change the matrix"""
        decay_data = co60_decay_data()
        decay_data[CO60].half_life = CO60_HALF_LIFE/2.0
        decay_data[CO60].modes[0].branching_ratio = 0.75
        decay_data[CO60].modes.append(core.DecayMode(rtyp=10, branching_ratio=0.25))
        self.assertAlmostEqual(decay_data[CO60].total_branching_ratio, 1.0, delta=1e-15)

        mat = core.decay_mat(decay_data).toarray()
        lam = 2.0*math.log(2)/CO60_HALF_LIFE
        self.assertAlmostEqual(mat[index(CO60), index(CO60)], -lam, delta=lam*1e-14)
        self.assertAlmostEqual(mat[index(CO60), index(NI60)], 0.75*lam, delta=lam*1e-14)
        self.assertAlmostEqual(mat[index(CO60), nis.data.SINK_INDEX], 0.25*lam, delta=lam*1e-14)

    def test_set_modes_from_a_list(self):
        decay_data = core.DecayData()
        decay_data[PO210] = core.DecayNuclide(half_life=1.195569e7)
        decay_data[PO210].modes = [core.DecayMode(4.0, 1.0)]
        mat = core.decay_mat(decay_data).toarray()
        lam = math.log(2)/1.195569e7
        self.assertAlmostEqual(mat[index(PO210), index(PB206)], lam, delta=lam*1e-14)
        self.assertAlmostEqual(mat[index(PO210), index(HE4)], lam, delta=lam*1e-14)

    def test_decay_energy_data(self):
        energy = core.decay_energy_data(co60_decay_data())
        self.assertAlmostEqual(energy.decay_constants[index(CO60)], math.log(2)/CO60_HALF_LIFE, delta=1e-20)
        self.assertEqual(energy.light_particle[index(CO60)], 9.577e4)
        self.assertEqual(energy.electromagnetic[index(CO60)], 2.5036e6)
        self.assertEqual(energy.decay_constants[index(NI60)], 0.0)

    def test_solve(self):
        """Co60 after one half-life, with the matrix from the data made by hand"""
        model = core.BasicModel()
        model.nuclides = {"Co60": 1e20}
        model.dt = CO60_HALF_LIFE
        result = core.compute_inventory(model, core.decay_mat(co60_decay_data()), core.solver_cram48)
        self.assertAlmostEqual(result.nuclides["Co60"]/1e20, 0.5, delta=1e-10)
        self.assertAlmostEqual(result.nuclides["Ni60"]/1e20, 0.5, delta=1e-10)

    def test_products(self):
        self.assertEqual(core.decay_products(PO210, 4.0), (PB206, [HE4]))
        self.assertEqual(core.expand_rtypes(1.5), [1, 5])


class ReactionDataTest(unittest.TestCase):

    flux = core.Spectrum([1.0, 10.0, 100.0], [1e12, 4e12])

    def fe56_data(self):
        reaction_data = core.ReactionData()
        reaction_data.append(core.CrossSection(parent_zai=FE56, mt_number=102,
                                               energies=[1.0, 10.0, 100.0], values=[2.0, 3.0, 3.0]))
        reaction_data.append(core.CrossSection(FE56, 103, [10.0, 100.0], [0.5, 0.5]))
        return reaction_data

    def test_cross_section(self):
        xs = core.CrossSection(FE56, 102, [1.0, 10.0], [1.0, 1.0])
        self.assertEqual(xs.final_state, core.ALL_STATES)
        self.assertEqual(xs.product_za, 0)
        xs.values = np.array([2.0, 2.0])
        self.assertEqual(xs.values, [2.0, 2.0])

    def test_reaction_matrix(self):
        mat = core.reaction_mat(self.fe56_data(), self.flux).toarray()
        ng = (2.0*1e12 + 3.0*4e12)*BARN_TO_CM2
        np_rate = 0.5*4e12*BARN_TO_CM2
        self.assertAlmostEqual(mat[index(FE56), index(FE57)], ng, delta=ng*1e-12)
        self.assertAlmostEqual(mat[index(FE56), index(MN56)], np_rate, delta=np_rate*1e-12)
        self.assertAlmostEqual(mat[index(FE56), index(FE56)], -(ng + np_rate), delta=ng*1e-12)

    def test_change_in_place(self):
        """reaction_data[i] is the cross section in the data, not a copy"""
        reaction_data = self.fe56_data()
        before = core.reaction_mat(reaction_data, self.flux).toarray()
        reaction_data[0].values = [2.0*v for v in reaction_data[0].values]
        for xs in reaction_data.find(FE56, 103):
            xs.values = np.array(xs.values)*3.0
        after = core.reaction_mat(reaction_data, self.flux).toarray()
        self.assertAlmostEqual(after[index(FE56), index(FE57)]/before[index(FE56), index(FE57)], 2.0, delta=1e-12)
        self.assertAlmostEqual(after[index(FE56), index(MN56)]/before[index(FE56), index(MN56)], 3.0, delta=1e-12)

    def test_find_and_remove(self):
        reaction_data = self.fe56_data()
        self.assertEqual(len(reaction_data.find(FE56)), 2)
        self.assertEqual([xs.mt_number for xs in reaction_data.find(FE56, 102)], [102])
        self.assertEqual(reaction_data.remove_reactions(FE56, 102), 1)
        self.assertEqual(len(reaction_data), 1)
        mat = core.reaction_mat(reaction_data, self.flux).toarray()
        self.assertEqual(mat[index(FE56), index(FE57)], 0.0)
        self.assertEqual(reaction_data.remove_reactions(FE56), 1)
        self.assertEqual(len(reaction_data), 0)

    def test_list_interface(self):
        reaction_data = self.fe56_data()
        self.assertEqual([xs.mt_number for xs in reaction_data], [102, 103])
        del reaction_data[0]
        self.assertEqual(reaction_data[-1].mt_number, 103)
        with self.assertRaises(IndexError):
            reaction_data[5]

    def test_fission_with_yields(self):
        reaction_data = core.ReactionData()
        reaction_data.append(core.CrossSection(U235, core.FISSION_MT, [1.0, 10.0, 100.0], [100.0, 100.0, 100.0]))
        fission_yields = core.FissionYieldData()
        fission_yields[U235] = core.FissionYields(U235, [0.0253], [{BA140: 1.2, MO100: 0.8}])

        rate = 100.0*5e12*BARN_TO_CM2
        to_sink = core.reaction_mat(reaction_data, self.flux).toarray()
        self.assertAlmostEqual(to_sink[index(U235), nis.data.SINK_INDEX], rate, delta=rate*1e-12)

        mat = core.reaction_mat(reaction_data, self.flux, fission_yields).toarray()
        self.assertAlmostEqual(mat[index(U235), index(BA140)], 1.2*rate, delta=rate*1e-12)
        self.assertAlmostEqual(mat[index(U235), index(MO100)], 0.8*rate, delta=rate*1e-12)
        self.assertEqual(mat[index(U235), nis.data.SINK_INDEX], 0.0)

        # change the yields in place
        fission_yields[U235].yields = [{BA140: 2.0}]
        mat = core.reaction_mat(reaction_data, self.flux, fission_yields).toarray()
        self.assertAlmostEqual(mat[index(U235), index(BA140)], 2.0*rate, delta=rate*1e-12)
        self.assertEqual(mat[index(U235), index(MO100)], 0.0)

    def test_yields_need_one_set_for_each_energy(self):
        with self.assertRaises(ValueError):
            core.FissionYields(U235, [0.0253, 1e6], [{BA140: 1.0}])

    def test_products(self):
        self.assertEqual(core.reaction_products(FE56, 102), (FE57, []))
        self.assertTrue(core.is_reaction_mt(102))
        self.assertFalse(core.is_reaction_mt(2))
