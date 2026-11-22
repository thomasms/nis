"""
Tests for compute_steps in python: irradiation and cooling steps, with data made by hand.
No data files are needed.
"""

import math
import unittest

import pynis as nis

core = nis.core

CO59 = 270590
CO60 = 270600
NI60 = 280600

CO60_HALF_LIFE = 1.663415e8
BARN_TO_CM2 = 1e-24

# 2 groups, 1e12 and 4e12 n/cm^2/s, a total of 5e12 n/cm^2/s
SPECTRUM = core.Spectrum([1.0, 10.0, 100.0], [1e12, 4e12])
SPECTRUM_TOTAL = 5e12
# Co59 (n,g) Co60 of 10 barns in both groups
CO59_RATE = 10.0*SPECTRUM_TOTAL*BARN_TO_CM2


def decay_mat():
    decay_data = core.DecayData()
    decay_data[CO60] = core.DecayNuclide(half_life=CO60_HALF_LIFE, modes=[core.DecayMode(1.0, 1.0)])
    return core.decay_mat(decay_data)


def reaction_mat():
    reaction_data = core.ReactionData()
    reaction_data.append(core.CrossSection(CO59, 102, [1.0, 10.0, 100.0], [10.0, 10.0, 10.0]))
    return core.reaction_mat(reaction_data, SPECTRUM)


def model(nuclides, steps, spectrum=None):
    m = core.BasicModel()
    m.nuclides = nuclides
    m.steps = steps
    m.dt = steps[0].dt
    if spectrum is not None:
        m.spectrum = spectrum
    return m


class ComputeStepsTest(unittest.TestCase):

    def test_time_steps(self):
        m = core.BasicModel()
        m.dt = 10.0
        steps = core.time_steps(m)
        self.assertEqual(len(steps), 1)
        self.assertEqual(steps[0].dt, 10.0)
        self.assertIsNone(steps[0].flux)

    def test_flux_scale(self):
        self.assertEqual(core.flux_scale(core.TimeStep(1.0), SPECTRUM_TOTAL), 1.0)
        self.assertEqual(core.flux_scale(core.TimeStep(1.0, 0.0), SPECTRUM_TOTAL), 0.0)
        # the step flux is in n/m^2/s, the spectrum in n/cm^2/s
        self.assertAlmostEqual(core.flux_scale(core.TimeStep(1.0, 1e17), SPECTRUM_TOTAL), 2.0, delta=1e-15)

    def test_decay_only(self):
        """two steps of one half-life, each from the result of the step before"""
        m = model({"Co60": 1e20}, [core.TimeStep(CO60_HALF_LIFE), core.TimeStep(CO60_HALF_LIFE)])
        results = core.compute_steps(m, decay_mat())
        self.assertEqual(len(results), 2)
        self.assertAlmostEqual(results[0].nuclides["Co60"]/1e20, 0.5, delta=1e-12)
        self.assertAlmostEqual(results[1].nuclides["Co60"]/1e20, 0.25, delta=1e-12)
        self.assertAlmostEqual(results[1].nuclides["Ni60"]/1e20, 0.75, delta=1e-12)

    def test_irradiation_then_cooling(self):
        irradiation = 1e6
        cooling = CO60_HALF_LIFE
        # 2 times the flux of the spectrum
        flux = 2.0*SPECTRUM_TOTAL*1e4
        m = model({"Co59": 1e20}, [core.TimeStep(irradiation, flux), core.TimeStep(cooling, 0.0)], SPECTRUM)
        results = core.compute_steps(m, decay_mat(), reaction_mat())
        self.assertEqual(len(results), 2)

        rate = 2.0*CO59_RATE
        lam = math.log(2)/CO60_HALF_LIFE
        co59 = 1e20*math.exp(-rate*irradiation)
        # Bateman for Co59 -> Co60 -> Ni60
        co60 = 1e20*rate/(lam - rate)*(math.exp(-rate*irradiation) - math.exp(-lam*irradiation))
        self.assertAlmostEqual(results[0].nuclides["Co59"]/co59, 1.0, delta=1e-10)
        self.assertAlmostEqual(results[0].nuclides["Co60"]/co60, 1.0, delta=1e-10)

        # no reactions in cooling: Co59 stays the same, Co60 decays for one half-life
        self.assertAlmostEqual(results[1].nuclides["Co59"]/co59, 1.0, delta=1e-10)
        self.assertAlmostEqual(results[1].nuclides["Co60"]/co60, 0.5, delta=1e-10)

    def test_spectrum_total_from_the_model(self):
        steps = [core.TimeStep(1e6, 1e17)]
        given = core.compute_steps(model({"Co59": 1e20}, steps), decay_mat(), reaction_mat(), SPECTRUM_TOTAL)
        from_model = core.compute_steps(model({"Co59": 1e20}, steps, SPECTRUM), decay_mat(), reaction_mat())
        self.assertAlmostEqual(given[0].nuclides["Co60"]/from_model[0].nuclides["Co60"], 1.0, delta=1e-15)

    def test_no_flux_is_the_spectrum_as_given(self):
        m = model({"Co59": 1e20}, [core.TimeStep(1e6)])
        results = core.compute_steps(m, decay_mat(), reaction_mat())
        one = core.compute_inventory(m, decay_mat() + reaction_mat(), core.solver_cram48)
        self.assertAlmostEqual(results[0].nuclides["Co60"]/one.nuclides["Co60"], 1.0, delta=1e-14)

    def test_same_as_compute_inventory_with_the_solver(self):
        m = model({"Co60": 1e20}, [core.TimeStep(1e7)])
        results = core.compute_steps(m, decay_mat(), solver=core.solver_cram48)
        one = core.compute_inventory(m, decay_mat(), core.solver_cram48)
        self.assertEqual(results[0].nuclides["Co60"], one.nuclides["Co60"])

    def test_flux_needs_a_spectrum(self):
        m = model({"Co59": 1e20}, [core.TimeStep(1e6, 1e17)])
        with self.assertRaises(ValueError):
            core.compute_steps(m, decay_mat(), reaction_mat())


if __name__ == "__main__":
    unittest.main()
