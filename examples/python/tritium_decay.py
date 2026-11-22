"""
Tritium decay, with decay data made by hand. No nuclear data files are needed.

H3 decays by beta- to He3, with a half-life of 12.32 years.
We make the decay data, make the decay matrix from it, then solve 5 steps of one half-life
with compute_steps. Then we compare the result with the analytic answer N(t) = N0 * exp(-lambda*t).

Run it with:
    uv run python examples/python/tritium_decay.py
"""

import math

import pynis as nis

core = nis.core

H3 = nis.zai.get_from_name("H3")    # 10030
HE3 = nis.zai.get_from_name("He3")  # 20030

YEAR = 365.0*24.0*3600.0
H3_HALF_LIFE = 12.32*YEAR
# mean beta energy of H3 decay in eV (the light particle energy)
H3_BETA_ENERGY = 5.68e3

N0 = 1e20


def tritium_decay_data():
    """H3 -> He3 by beta- (RTYP=1), all decays (branching ratio 1). He3 is stable."""
    decay_data = core.DecayData()
    decay_data[H3] = core.DecayNuclide(half_life=H3_HALF_LIFE,
                                       modes=[core.DecayMode(rtyp=1.0, branching_ratio=1.0)],
                                       light_particle=H3_BETA_ENERGY)
    decay_data[HE3] = core.DecayNuclide()
    return decay_data


def solve(decay_mat, nuclides, dt):
    model = core.BasicModel()
    model.nuclides = nuclides
    model.dt = dt
    return core.compute_inventory(model, decay_mat, core.solver_cram48).nuclides


def main():
    decay_data = tritium_decay_data()
    decay_mat = core.decay_mat(decay_data)
    energy_data = core.decay_energy_data(decay_data)
    lam = decay_data[H3].decay_constant

    print(f"H3 half-life {decay_data[H3].half_life/YEAR:.2f} y, decay constant {lam:.6e} 1/s\n")
    print(f"{'time (y)':>9} {'H3':>13} {'He3':>13} {'H3 exact':>13} {'rel diff':>10} {'activity (Bq)':>14} {'heat (W)':>11}")

    # 5 steps of one half-life, each step starts from the result of the step before
    model = core.BasicModel()
    model.nuclides = {"H3": N0}
    model.steps = [core.TimeStep(H3_HALF_LIFE) for _ in range(5)]
    results = core.compute_steps(model, decay_mat)

    for n, nuclides in enumerate([model.nuclides] + [r.nuclides for r in results]):
        t = n*H3_HALF_LIFE
        h3 = nuclides.get("H3", 0.0)
        he3 = nuclides.get("He3", 0.0)
        exact = N0*math.exp(-lam*t)

        inventory = core.vector_from_inventory(nuclides)
        activity = core.compute_activity(energy_data, inventory)
        heat = core.compute_decay_heat(energy_data, inventory).total

        print(f"{t/YEAR:9.2f} {h3:13.6e} {he3:13.6e} {exact:13.6e} {abs(h3 - exact)/exact:10.1e} {activity:14.6e} {heat:11.4e}")

    # change the data in place: a half-life of 1 year, then solve again for 1 year
    decay_data[H3].half_life = YEAR
    nuclides = solve(core.decay_mat(decay_data), {"H3": N0}, YEAR)
    print(f"\nwith a half-life of 1 y, H3 after 1 y: {nuclides['H3']/N0:.6f} of N0 (exact 0.5)")


if __name__ == "__main__":
    main()
