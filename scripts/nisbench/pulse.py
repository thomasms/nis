"""
Decay heat after a fission pulse: the fission products of one fission (the independent yields) at t = 0,
then decay only. This is the quantity of the fission pulse experiments (Dickens, Akiyama, Tobias, ...),
usually given as t * P(t) in MeV per fission.

A pulse case in the benchmark format (a short irradiation with a flux) gives the same result for cooling
times much longer than the irradiation. This module does not need the reaction data.
"""

import math
from dataclasses import dataclass

MEV_TO_J = 1.602176634e-13

# the decay heat of the default times goes from 0.1 s to ~30 years
DEFAULT_TIMES = [10 ** (k / 10) for k in range(-10, 91)]


@dataclass
class PulseHeat:
    """The decay heat per fission at a time after the pulse, in MeV/s."""
    time_s: float
    beta: float
    gamma: float
    alpha: float

    @property
    def total(self) -> float:
        return self.beta + self.gamma + self.alpha


def pulse_inventory(fission_yields, parent: str, energy_ev: float = 0.0253) -> dict:
    """The atoms of each fission product made by one fission of the parent (a name, e.g. U235) by a neutron
    of energy_ev, from a fission yield file (pynis.core.fission_yields_from_file)."""
    import pynis as nis
    zai = nis.zai.get_from_name(parent)
    if zai not in fission_yields:
        raise ValueError(f"{parent}: no fission yields in the file")
    products = nis.core.fission_yields_at_energy(fission_yields[zai], energy_ev)
    return {nis.zai.get_name_from(z): y for z, y in products.items() if y > 0}


def pulse_decay_heat(decay_mat, energy_data, inventory: dict, times=None) -> list:
    """The decay heat (MeV/s per fission) at each time (s, going up) after the pulse.
    Each step goes from the time before, so the result does not depend on how the times are spaced."""
    import pynis as nis
    times = DEFAULT_TIMES if times is None else list(times)
    if any(b <= a for a, b in zip(times, times[1:])) or times[0] < 0:
        raise ValueError("the times must be 0 or more, and going up")
    results = []
    nuclides, previous = dict(inventory), 0.0
    for t in times:
        if t > previous:
            model = nis.core.BasicModel()
            model.nuclides = nuclides
            model.dt = t - previous
            nuclides = nis.core.compute_inventory(model, decay_mat, nis.core.solver_cram48).nuclides
        previous = t
        heat = nis.core.compute_decay_heat(energy_data, nis.core.vector_from_inventory(nuclides))
        results.append(PulseHeat(t, heat.light_particle / MEV_TO_J, heat.electromagnetic / MEV_TO_J,
                                 heat.heavy_particle / MEV_TO_J))
    return results


def released_energy(decay_mat, energy_data, inventory: dict, time_s: float = math.inf,
                    max_half_life_s: float = 1e15) -> tuple:
    """The (beta, gamma) energy in MeV that the decays give from t = 0 to time_s (all time if inf).

    The integral of the atoms of the radionuclides A is solved directly: A^T x = n(t) - n(0), with n(t) = 0 for
    all time. Nuclides with a half-life longer than max_half_life_s are left out: some libraries give the
    double-beta emitters (Zr96, Mo100, Xe136, ...) a half-life of ~1e27 s, which adds ~0.9 MeV of beta energy
    over all time that is never seen."""
    import numpy as np
    import scipy.sparse.linalg as spl
    import pynis as nis
    lam = np.asarray(energy_data.decay_constants)
    radioactive = np.nonzero(lam > math.log(2) / max_half_life_s)[0]
    n0 = np.asarray(nis.core.vector_from_inventory(inventory))
    change = -n0[radioactive]
    if math.isfinite(time_s):
        change = change + np.asarray(nis.core.solver_cram48(decay_mat, n0, time_s))[radioactive]
    a = decay_mat.tocsc().T.tocsr()[radioactive][:, radioactive].tocsc()
    decays = lam[radioactive] * spl.spsolve(a, change)
    beta = decays @ np.asarray(energy_data.light_particle)[radioactive] / 1e6
    gamma = decays @ np.asarray(energy_data.electromagnetic)[radioactive] / 1e6
    return float(beta), float(gamma)


def delayed_energy(decay_mat, energy_data, inventory: dict, max_half_life_s: float = 1e15) -> tuple:
    """The (beta, gamma) energy in MeV that the decays of the fission products give in all time after the pulse:
    the integral of the decay heat. This is EB and EGD of ENDF MF=1 MT=458."""
    return released_energy(decay_mat, energy_data, inventory, math.inf, max_half_life_s)
