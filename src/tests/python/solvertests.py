"""
The CRAM48 solver against a matrix exponential computed with 60 digits (mpmath), for random decay
matrices with and without loops. No nuclear data is needed, so these run in CI.
"""

import unittest

import numpy as np
import scipy.sparse as sp

import pynis as nis

try:
    import mpmath
except ImportError:  # mpmath is in the dev dependencies
    mpmath = None


def random_matrix(rng: np.random.Generator, size: int, loops: bool) -> np.ndarray:
    """A random (parent, daughter) matrix. Each nuclide decays to 1 to 3 of the nuclides after it, with
    half-lives from 1 us to 300 years, the last nuclide is stable. With loops, some nuclides also go back to
    a nuclide before them at a rate like a reaction rate (1e-8 to 1e-3 /s)."""
    mat = np.zeros((size, size))
    for i in range(size - 1):
        lam = np.log(2) / 10 ** rng.uniform(-6, 10)
        mat[i, i] -= lam
        daughters = rng.choice(np.arange(i + 1, size), size=min(int(rng.integers(1, 4)), size - 1 - i),
                               replace=False)
        weights = rng.uniform(0.01, 1, len(daughters))
        for d, w in zip(daughters, weights / weights.sum()):
            mat[i, d] += lam * w
    if loops:
        for _ in range(3):
            i, j = int(rng.integers(1, size)), int(rng.integers(0, size))
            if i > j:
                rate = 10 ** rng.uniform(-8, -3)
                mat[i, i] -= rate
                mat[i, j] += rate
    return mat


def random_atoms(rng: np.random.Generator, size: int) -> np.ndarray:
    atoms = 10 ** rng.uniform(0, 22, size)
    atoms[rng.uniform(size=size) < 0.3] = 0.0
    if atoms.sum() == 0:
        atoms[0] = 1e20
    return atoms


def exact(mat: np.ndarray, n0: np.ndarray, dt: float) -> np.ndarray:
    """n(dt) = exp(mat^T dt) n0, with 60 digits"""
    with mpmath.workdps(60):
        expm = mpmath.expm(mpmath.matrix((mat.T * dt).tolist()))
        return np.array([float(x) for x in expm * mpmath.matrix(n0.tolist())])


@unittest.skipIf(mpmath is None, "needs mpmath (uv sync installs it with the dev dependencies)")
class Cram48AgainstMatrixExponentialTests(unittest.TestCase):
    STEPS = (1.0, 1e3, 1e6)

    def compare(self, loops: bool, seed: int, trials: int = 10, size: int = 10):
        rng = np.random.default_rng(seed)
        for trial in range(trials):
            mat = random_matrix(rng, size, loops)
            n0 = random_atoms(rng, size)
            total = n0.sum()
            for dt in self.STEPS:
                cram = nis.core.solver_cram48(sp.csc_matrix(mat), n0, dt)
                reference = exact(mat, n0, dt)
                for i in range(size):
                    with self.subTest(trial=trial, dt=dt, nuclide=i):
                        if loops:
                            # a loop mixes the nuclides, so a very small value has the rounding noise of the others
                            tol = max(1e-12 * abs(reference[i]), 1e-14 * total)
                        else:
                            # with no loops each value is precise to its own size, even 1e-40 of the total
                            tol = max(1e-12 * abs(reference[i]), 1e-45 * total)
                        self.assertLessEqual(abs(cram[i] - reference[i]), tol,
                                             f"cram {cram[i]:.17g}, exact {reference[i]:.17g}")

    def test_decay_chains(self) -> None:
        self.compare(loops=False, seed=1)

    def test_with_loops(self) -> None:
        self.compare(loops=True, seed=2)

    def test_larger_matrix(self) -> None:
        self.compare(loops=False, seed=3, trials=2, size=25)


__all__ = ["Cram48AgainstMatrixExponentialTests"]
