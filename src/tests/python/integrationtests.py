from pathlib import Path
import math
import random as rd
import unittest

import pynis as nis


_data_file = Path(__file__).parent.parent.parent.parent / "data" / "jeff311.dec"


def generate_random_nuclide() -> str:
    zai = rd.sample(nis.data.ALL_ZAIS(), 1).pop()
    return nis.zai.get_name_from(zai)


def generate_random_inv() -> dict[str, float]:
    inv = {
        generate_random_nuclide(): math.exp(rd.random() * 65)
        for _ in range(rd.randint(1, 100))
    }
    return inv


class IntegrationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.mat = nis.core.decay_mat_from_file(str(_data_file))
        self.inp = nis.core.InputData()

    def test_h3_halflife(self) -> None:
        self.inp.nuclides = {"h3": 1e20}
        self.inp.dt = 388523520

        y = nis.core.compute_inventory(self.inp, self.mat)
        nuclides = y.nuclides

        # roughly half
        self.assertTrue(nuclides.get("H3") > 4.9e19)
        self.assertTrue(nuclides.get("He3") > 4.9e19)

    def test_random_inv(self) -> None:
        for _ in range(100):
            inv = generate_random_inv()
            # 22 keeps this to a max of ~41.8 years
            dt = math.exp(rd.random() * 22)
            self.inp.nuclides = inv
            self.inp.dt = dt
            y = nis.core.compute_inventory(self.inp, self.mat)
            self.assertTrue(len(y.nuclides) > 0)


__all__ = ["IntegrationTests"]
