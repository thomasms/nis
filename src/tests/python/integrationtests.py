from pathlib import Path
import math
import random as rd
import unittest

import pynis as nis


_jeff311_data_file = Path(__file__).parent.parent.parent.parent / "data" / "jeff311.dec"
_endfbvii_data_file = Path(__file__).parent.parent.parent.parent / "data" / "endfbvii1.dec"


def generate_random_nuclide() -> str:
    zai = rd.sample(nis.data.ALL_ZAIS(), 1).pop()
    return nis.zai.get_name_from(zai)


def generate_random_inv() -> dict[str, float]:
    inv = {
        generate_random_nuclide(): math.exp(rd.random() * 65)
        for _ in range(rd.randint(1, 100))
    }
    return inv


class IntegrationBaseTests(unittest.TestCase):
    def setUp(self) -> None:
        self.mat = nis.core.decay_mat_from_file(str(_jeff311_data_file))
        self.model = nis.core.BasicModel()

    def test_h3_halflife(self) -> None:
        self.model.nuclides = {"h3": 1e20}
        self.model.dt = 388523520

        y = nis.core.compute_inventory(self.model, self.mat)
        nuclides = y.nuclides

        # roughly half
        self.assertTrue(nuclides.get("H3") > 4.9e19)
        self.assertTrue(nuclides.get("He3") > 4.9e19)

    def test_random_inv(self) -> None:
        for _ in range(100):
            inv = generate_random_inv()
            # 22 keeps this to a max of ~41.8 years
            dt = math.exp(rd.random() * 22)
            self.model.nuclides = inv
            self.model.dt = dt
            y = nis.core.compute_inventory(self.model, self.mat)
            self.assertTrue(len(y.nuclides) > 0)


class IntegrationJeff311Tests(IntegrationBaseTests):
    def setUp(self) -> None:
        self.mat = nis.core.decay_mat_from_file(str(_jeff311_data_file))
        self.inp = nis.core.BasicModel()


class IntegrationEndfbVIITests(IntegrationBaseTests):
    def setUp(self) -> None:
        self.mat = nis.core.decay_mat_from_file(str(_endfbvii_data_file))
        self.inp = nis.core.BasicModel()


__all__ = ["IntegrationJeff311Tests", "IntegrationEndfbVIITests"]
