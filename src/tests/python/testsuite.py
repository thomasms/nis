import os
import unittest

from .integrationtests import *
from .reactiontests import *


def main():
    unittest.TextTestRunner(verbosity=5).run(unittest.TestSuite())


if __name__ == "__main__":
    main()
