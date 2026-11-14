import os
import unittest

from .integrationtests import *
from .solvertests import *
from .datatests import *
from .pulsetests import *
from .reactiontests import *
from .fissiontests import *
from .scripttests import *
from .processingtests import *
from .fnstests import *
from .nisprogramtests import *


def main():
    unittest.TextTestRunner(verbosity=5).run(unittest.TestSuite())


if __name__ == "__main__":
    main()
