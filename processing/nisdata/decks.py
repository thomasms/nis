"""
The inputs of NJOY and of each PREPRO code. They are the inputs of the UKAEA processing for FISPACT-II
(from the dpp prototype), only the file names are changed so that all the steps of one nuclide can run in
one folder.

Fixed format: do not change the spaces.
"""

from pathlib import Path

GROUPS_FILE = Path(__file__).with_name("ukaea1102.txt")
# the number of groups (the file has the 1103 bounds)
NGROUPS = 1102

DEFAULT_TEMPERATURE = 293.6


def njoy(mat: int, temperature: float = DEFAULT_TEMPERATURE) -> str:
    """tape20: the ENDF file. tape35: the PENDF (reconstructed, broadened, with heating and gas production)."""
    return f"""moder
20 -21
reconr
-21 -22
'pendf for material {mat}'/
{mat}/
0.001/
0/
moder
-22 36
broadr
-21 -22 -23
{mat} 1/
0.001/
{temperature}/
0/
unresr
-21 -23 -24
{mat} 1 1 0/
{temperature}/
1.e+10/
0/
moder
-24 32
heatr
-21 -24 -26 30
{mat} 7 0 0 1 2/
302 303 304 318 401 403 407/
heatr
-21 -26 -27 33
{mat} 6 0 0 1 2/
442 443 444 445 446 447/
gaspr
-21 -27 -29
moder
-29 35/
stop
"""


def linear(source: str, output: str, tolerance: str) -> str:
    return f"""          0          0 1.0000d-10          1
{source}
{output}
     0 0  0  999999999
                        (BLANK CARD TERMINATES MAT REQUEST RANGES)
 0.00000-00 {tolerance}
                        (BLANK CARD TERMINATES FILE 3 ERROR LAW)
"""


def recent(source: str, output: str) -> str:
    return f"""          0 1.00000-10          1          1          1          1
{source}
{output}
          1       9999
                        (BLANK CARD TERMINATES MAT REQUEST RANGES)
 0.00000-00 1.00000-04
 1.00000+00 1.00000-04
 2.00000+00 1.00000-03
 2.00000+08 1.00000-03
                        (BLANK CARD TERMINATES FILE 3 ERROR LAW)
"""


def sigma1(source: str, output: str, temperature: float = DEFAULT_TEMPERATURE) -> str:
    return f"""          0          0 {temperature:<11}1.00000-10          1          0
{source}
{output}
          1       9999
                       (BLANK CARD TERMINATES MAT RANGE REQUESTS)
 0.00000+ 0 1.00000-04
 1.00000+ 0 1.00000-04
 2.00000+ 0 1.00000-03
 2.00000+ 8 1.00000-03
                       (BLANK CARD TERMINATES FILE 3 ERROR LAW)
"""


def sixpak(source: str) -> str:
    """Writes ENDFB.MF3, ENDFB.MF4, ENDFB.MF5, ENDFB.MF9, ... (the files that the data has)."""
    return f"""{source}
1    1  9999  999
"""


def merger(output: str, sources: list, requests: list) -> str:
    files = "\n".join(sources)
    lines = "\n".join(requests)
    return f"""{output}
                                                                     0 0
{files}
END
{lines}
                        (blank line terminate request list)
"""


def activate(source: str, output: str) -> str:
    return f"""{source}
{output}
"""


def fixup(source: str, output: str) -> str:
    return f"""10001100110001
{source}
{output}

"""


def dictin(source: str, output: str) -> str:
    return f"""{source}
{output}
"""


def groupie(source: str, output: str) -> str:
    bounds = GROUPS_FILE.read_text().rstrip("\n")
    return f"""         0       {NGROUPS}         0          -1 1.00000-03           0   1
{source}
{output}
          1          1         1           1          1
Groupie UKAEA-1102 groups
     1 1  1  999999999
                    (blank line terminate request list)
{bounds}
"""


# the MERGER requests (MAT, MF and MT ranges) of each merge
MERGE_PREPRO = ["000000 1  100000099999", "000000 1  100000099999"]
MERGE_NJOY = ["     1 1  1  999999220", "     1 1222  999999999"]
MERGE_BOTH = ["000000 1  100000099999", "000000 320300000003447"]
MERGE_SELECT = ["     1 1  1  9999 2999", "     1 3  1  9999 3  2", "     1 3  4  9999 3999", "     1 4  1  999999999"]
MERGE_FINAL = ["000000 1  1000000 3999", "000000 8  1000000 8999", "00000010  100000010999",
               "00000033  100000033999", "00000040  100000040999"]
