import argparse
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import math
from pypact.library.nuclidelib import (
    get_zai,
    get_name,
    get_zai_props,
    NUCLIDE_DICTIONARY,
)
from typing import Callable, Union

import ENDFtk as tk

LN_2 = math.log(2)
ZAI_ALPHA = 20040
ZAI_SINK = -1

# CRAM coefficients
theta_r = np.array(
    [
        -4.465731934165702e1,
        -5.284616241568964e0,
        -8.867715667624458e0,
        +3.493013124279215e0,
        +1.564102508858634e1,
        +1.742097597385893e1,
        -2.834466755180654e1,
        +1.661569367939544e1,
        +8.011836167974721e0,
        -2.056267541998229e0,
        +1.449208170441839e1,
        +1.853807176907916e1,
        +9.932562704505182e0,
        -2.244223871767187e1,
        +8.590014121680897e-1,
        -1.286192925744479e1,
        +1.164596909542055e1,
        +1.806076684783089e1,
        +5.870672154659249e0,
        -3.542938819659747e1,
        +1.901323489060250e1,
        +1.885508331552577e1,
        -1.734689708174982e1,
        +1.316284237125190e1,
    ]
)

theta_i = np.array(
    [
        +6.233225190695437e1,
        +4.057499381311059e1,
        +4.325515754166724e1,
        +3.281615453173585e1,
        +1.558061616372237e1,
        +1.076629305714420e1,
        +5.492841024648724e1,
        +1.316994930024688e1,
        +2.780232111309410e1,
        +3.794824788914354e1,
        +1.799988210051809e1,
        +5.974332563100539e0,
        +2.532823409972962e1,
        +5.179633600312162e1,
        +3.536456194294350e1,
        +4.600304902833652e1,
        +2.287153304140217e1,
        +8.368200580099821e0,
        +3.029700159040121e1,
        +5.834381701800013e1,
        +1.194282058271408e0,
        +3.583428564427879e0,
        +4.883941101108207e1,
        +2.042951874827759e1,
    ]
)

c48_theta = np.array(theta_r + theta_i * 1j, dtype=np.complex128)

alpha_r = np.array(
    [
        +6.387380733878774e2,
        +1.909896179065730e2,
        +4.236195226571914e2,
        +4.645770595258726e2,
        +7.765163276752433e2,
        +1.907115136768522e3,
        +2.909892685603256e3,
        +1.944772206620450e2,
        +1.382799786972332e5,
        +5.628442079602433e3,
        +2.151681283794220e2,
        +1.324720240514420e3,
        +1.617548476343347e4,
        +1.112729040439685e2,
        +1.074624783191125e2,
        +8.835727765158191e1,
        +9.354078136054179e1,
        +9.418142823531573e1,
        +1.040012390717851e2,
        +6.861882624343235e1,
        +8.766654491283722e1,
        +1.056007619389650e2,
        +7.738987569039419e1,
        +1.041366366475571e2,
    ]
)

alpha_i = np.array(
    [
        -6.743912502859256e2,
        -3.973203432721332e2,
        -2.041233768918671e3,
        -1.652917287299683e3,
        -1.783617639907328e4,
        -5.887068595142284e4,
        -9.953255345514560e3,
        -1.427131226068449e3,
        -3.256885197214938e6,
        -2.924284515884309e4,
        -1.121774011188224e3,
        -6.370088443140973e4,
        -1.008798413156542e6,
        -8.837109731680418e1,
        -1.457246116408180e2,
        -6.388286188419360e1,
        -2.195424319460237e2,
        -6.719055740098035e2,
        -1.693747595553868e2,
        -1.177598523430493e1,
        -4.596464999363902e3,
        -1.738294585524067e3,
        -4.311715386228984e1,
        -2.777743732451969e2,
    ]
)

c48_alpha = np.array(alpha_r + alpha_i * 1j, dtype=np.complex128)

c48_alpha0 = 2.258038182743983e-47


def compute_zai(z: int, a: int, i: int) -> int:
    return z * 10000 + a * 10 + i


def change_zai(zai: int, z_diff: int, a_diff: int, i_diff: int) -> int:
    z, a, i = get_zai_props(zai)
    n_z = int(z + z_diff)
    n_a = int(a + a_diff)
    n_i = int(i + i_diff)
    new_zai = compute_zai(n_z, n_a, n_i)
    return new_zai


def parse_input_file(file_path):
    data = {}
    with open(file_path, "r") as file:
        lines = file.readlines()
        for iline, line in enumerate(lines):
            line = line.strip().strip("\n")
            if line.startswith("#"):
                continue

            if line.startswith("dt"):
                dt = float(line.split()[1])
                data["dt"] = dt
            elif line == "nuclides":
                nuclides = {}
                for nuclide_line in lines[iline + 1 :]:
                    nuclide_line = nuclide_line.strip()
                    if not nuclide_line or nuclide_line.startswith("dt"):
                        break
                    if nuclide_line.startswith("#"):
                        continue

                    nuclide, amount = nuclide_line.split()
                    nuclides[nuclide] = float(amount)
                data["nuclides"] = nuclides
    return data


reaction_product_lookup_func: dict[int, Callable] = {
    0: lambda parent_zai: (
        parent_zai,
        [],
    ),
    1: lambda parent_zai: (
        change_zai(parent_zai, 1, 0, 0),
        [],
    ),
    2: lambda parent_zai: (
        change_zai(parent_zai, -1, 0, 0),
        [],
    ),
    # need to handle isomeric state for IT transistion
    3: lambda parent_zai: (
        change_zai(parent_zai, 0, 0, 0),
        [],
    ),
    4: lambda parent_zai: (
        change_zai(parent_zai, -2, -4, 0),
        [ZAI_ALPHA],
    ),
    5: lambda parent_zai: (
        change_zai(parent_zai, 0, -1, 0),
        [],
    ),
    # SF is tricky: need to implement
    6: lambda parent_zai: (
        ZAI_SINK,
        [],
    ),
    7: lambda parent_zai: (
        change_zai(parent_zai, -1, 0, 0),
        [],
    ),
    8: lambda parent_zai: (
        change_zai(parent_zai, 1, 0, 0),
        [],
    ),
    9: lambda parent_zai: (
        parent_zai,
        [],
    ),
    10: lambda parent_zai: (
        ZAI_SINK,
        [],
    ),
}


def get_digits(number: Union[float, int]) -> list[int]:
    if number == 0:
        return [0]

    sign_removed = max(number, -number) if number < 0 else number
    value = sign_removed

    digits = []
    while value > 0:
        digit = value % 10
        digits.append(digit)
        value /= 10

    return list(reversed(digits))


def expand_rtypes(rtyp: float) -> list[int]:
    _rtyp_factor = 1_000_000
    initial_rtype = int(math.floor(rtyp))

    # if we are zero first then we are gamma
    if initial_rtype == 0:
        return [0]

    # if the initial rtype is not valid
    # then return an empty vector
    # regardless of what else comes next
    if initial_rtype not in reaction_product_lookup_func.keys():
        return []

    # all other cases we process the subsequent values
    rtyp_as_int = int(rtyp * _rtyp_factor)
    expanded_rtype = rtyp_as_int - (initial_rtype * _rtyp_factor)
    rtypes = get_digits(expanded_rtype)
    # we also need to remove trailing zeros,
    # as this does not mean GammaRay emission
    # it is the end
    reaction_types = [initial_rtype]
    for v in rtypes:
        if (v == 0) or (v not in reaction_product_lookup_func.keys()):
            break

        reaction_types.append(v)

    return reaction_types


def compute_products(
    parent_zai: int, rtyp: float, final_isomeric_state: int
) -> tuple[int, list[int]]:
    rtypes = expand_rtypes(rtyp)

    # we need to trim remaining zeros, as there is no way to tell between
    # 1.0 (beta -> gamma) and 1.0 (beta)
    primary_zai = parent_zai
    secondaries = []
    for reaction_type in rtypes:
        if reaction_type in reaction_product_lookup_func.keys():
            dpc = reaction_product_lookup_func.get(reaction_type)
            # we assume - and because the data tells us this -
            # that SF (6) doesn't have any modes following it
            # to make this simpler with the primary
            # Perhaps a tree is a better data structure here....
            prods = dpc(primary_zai)
            primary_zai = prods[0]
            for secondary in prods[1]:
                secondaries.append(secondary)

    # finally check the final isomeric state
    if primary_zai > 0:
        z, a, i = get_zai_props(primary_zai)
        primary_zai = compute_zai(z, a, final_isomeric_state)

    return (
        primary_zai,
        secondaries,
    )


def get_zai_lookups(all_zais: list[int]) -> tuple[dict[int, int], dict[int, int]]:
    zai_lookup = {zai: i for i, zai in enumerate(all_zais)}
    mat_n = len(all_zais) + 1
    sink_index = mat_n - 1
    zai_lookup[ZAI_SINK] = sink_index

    zai_lookup_r = {v: k for k, v in zai_lookup.items()}

    return zai_lookup, zai_lookup_r


def make_dense_matrix(tape: tk.tree.Tape, all_zais: list[int]) -> np.array:
    zai_lookup, _ = get_zai_lookups(all_zais)
    mat_n = len(all_zais) + 1
    sink_index = mat_n - 1

    matrix = np.zeros((mat_n, mat_n))

    for mat in tape.materials:
        if mat.has_MF(8):
            mf8 = mat.file(8)
            if mf8.has_MT(457):
                parsed = mf8.MT(457).parse()
                for decay_mode in parsed.decay_modes.decay_modes:
                    rtype = int(decay_mode.RTYP * 1_000_000)
                    za = parsed.ZA
                    i = parsed.isomeric_state
                    zai = za * 10 + i
                    is_stable = parsed.is_stable
                    if not is_stable and (zai in zai_lookup):
                        za_index = zai_lookup.get(zai, sink_index)
                        half_life = parsed.half_life[0]
                        decay_const = LN_2 / half_life

                        # add on diagonal
                        matrix[za_index, za_index] -= decay_const

                        for decay_mode in parsed.decay_modes.decay_modes:
                            reaction_type = decay_mode.RTYP
                            final_isomeric_state = math.floor(
                                decay_mode.final_isomeric_state
                            )
                            products = compute_products(
                                zai, reaction_type, final_isomeric_state
                            )
                            primary, secondaries = products
                            secondaries.append(primary)
                            br = decay_mode.branching_ratio[0]
                            # print(zai, get_name(zai), half_life, reaction_type, products)
                            for daughter_zai in secondaries:
                                dzai_index = zai_lookup.get(daughter_zai, sink_index)
                                matrix[za_index, dzai_index] += br * decay_const

    return matrix


def cram_solve(nvec: np.array, matrix: np.array, dt: float) -> np.array:
    mat_sparse = sp.csc_matrix(matrix.T, dtype=np.float64)
    A = dt * mat_sparse
    y = nvec.copy()
    ident = sp.eye(A.shape[0], format="csc")
    for alpha, theta in zip(c48_alpha, c48_theta):
        lhs = A - theta * ident
        solved = sla.spsolve(lhs, y)
        y += 2 * np.real(alpha * solved)
    return y * c48_alpha0


def main():
    parser = argparse.ArgumentParser()
    # should be ../../data/serpent/transmutation/transmutation/s2v0_jeff311.dec
    parser.add_argument(
        "-d",
        "--data",
    )
    # see example files
    parser.add_argument("-i", "--input_file")
    args = parser.parse_args()

    all_zais = [
        get_zai(f"{entry['element']}{isotope}")
        for entry in NUCLIDE_DICTIONARY
        for isotope in entry["isotopes"]
    ]

    zai_lookup, zai_lookup_r = get_zai_lookups(all_zais)
    tape = tk.tree.Tape.from_file(args.data)
    matrix = make_dense_matrix(tape, all_zais)

    input_data = parse_input_file(args.input_file)
    dt = input_data.get("dt")

    min_natoms = 1e12
    n0 = np.zeros((matrix.shape[0]))
    for nname, natoms in input_data.get("nuclides", {}).items():
        zai = get_zai(nname)
        if zai is None:
            zai = ZAI_SINK
        n_index = zai_lookup.get(zai)
        n0[n_index] = natoms

    # timestepping not so great after 100 years
    # limit to 5e9 ~ 158 years
    # anything above this length is simply infeasible
    if abs(dt) > 5e9:
        raise Exception("Cannot go that far ahead! Limit to 5e9 seconds (158 years)")

    nt = cram_solve(n0, matrix, dt)

    print("=====================")
    print("       forward       ")
    print("=====================")
    for i, nx in enumerate(nt):
        if nx > min_natoms:
            zai = zai_lookup_r.get(i)
            print(f"{get_name(zai):<10}", f"{nx:>10.4g}")

    # run the output back again
    print("=====================")
    print("       reverse       ")
    print("=====================")
    n0_2 = cram_solve(nt, matrix, -dt)
    for i, nx in enumerate(n0_2):
        if nx > min_natoms:
            zai = zai_lookup_r.get(i)
            print(f"{get_name(zai):<10}", f"{nx:>10.4g}")


if __name__ == "__main__":
    main()
