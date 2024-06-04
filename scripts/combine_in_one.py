"""
Combine separate ENDF files into one big single file
"""
import argparse
import os
from pathlib import Path
from typing import Optional, Union

import ENDFtk as tk


def _try_tape_or_none(file_path: Path) -> Union[None,tk.tree.Tape]:
    tape = None
    try:
        tape = tk.tree.Tape.from_file(str(file_path))
    except Exception as e:
        pass
    return tape

def _try_parse_material_or_none(material: tk.tree.Material) -> Union[None, tk.Material]:
    parsed_mat = None
    try:
        parsed_mat = material.parse()
    except Exception as e:
        pass
    return parsed_mat


def _combine(root_dir: Path) -> tuple[Union[None,tk.tree.Tape], dict]:
    results = {}
    base_tape = None
    for rdir, _, files in os.walk(root_dir):
        for file in files:
            file_path = Path(rdir) / file
            tape = _try_tape_or_none(file_path)
            if tape is None:
                results[file] = 0
                continue
            else:
                if base_tape is None:
                    base_tape = tk.tree.Tape(tape)
                else:
                    for material in tape.materials:
                        parsed_mat = _try_parse_material_or_none(material)
                        if parsed_mat is None:
                            results[material.MAT] = -1
                            continue
                        else:
                            base_tape.insert(material)
                            results[material.MAT] = 1

    return (base_tape, results,)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file_dir")
    parser.add_argument("-o", "--output-file", default="./all.endf")

    args = parser.parse_args()

    file_dir = args.file_dir
    output_file = args.output_file
    full_tape, results = _combine(file_dir)
    if full_tape is not None:
        full_tape.to_file(output_file)

    for k, v in results.items():
        if v != 1:
            print(k)


if __name__ == "__main__":
    main()
