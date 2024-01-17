from collections import Counter

import ENDFtk as tk

tape = tk.tree.Tape.from_file('../../data/serpent/transmutation/transmutation/s2v0_jeff311.dec')

types = Counter()
initial_isomeric_states = Counter()
final_isomeric_states = Counter()

for mat in tape.materials:
    if mat.has_MF(8):
        mf8 = mat.file(8)
        if mf8.has_MT(457):
            parsed = mf8.MT(457).parse()
            for decay_mode in parsed.decay_modes.decay_modes:
                rtype = int(decay_mode.RTYP*1_000_000)
                types[rtype] += 1
                initial_isomeric_states[parsed.isomeric_state] += 1
                final_isomeric_states[decay_mode.final_isomeric_state] += 1

print(initial_isomeric_states)