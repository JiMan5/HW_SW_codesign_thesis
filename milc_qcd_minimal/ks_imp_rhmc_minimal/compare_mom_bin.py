import os
import filecmp
import numpy as np

#paths
old_dir = "inputs_outputs_bigone_v1"
new_dir = "inputs_outputs_bigone_v2"

#momenta compare before and after
files = [
    "ff_mom_before.bin",
    "ff_mom_after.bin"
]

for fname in files:
    old_path = os.path.join(old_dir, fname)
    new_path = os.path.join(new_dir, fname)

    if not (os.path.exists(old_path) and os.path.exists(new_path)):
        print(f"Missing file: {fname}")
        continue

    #comparison
    if filecmp.cmp(old_path, new_path, shallow=False):
        print(f"{fname}: identical")
    else:
        print(f"{fname}: DIFFERENT")

