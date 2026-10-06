import argparse
import uproot3 as uproot
import pandas as pd
import numpy as np
import os

parser = argparse.ArgumentParser(
    description="Extract unique run/subrun pairs from a ROOT file."
)
parser.add_argument(
    "path",
    help="Directory containing the ROOT file."
)
parser.add_argument(
    "name",
    help="ROOT filename without the .root extension."
)
parser.add_argument(
    "-o", "--output",
    default=None,
    help="Output text filename (default: <name>.txt)"
)

args = parser.parse_args()

path = args.path
name = args.name

# Ensure the path ends with a slash
if not path.endswith("/"):
    path += "/"

output_file = args.output if args.output else f"{name}.txt"

df_file = uproot.open(os.path.join(path, name + ".root"))["wcpselection"]
df_pot = df_file["T_pot"].pandas.df(
    ["pot_tor875good", "runNo", "subRunNo"],
    flatten=False,
)

df_pot["rs_num"] = (
    df_pot["runNo"].to_numpy() * 100_000_000_000
    + df_pot["subRunNo"].to_numpy() * 1_000_000
)

duplicate_rows = df_pot[df_pot.duplicated(keep=False)]
print(duplicate_rows)

print(df_pot.shape[0])
df_pot = df_pot.drop_duplicates(subset=["rs_num"])
print(df_pot.shape[0])

run = df_pot["runNo"].to_numpy()
subrun = df_pot["subRunNo"].to_numpy()
run_subrun_list = {}

for event in range(len(run)):
    if run[event] in run_subrun_list:
        if subrun[event] in run_subrun_list[run[event]]:
            continue
        run_subrun_list[run[event]].append(subrun[event])
    else:
        run_subrun_list[run[event]] = [subrun[event]]

run_list = sorted(run_subrun_list)

with open(output_file, "w") as f:
    for run_no in run_list:
        for subrun_no in run_subrun_list[run_no]:
            f.write(f"{run_no} {subrun_no}\n")
