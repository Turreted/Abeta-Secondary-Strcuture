import os
import argparse
import pandas as pd

des="""
lists the residues in each chain of a strcuture csv file. This can be used to overwrite
the values in constants.py
"""


parser = argparse.ArgumentParser(
    description=des
)
parser.add_argument(
    "-i",
    "--input",
    help="Input Directory (must be the output generate_strcuture.py and contain .csv files)",
    required=True,
)
args = parser.parse_args()

def print_chains(input_dir: str):
    unique_chains = set()
        
    for f in os.listdir(input_dir):
        infile = os.path.join(input_dir, f)
        df = pd.read_csv(infile)
        for cname in set(df["Chain"]):
            residue_list = tuple(df.loc[df["Chain"] == cname]["Residue Name"])
            unique_chains.add(residue_list)

    for chain in unique_chains:
        print(chain)
if __name__ == "__main__":
    print_chains(args.input)
