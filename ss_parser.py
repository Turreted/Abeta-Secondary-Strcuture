#!/usr/bin/env python

from ssparser.plotter import plot_monomer_hits
from ssparser.ssparser import parse_dataset_freq
from ssparser.scorer import generate_score_table
from ssparser.scorer import generate_raw_pdb_score
from ssparser.utils import get_run_pdb, get_all_pdb
from ssparser.stride import parse_with_stride
from ssparser.constants import *

import pandas as pd

import argparse
from pathlib import Path

"""
Remove all instances of Amyloid-Beta strcutures from the MSA 
step of AlphaFold so we can oberserve the results in isolation
"""


def t_or_f(arg):
    ua = str(arg).upper()
    if "TRUE".startswith(ua):
        return True
    elif "FALSE".startswith(ua):
        return False
    else:
        raise BaseException(f"Argument {arg} is not a boolean value")


parser = argparse.ArgumentParser(
    description="Aligns AF-generated strcuture against the aggregate in PDB"
)
parser.add_argument(
    "-i",
    "--input",
    help="Input Directory (must be the output of a AlphaFold run or contain PDB files)",
    required=True,
)
parser.add_argument("-o", "--output", help="Output Directory", required=True)
parser.add_argument(
    "--alphafold",
    "-a",
    default=True,
    type=t_or_f,
    help="Specifies whether the input directory is the output of an AlphaFold run, which affects if the .pkl metadata is loaded in by the program",
)
parser.add_argument(
    "--multimer",
    "-m",
    default=True,
    type=t_or_f,
    help="Specifies whether the input directory is the output of a multimer run",
)
parser.add_argument(
    "--relaxed",
    "-r",
    default=True,
    type=t_or_f,
    help="Specifies whether the outputs of alphafold have been relaxed (which have different file names)",
)

# TODO
parser.add_argument(
    "--dataset",
    "-d",
    default=True,
    type=t_or_f,
    help="Specifies whether to use full of fibril dataset",
)
parser.add_argument(
    "--output-csv",
    "-c",
    default="table.csv",
    help="Specifies the name of the output data table",
)
args = parser.parse_args()


# extracts monomer secondary structure from dataframe
def extract_monomer_ss(chain: pd.DataFrame) -> list:
    # pad with coils if the monomer is 40 residues instead of 42
    ss = list(chain["SS Code"])
    if len(chain["SS Code"]) == 40:
        ss += ["C", "C"]

    return ss

def run_analysis(
    input_dir: str, output_dir: str, is_alphafold=True, multimer=True, relaxed=True
):
    print("Generating scoring table...")
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    if is_alphafold:
        generate_score_table(input_dir, output_dir)
    else:
        generate_raw_pdb_score(input_dir, output_dir)

    print("Generating frequency table from fibril structures...")
    structure_freq = parse_dataset_freq(dataset_dir=DATASET_DIR)

    print("Generating monomer secondary strcuture charts...")
    # input directory is the result of an alphafold run so contains .pkl metadata files
    if is_alphafold:
        for pdb in get_run_pdb(input_dir, multimer=multimer, relaxed=relaxed):
            fullpath = os.path.join(input_dir, pdb)
            df = parse_with_stride(fullpath)

            # output for each chain in complex, pad with coils
            if multimer:
                for chain_iter, chain_code in enumerate(set(df["Chain"])):
                    chain = df.loc[df["Chain"] == chain_code]
                    ss = extract_monomer_ss(chain)

                    outfile = f"{os.path.join(output_dir, pdb.strip('.pdb'))}-chain-{chain_iter}.png"
                    plot_monomer_hits(ss, structure_freq, output_file=outfile)
                    chain_iter += 1

                    print(f"Output {outfile}")

            # output for single chain in monomer
            else:
                ss = extract_monomer_ss(df)

                outfile = os.path.join(output_dir, pdb.replace("pdb", "png"))
                plot_monomer_hits(ss, structure_freq, output_file=outfile)
                print(f"Output {outfile}")

    # input directory is not the result of an alphafold run and so does not
    # contain .pkl metadata files
    else:
        for pdb in get_all_pdb(input_dir):
            fullpath = os.path.join(input_dir, pdb)
            df = parse_with_stride(fullpath)

            for chain_iter, chain_code in enumerate(set(df["Chain"])):
                chain = df.loc[df["Chain"] == chain_code]
                ss = extract_monomer_ss(chain)

                outfile = f"{os.path.join(output_dir, pdb.strip('.pdb'))}-chain-{chain_iter}.png"
                plot_monomer_hits(ss, structure_freq, output_file=outfile)

                print(f"Output {outfile}")


if __name__ == "__main__":
    run_analysis(
        args.input,
        args.output,
        is_alphafold=args.alphafold,
        multimer=args.multimer,
        relaxed=args.relaxed,
    )
