#!/usr/bin/env python

from ssparser.plotter import plot_monomer_hits, plot_muiltiple, plot_complex_hits
from ssparser.ssparser import parse_ss_freq
from ssparser.scorer import generate_score_table
from ssparser.scorer import generate_raw_pdb_score
from ssparser.utils import get_run_pdb, get_all_pdb
from ssparser.stride import parse_with_stride
from ssparser.constants import *
from ssparser.utils import *

import pandas as pd
import os
import glob

import argparse
from pathlib import Path

"""
Aligns the strcuture of an individual momomer against an aggregate
"""

parser = argparse.ArgumentParser(
    description="Aligns AF-generated strcuture against the aggregate in PDB"
)
parser.add_argument(
    "-i",
    "--input",
    help="Input Directory (must be the output of a AlphaFold run or contain PDB files)",
    required=True,
)
parser.add_argument("-o", "--output", help="Output Directory")
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
parser.add_argument(
    "--plot_freq",
    "-pf",
    default=False,
    type=t_or_f,
    help="Specifies whether to plot a graph of the frequency of secondary strcutures at each position, or to plot the regular graph (of monomer against aggregate)",
)
parser.add_argument(
    "--output_csv",
    "-c",
    default=True,
    type=t_or_f,
    help="Specifies whether or not to generate a scoring table",
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
    input_dir: str,
    output_dir: str,
    is_alphafold=True,
    multimer=True,
    relaxed=True,
    plot_freq=False,
    csv=False,
):
    if not output_dir: output_dir = os.path.join("./output", input_dir)
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    print("Generating frequency table from fibril structures...")
    template_dict = parse_ss_freq(struct_files=absolute_file_paths(DATASET_DIR))

    print("Generating monomer secondary strcuture charts...")
    pdb_list = [os.path.join(input_dir, f) for f in (
        get_run_pdb(input_dir, multimer=multimer, relaxed=relaxed)
        if is_alphafold
        else get_all_pdb(input_dir)
    )]

    # plot frequency dictionary
    if plot_freq:
        # plot ss frequency of individual cimplex
        for pdb in pdb_list:
            d = parse_ss_freq(struct_files=[pdb], filetype="pdb")
            outfile = os.path.join(output_dir, os.path.basename(pdb).replace(".pdb", ".png"))
            print(f"Generating {outfile}...")
            plot_muiltiple(d, output_file=outfile)
    
    # plot individual complexes
    else:
        for pdb in pdb_list:
            df = parse_with_stride(pdb)

            chain_ss_list = []

            for chain_code in set(df["Chain"]):
                chain = df.loc[df["Chain"] == chain_code]
                ss = extract_monomer_ss(chain)
                chain_ss_list.append(ss)

            outfile = os.path.join(output_dir, os.path.basename(pdb).replace(".pdb", ".png"))
            print(f"Generating {outfile}...")
            plot_complex_hits(chain_ss_list, template_dict, output_file=outfile)


if __name__ == "__main__":
    run_analysis(
        args.input,
        args.output,
        is_alphafold=args.alphafold,
        multimer=args.multimer,
        relaxed=args.relaxed,
        plot_freq=args.plot_freq,
        csv=args.output_csv,
    )
