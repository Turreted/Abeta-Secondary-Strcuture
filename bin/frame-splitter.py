#! /usr/bin/env python

import os
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(
    description="Break a single PDB file into multiple files, one for each frame"
)
parser.add_argument(
    "-i",
    "--input",
    help="Input file (must be a PDB file with multiple frames)",
    required=True,
)
parser.add_argument("-o", "--output", help="Output Directory", required=True)
args = parser.parse_args()

def main():
    Path(args.output).mkdir(parents=True, exist_ok=True)
    ofn = args.input.split(".")[0]
    output_file_contents = ""
    i = 1

    for line in open(args.input, "r"):
        if "END" not in line:
            output_file_contents += line
        else:
            with open(os.path.join(args.output, f"{ofn}_{i}.pdb"), "w") as f:
                f.write(output_file_contents)
            i += 1
            output_file_contents = ""

if __name__ == "__main__":
    main()
