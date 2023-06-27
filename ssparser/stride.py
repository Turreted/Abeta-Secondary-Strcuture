import subprocess
import tempfile
import re
import pandas
from .constants import STRIDE_EXE

"""
Generates the structure of a given PDB file using the STRIDE algorithm, which must
be installed to the user's $PATH. Instructions for installing the STRIDE executable
can be found at https://webclu.bio.wzw.tum.de/stride/install.html.
"""

# overwrite with path to your own binary, if necessary
# STRIDE_EXE = "../bin/stride-macos"


def parse_with_stride(infile: str) -> pandas.DataFrame:
    """
    Output Format:
        6-8  Residue name
        10-10 Protein chain identifier
        12-15 PDB	residue	number
        17-20 Ordinal residue number
        25-25 One	letter secondary structure code	**)
        27-39 Full secondary structure name
        43-49 Phi	angle
        53-59 Psi	angle
        65-69 Residue solvent accessible area

        I am so sorry for this code.
    """

    csv_data: pandas.DataFrame

    with tempfile.NamedTemporaryFile(mode="w", delete=True) as tfile:
        # pipe output to tmpfile, wait for it to complete
        p = subprocess.call(
            [STRIDE_EXE, infile], stdout=tfile, stderr=subprocess.STDOUT
        )

        csv_contents = [
            "Residue Name,Chain,PDB ResID,Ordinal ResID,SS Code,SS Full,Phi Angle,Psi Angle,Solvent Area,PDB Code"
        ]
        with open(tfile.name, "r") as f:
            for line in f.readlines():
                if line.startswith("ASG"):
                    csv_contents.append(re.sub(r"\s+", ",", line[5:-1]))

        # lazy hack. 
        with open(tfile.name, "w") as f:
            for line in csv_contents:
                f.write(f"{line}\n")

        csv_data = pandas.read_csv(tfile.name)
    
    return csv_data
