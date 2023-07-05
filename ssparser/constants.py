import pathlib
import os

"""
Contains constants for parsing and plotting each ABeta monomer
"""

file_cwd = pathlib.Path(__file__).parent.resolve()
DATASET_DIR = os.path.join(file_cwd, "../datasets/M1_beus")
STRIDE_EXE = os.path.join(file_cwd, "../bin/stride-macos")

ABETA_LONG_SEQUENCE = ('ASN', 'ASN', 'PRO', 'ALA', 'ILE', 'LYS', 'ARG', 'ILE', 'GLY', 'ASN', 'HSD', 'ILE', 'THR', 'LYS', 'SER', 'PRO', 'GLU', 'ASP', 'LYS', 'ARG', 'GLU', 'TYR', 'ARG', 'GLY', 'LEU', 'GLU', 'LEU', 'ALA', 'ASN', 'GLY', 'ILE', 'LYS', 'VAL', 'LEU', 'LEU', 'ILE', 'SER', 'ASP', 'PRO', 'THR', 'THR', 'ASP', 'LYS', 'SER', 'SER', 'ALA', 'ALA', 'LEU', 'ASP', 'VAL', 'HSD', 'ILE', 'GLY', 'SER', 'LEU', 'SER', 'ASP', 'PRO', 'PRO', 'ASN', 'ILE', 'ALA', 'GLY', 'LEU', 'SER', 'HSD', 'PHE', 'LEU', 'GLU', 'HSD', 'MET', 'LEU', 'PHE', 'LEU', 'GLY', 'THR', 'LYS', 'LYS', 'TYR', 'PRO', 'LYS', 'GLU', 'ASN', 'GLU', 'TYR', 'SER', 'GLN', 'PHE', 'LEU', 'SER', 'GLU', 'HSD', 'ALA', 'GLY', 'SER', 'SER', 'ASN', 'ALA', 'PHE', 'THR', 'SER', 'GLY', 'GLU', 'HSD', 'THR', 'ASN', 'TYR', 'TYR', 'PHE', 'ASP', 'VAL', 'SER', 'HSD', 'GLU', 'HSD', 'LEU', 'GLU', 'GLY', 'ALA', 'LEU', 'ASP', 'ARG', 'PHE', 'ALA', 'GLN', 'PHE', 'PHE', 'LEU', 'SER', 'PRO', 'LEU', 'PHE', 'ASP', 'GLU', 'SER', 'ALA', 'LYS', 'ASP', 'ARG', 'GLU', 'VAL', 'ASN', 'ALA', 'VAL', 'ASP', 'SER', 'GLU', 'HSD', 'GLU', 'LYS', 'ASN', 'VAL', 'MET', 'ASN', 'ASP', 'ALA', 'TRP', 'ARG', 'LEU', 'PHE', 'GLN', 'LEU', 'GLU', 'LYS', 'ALA', 'THR', 'GLY', 'ASN', 'PRO', 'LYS', 'HSD', 'PRO', 'PHE', 'SER', 'LYS', 'PHE', 'GLY', 'THR', 'GLY', 'ASN', 'LYS', 'TYR', 'THR', 'LEU', 'GLU', 'THR', 'ARG', 'PRO', 'ASN', 'GLN', 'GLU', 'GLY', 'ILE', 'ASP', 'VAL', 'ARG', 'GLN', 'GLU', 'LEU', 'LEU', 'LYS', 'PHE', 'HSD', 'SER', 'ALA', 'TYR', 'TYR', 'SER', 'SER', 'ASN', 'LEU', 'MET', 'ALA', 'VAL', 'VAL', 'VAL', 'LEU', 'GLY', 'ARG', 'GLU', 'SER', 'LEU', 'ASP', 'ASP', 'LEU', 'THR', 'ASN', 'LEU', 'VAL', 'VAL', 'LYS', 'LEU', 'PHE', 'SER', 'GLU', 'VAL', 'GLU', 'ASN', 'LYS', 'ASN', 'VAL', 'PRO', 'LEU', 'PRO', 'GLU', 'PHE', 'PRO', 'GLU', 'HSD', 'PRO', 'PHE', 'GLN', 'GLU', 'GLU', 'HSD', 'LEU', 'LYS', 'GLN', 'LEU', 'TYR', 'LYS', 'ILE', 'VAL', 'PRO', 'ILE', 'LYS', 'ASP', 'ILE', 'ARG', 'ASN', 'LEU', 'TYR', 'VAL', 'THR', 'PHE', 'PRO', 'ILE', 'PRO', 'ASP', 'LEU', 'GLN', 'LYS', 'TYR', 'TYR', 'LYS', 'SER', 'ASN', 'PRO', 'GLY', 'HSD', 'TYR', 'LEU', 'GLY', 'HSD', 'LEU', 'ILE', 'GLY', 'HSD', 'GLU', 'GLY', 'PRO', 'GLY', 'SER', 'LEU', 'LEU', 'SER', 'GLU', 'LEU', 'LYS', 'SER', 'LYS', 'GLY', 'TRP', 'VAL', 'ASN', 'THR', 'LEU', 'VAL', 'GLY', 'GLY', 'GLN', 'LYS', 'GLU', 'GLY', 'ALA', 'ARG', 'GLY', 'PHE', 'MET', 'PHE', 'PHE', 'ILE', 'ILE', 'ASN', 'VAL', 'ASP', 'LEU', 'THR', 'GLU', 'GLU', 'GLY', 'LEU', 'LEU', 'HSD', 'VAL', 'GLU', 'ASP', 'ILE', 'ILE', 'LEU', 'HSD', 'MET', 'PHE', 'GLN', 'TYR', 'ILE', 'GLN', 'LYS', 'LEU', 'ARG', 'ALA', 'GLU', 'GLY', 'PRO', 'GLN', 'GLU', 'TRP', 'VAL', 'PHE', 'GLN', 'GLU', 'LEU', 'LYS', 'ASP', 'LEU', 'ASN', 'ALA', 'VAL', 'ALA', 'PHE', 'ARG', 'PHE', 'LYS', 'ASP', 'LYS', 'GLU', 'ARG', 'PRO', 'ARG', 'GLY', 'TYR', 'THR', 'SER', 'LYS', 'ILE', 'ALA', 'GLY', 'ILE', 'LEU', 'HSD', 'TYR', 'TYR', 'PRO', 'LEU', 'GLU', 'GLU', 'VAL', 'LEU', 'THR', 'ALA', 'GLU', 'TYR', 'LEU', 'LEU', 'GLU', 'GLU', 'PHE', 'ARG', 'PRO', 'ASP', 'LEU', 'ILE', 'GLU', 'MET', 'VAL', 'LEU', 'ASP', 'LYS', 'LEU', 'ARG', 'PRO', 'GLU', 'ASN', 'VAL', 'ARG', 'VAL', 'ALA', 'ILE', 'VAL', 'SER', 'LYS', 'SER', 'PHE', 'GLU', 'GLY', 'LYS', 'THR', 'ASP', 'ARG', 'THR', 'GLU', 'GLU', 'TRP', 'TYR', 'GLY', 'THR', 'GLN', 'TYR', 'LYS', 'GLN', 'GLU', 'ALA', 'ILE', 'PRO', 'ASP', 'GLU', 'VAL', 'ILE', 'LYS', 'LYS', 'TRP', 'GLN', 'ASN', 'ALA', 'ASP', 'LEU', 'ASN', 'GLY', 'LYS', 'PHE', 'LYS', 'LEU', 'PRO', 'THR', 'LYS', 'ASN', 'GLU', 'PHE', 'ILE', 'PRO', 'THR', 'ASN', 'PHE', 'GLU', 'ILE', 'LEU', 'PRO', 'LEU', 'GLU', 'LYS', 'GLU', 'ALA', 'THR', 'PRO', 'TYR', 'PRO', 'ALA', 'LEU', 'ILE', 'LYS', 'ASP', 'THR', 'ALA', 'MET', 'SER', 'LYS', 'LEU', 'TRP', 'PHE', 'LYS', 'GLN', 'ASP', 'ASP', 'LYS', 'PHE', 'PHE', 'LEU', 'PRO', 'LYS', 'ALA', 'ASN', 'LEU', 'ASN', 'PHE', 'GLU', 'PHE', 'PHE', 'SER', 'PRO', 'PHE', 'ALA', 'TYR', 'VAL', 'ASP', 'PRO', 'LEU', 'HSD', 'SER', 'ASN', 'MET', 'ALA', 'TYR', 'LEU', 'TYR', 'LEU', 'GLU', 'LEU', 'LEU', 'LYS', 'ASP', 'SER', 'LEU', 'ASN', 'GLU', 'TYR', 'ALA', 'TYR', 'ALA', 'ALA', 'GLU', 'LEU', 'ALA', 'GLY', 'LEU', 'SER', 'TYR', 'ASP', 'LEU', 'GLN', 'ASN', 'THR', 'ILE', 'TYR', 'GLY', 'MET', 'TYR', 'LEU', 'SER', 'VAL', 'LYS', 'GLY', 'TYR', 'ASN', 'ASP', 'LYS', 'GLN', 'PRO', 'ILE', 'LEU', 'LEU', 'LYS', 'LYS', 'ILE', 'ILE', 'GLU', 'LYS', 'MET', 'ALA', 'THR', 'PHE', 'GLU', 'ILE', 'ASP', 'GLU', 'LYS', 'ARG', 'PHE', 'GLU', 'ILE', 'ILE', 'LYS', 'GLU', 'ALA', 'TYR', 'MET', 'ARG', 'SER', 'LEU', 'ASN', 'ASN', 'PHE', 'ARG', 'ALA', 'GLU', 'GLN', 'PRO', 'HSD', 'GLN', 'HSD', 'ALA', 'MET', 'TYR', 'TYR', 'LEU', 'ARG', 'LEU', 'LEU', 'MET', 'THR', 'GLU', 'VAL', 'ALA', 'TRP', 'THR', 'LYS', 'ASP', 'GLU', 'LEU', 'LYS', 'GLU', 'ALA', 'LEU', 'ASP', 'ASP', 'VAL', 'THR', 'LEU', 'PRO', 'ARG', 'LEU', 'LYS', 'ALA', 'PHE', 'ILE', 'PRO', 'GLN', 'LEU', 'LEU', 'SER', 'ARG', 'LEU', 'HSD', 'ILE', 'GLU', 'ALA', 'LEU', 'LEU', 'HSD', 'GLY', 'ASN', 'ILE', 'THR', 'LYS', 'GLN', 'ALA', 'ALA', 'LEU', 'GLY', 'ILE', 'MET', 'GLN', 'MET', 'VAL', 'GLU', 'ASP', 'THR', 'LEU', 'ILE', 'GLU', 'HSD', 'ALA', 'HSD', 'THR', 'LYS', 'PRO', 'LEU', 'LEU', 'PRO', 'SER', 'GLN', 'LEU', 'VAL', 'ARG', 'TYR', 'ARG', 'GLU', 'VAL', 'GLN', 'LEU', 'PRO', 'ASP', 'ARG', 'GLY', 'TRP', 'PHE', 'VAL', 'TYR', 'GLN', 'GLN', 'ARG', 'ASN', 'GLU', 'VAL', 'HSD', 'ASN', 'ASN', 'SER', 'GLY', 'ILE', 'GLU', 'ILE', 'TYR', 'TYR', 'GLN', 'THR', 'ASP', 'MET', 'GLN', 'SER', 'THR', 'SER', 'GLU', 'ASN', 'MET', 'PHE', 'LEU', 'GLU', 'LEU', 'PHE', 'ALA', 'GLN', 'ILE', 'ILE', 'SER', 'GLU', 'PRO', 'ALA', 'PHE', 'ASN', 'THR', 'LEU', 'ARG', 'THR', 'LYS', 'GLU', 'GLN', 'LEU', 'GLY', 'TYR', 'ILE', 'VAL', 'PHE', 'SER', 'GLY', 'PRO', 'ARG', 'ARG', 'ALA', 'ASN', 'GLY', 'ILE', 'GLN', 'GLY', 'LEU', 'ARG', 'PHE', 'ILE', 'ILE', 'GLN', 'SER', 'GLU', 'LYS', 'PRO', 'PRO', 'HSD', 'TYR', 'LEU', 'GLU', 'SER', 'ARG', 'VAL', 'GLU', 'ALA', 'PHE', 'LEU', 'ILE', 'THR', 'MET', 'GLU', 'LYS', 'SER', 'ILE', 'GLU', 'ASP', 'MET', 'THR', 'GLU', 'GLU', 'ALA', 'PHE', 'GLN', 'LYS', 'HSD', 'ILE', 'GLN', 'ALA', 'LEU', 'ALA', 'ILE', 'ARG', 'ARG', 'LEU', 'ASP', 'LYS', 'PRO', 'LYS', 'LYS', 'LEU', 'SER', 'ALA', 'GLU', 'SER', 'ALA', 'LYS', 'TYR', 'TRP', 'GLY', 'GLU', 'ILE', 'ILE', 'SER', 'GLN', 'GLN', 'TYR', 'ASN', 'PHE', 'ASP', 'ARG', 'ASP', 'ASN', 'THR', 'GLU', 'VAL', 'ALA', 'TYR', 'LEU', 'LYS', 'THR', 'LEU', 'THR', 'LYS', 'GLU', 'ASP', 'ILE', 'ILE', 'LYS', 'PHE', 'TYR', 'LYS', 'GLU', 'MET', 'LEU', 'ALA', 'VAL', 'ASP', 'ALA', 'PRO', 'ARG', 'ARG', 'HSD', 'LYS', 'VAL', 'SER', 'VAL', 'HSD', 'VAL', 'LEU', 'ALA', 'ARG', 'GLU', 'MET', 'ASP', 'SER', 'ASN', 'PRO', 'VAL', 'VAL', 'GLY', 'GLU', 'PHE', 'PRO', 'ALA', 'GLN', 'ASN', 'ASP', 'ILE', 'ASN', 'LEU', 'SER', 'GLN', 'ALA', 'PRO', 'ALA', 'LEU', 'PRO', 'GLN', 'PRO', 'GLU', 'VAL', 'ILE', 'GLN', 'ASN', 'MET', 'THR', 'GLU', 'PHE', 'LYS', 'ARG', 'GLY', 'LEU', 'PRO', 'LEU', 'PHE', 'PRO', 'LEU', 'VAL', 'LYS', 'PRO', 'HSD', 'ILE', 'ASN', 'PHE', 'MET', 'ALA', 'ALA', 'LYS', 'LEU')
ABETA_LONG_STRING = "".join(ABETA_LONG_SEQUENCE)
ABETA_LONG_RESIDUE_COUNT = len(ABETA_LONG_SEQUENCE)

ABETA_SHORT_SEQUENCE = ('ASN', 'ASN', 'PRO', 'ALA', 'ILE', 'LYS', 'ARG', 'ILE', 'GLY', 'ASN', 'HSD', 'ILE', 'THR', 'LYS', 'SER', 'PRO', 'GLU', 'ASP', 'LYS', 'ARG', 'GLU', 'TYR', 'ARG', 'GLY', 'LEU', 'GLU', 'LEU', 'ALA', 'ASN', 'GLY', 'ILE', 'LYS', 'VAL', 'LEU', 'LEU', 'ILE', 'SER', 'ASP', 'PRO', 'THR', 'THR', 'ASP', 'LYS', 'SER', 'SER', 'ALA', 'ALA', 'LEU', 'ASP', 'VAL', 'HSD', 'ILE', 'GLY', 'SER', 'LEU', 'SER', 'ASP', 'PRO', 'PRO', 'ASN', 'ILE', 'ALA', 'GLY', 'LEU', 'SER', 'HSD', 'PHE', 'LEU', 'GLU', 'HSD', 'MET', 'LEU', 'PHE', 'LEU', 'GLY', 'THR', 'LYS', 'LYS', 'TYR', 'PRO', 'LYS', 'GLU', 'ASN', 'GLU', 'TYR', 'SER', 'GLN', 'PHE', 'LEU', 'SER', 'GLU', 'HSD', 'ALA', 'GLY', 'SER', 'SER', 'ASN', 'ALA', 'PHE', 'THR', 'SER', 'GLY', 'GLU', 'HSD', 'THR', 'ASN', 'TYR', 'TYR', 'PHE', 'ASP', 'VAL', 'SER', 'HSD', 'GLU', 'HSD', 'LEU', 'GLU', 'GLY', 'ALA', 'LEU', 'ASP', 'ARG', 'PHE', 'ALA', 'GLN', 'PHE', 'PHE', 'LEU', 'SER', 'PRO', 'LEU', 'PHE', 'ASP', 'GLU', 'SER', 'ALA', 'LYS', 'ASP', 'ARG', 'GLU', 'VAL', 'ASN', 'ALA', 'VAL', 'ASP', 'SER', 'GLU', 'HSD', 'GLU', 'LYS', 'ASN', 'VAL', 'MET', 'ASN', 'ASP', 'ALA', 'TRP', 'ARG', 'LEU', 'PHE', 'GLN', 'LEU', 'GLU', 'LYS', 'ALA', 'THR', 'GLY', 'ASN', 'PRO', 'LYS', 'HSD', 'PRO', 'PHE', 'SER', 'LYS', 'PHE', 'GLY', 'THR', 'GLY', 'ASN', 'LYS', 'TYR', 'THR', 'LEU', 'GLU', 'THR', 'ARG', 'PRO', 'ASN', 'GLN', 'GLU', 'GLY', 'ILE', 'ASP', 'VAL', 'ARG', 'GLN', 'GLU', 'LEU', 'LEU', 'LYS', 'PHE', 'HSD', 'SER', 'ALA', 'TYR', 'TYR', 'SER', 'SER', 'ASN', 'LEU', 'MET', 'ALA', 'VAL', 'VAL', 'VAL', 'LEU', 'GLY', 'ARG', 'GLU', 'SER', 'LEU', 'ASP', 'ASP', 'LEU', 'THR', 'ASN', 'LEU', 'VAL', 'VAL', 'LYS', 'LEU', 'PHE', 'SER', 'GLU', 'VAL', 'GLU', 'ASN', 'LYS', 'ASN', 'VAL', 'PRO', 'LEU', 'PRO', 'GLU', 'PHE', 'PRO', 'GLU', 'HSD', 'PRO', 'PHE', 'GLN', 'GLU', 'GLU', 'HSD', 'LEU', 'LYS', 'GLN', 'LEU', 'TYR', 'LYS', 'ILE', 'VAL', 'PRO', 'ILE', 'LYS', 'ASP', 'ILE', 'ARG', 'ASN', 'LEU', 'TYR', 'VAL', 'THR', 'PHE', 'PRO', 'ILE', 'PRO', 'ASP', 'LEU', 'GLN', 'LYS', 'TYR', 'TYR', 'LYS', 'SER', 'ASN', 'PRO', 'GLY', 'HSD', 'TYR', 'LEU', 'GLY', 'HSD', 'LEU', 'ILE', 'GLY', 'HSD', 'GLU', 'GLY', 'PRO', 'GLY', 'SER', 'LEU', 'LEU', 'SER', 'GLU', 'LEU', 'LYS', 'SER', 'LYS', 'GLY', 'TRP', 'VAL', 'ASN', 'THR', 'LEU', 'VAL', 'GLY', 'GLY', 'GLN', 'LYS', 'GLU', 'GLY', 'ALA', 'ARG', 'GLY', 'PHE', 'MET', 'PHE', 'PHE', 'ILE', 'ILE', 'ASN', 'VAL', 'ASP', 'LEU', 'THR', 'GLU', 'GLU', 'GLY', 'LEU', 'LEU', 'HSD', 'VAL', 'GLU', 'ASP', 'ILE', 'ILE', 'LEU', 'HSD', 'MET', 'PHE', 'GLN', 'TYR', 'ILE', 'GLN', 'LYS', 'LEU', 'ARG', 'ALA', 'GLU', 'GLY', 'PRO', 'GLN', 'GLU', 'TRP', 'VAL', 'PHE', 'GLN', 'GLU', 'LEU', 'LYS', 'ASP', 'LEU', 'ASN', 'ALA', 'VAL', 'ALA', 'PHE', 'ARG', 'PHE', 'LYS', 'ASP', 'LYS', 'GLU', 'ARG', 'PRO', 'ARG', 'GLY', 'TYR', 'THR', 'SER', 'LYS', 'ILE', 'ALA', 'GLY', 'ILE', 'LEU', 'HSD', 'TYR', 'TYR', 'PRO', 'LEU', 'GLU', 'GLU', 'VAL', 'LEU', 'THR', 'ALA', 'GLU', 'TYR', 'LEU', 'LEU', 'GLU', 'GLU', 'PHE', 'ARG', 'PRO', 'ASP', 'LEU', 'ILE', 'GLU', 'MET', 'VAL', 'LEU', 'ASP', 'LYS', 'LEU', 'ARG', 'PRO', 'GLU', 'ASN', 'VAL', 'ARG', 'VAL', 'ALA', 'ILE', 'VAL', 'SER', 'LYS', 'SER', 'PHE', 'GLU', 'GLY', 'LYS', 'THR', 'ASP', 'ARG', 'THR', 'GLU', 'GLU', 'TRP', 'TYR', 'GLY', 'THR', 'GLN', 'TYR', 'LYS', 'GLN', 'GLU', 'ALA', 'ILE', 'PRO', 'ASP', 'GLU', 'VAL', 'ILE', 'LYS', 'LYS', 'TRP', 'GLN', 'ASN', 'ALA', 'ASP', 'LEU', 'ASN', 'GLY', 'LYS', 'PHE', 'LYS', 'LEU', 'PRO', 'THR', 'LYS', 'ASN', 'GLU', 'PHE', 'ILE', 'PRO', 'THR', 'ASN', 'PHE', 'GLU', 'ILE', 'LEU', 'PRO', 'LEU', 'GLU', 'LYS', 'GLU', 'ALA', 'THR', 'PRO', 'TYR', 'PRO', 'ALA', 'LEU', 'ILE', 'LYS', 'ASP', 'THR', 'ALA', 'MET', 'SER', 'LYS', 'LEU', 'TRP', 'PHE', 'LYS', 'GLN', 'ASP', 'ASP', 'LYS', 'PHE', 'PHE', 'LEU', 'PRO', 'LYS', 'ALA', 'ASN', 'LEU', 'ASN', 'PHE', 'GLU', 'PHE', 'PHE', 'SER', 'PRO', 'PHE', 'ALA', 'TYR', 'VAL', 'ASP', 'PRO', 'LEU', 'HSD', 'SER', 'ASN', 'MET', 'ALA', 'TYR', 'LEU', 'TYR', 'LEU', 'GLU', 'LEU', 'LEU', 'LYS', 'ASP', 'SER', 'LEU', 'ASN', 'GLU', 'TYR', 'ALA', 'TYR', 'ALA', 'ALA', 'GLU', 'LEU', 'ALA', 'GLY', 'LEU', 'SER', 'TYR', 'ASP', 'LEU', 'GLN', 'ASN', 'THR', 'ILE', 'TYR', 'GLY', 'MET', 'TYR', 'LEU', 'SER', 'VAL', 'LYS', 'GLY', 'TYR', 'ASN', 'ASP', 'LYS', 'GLN', 'PRO', 'ILE', 'LEU', 'LEU', 'LYS', 'LYS', 'ILE', 'ILE', 'GLU', 'LYS', 'MET', 'ALA', 'THR', 'PHE', 'GLU', 'ILE', 'ASP', 'GLU', 'LYS', 'ARG', 'PHE', 'GLU', 'ILE', 'ILE', 'LYS', 'GLU', 'ALA', 'TYR', 'MET', 'ARG', 'SER', 'LEU', 'ASN', 'ASN', 'PHE', 'ARG', 'ALA', 'GLU', 'GLN', 'PRO', 'HSD', 'GLN', 'HSD', 'ALA', 'MET', 'TYR', 'TYR', 'LEU', 'ARG', 'LEU', 'LEU', 'MET', 'THR', 'GLU', 'VAL', 'ALA', 'TRP', 'THR', 'LYS', 'ASP', 'GLU', 'LEU', 'LYS', 'GLU', 'ALA', 'LEU', 'ASP', 'ASP', 'VAL', 'THR', 'LEU', 'PRO', 'ARG', 'LEU', 'LYS', 'ALA', 'PHE', 'ILE', 'PRO', 'GLN', 'LEU', 'LEU', 'SER', 'ARG', 'LEU', 'HSD', 'ILE', 'GLU', 'ALA', 'LEU', 'LEU', 'HSD', 'GLY', 'ASN', 'ILE', 'THR', 'LYS', 'GLN', 'ALA', 'ALA', 'LEU', 'GLY', 'ILE', 'MET', 'GLN', 'MET', 'VAL', 'GLU', 'ASP', 'THR', 'LEU', 'ILE', 'GLU', 'HSD', 'ALA', 'HSD', 'THR', 'LYS', 'PRO', 'LEU', 'LEU', 'PRO', 'SER', 'GLN', 'LEU', 'VAL', 'ARG', 'TYR', 'ARG', 'GLU', 'VAL', 'GLN', 'LEU', 'PRO', 'ASP', 'ARG', 'GLY', 'TRP', 'PHE', 'VAL', 'TYR', 'GLN', 'GLN', 'ARG', 'ASN', 'GLU', 'VAL', 'HSD', 'ASN', 'ASN', 'SER', 'GLY', 'ILE', 'GLU', 'ILE', 'TYR', 'TYR', 'GLN', 'THR', 'ASP', 'MET', 'GLN', 'SER', 'THR', 'SER', 'GLU', 'ASN', 'MET', 'PHE', 'LEU', 'GLU', 'LEU', 'PHE', 'ALA', 'GLN', 'ILE', 'ILE', 'SER', 'GLU', 'PRO', 'ALA', 'PHE', 'ASN', 'THR', 'LEU', 'ARG', 'THR', 'LYS', 'GLU', 'GLN', 'LEU', 'GLY', 'TYR', 'ILE', 'VAL', 'PHE', 'SER', 'GLY', 'PRO', 'ARG', 'ARG', 'ALA', 'ASN', 'GLY', 'ILE', 'GLN', 'GLY', 'LEU', 'ARG', 'PHE', 'ILE', 'ILE', 'GLN', 'SER', 'GLU', 'LYS', 'PRO', 'PRO', 'HSD', 'TYR', 'LEU', 'GLU', 'SER', 'ARG', 'VAL', 'GLU', 'ALA', 'PHE', 'LEU', 'ILE', 'THR', 'MET', 'GLU', 'LYS', 'SER', 'ILE', 'GLU', 'ASP', 'MET', 'THR', 'GLU', 'GLU', 'ALA', 'PHE', 'GLN', 'LYS', 'HSD', 'ILE', 'GLN', 'ALA', 'LEU', 'ALA', 'ILE', 'ARG', 'ARG', 'LEU', 'ASP', 'LYS', 'PRO', 'LYS', 'LYS', 'LEU', 'SER', 'ALA', 'GLU', 'SER', 'ALA', 'LYS', 'TYR', 'TRP', 'GLY', 'GLU', 'ILE', 'ILE', 'SER', 'GLN', 'GLN', 'TYR', 'ASN', 'PHE', 'ASP', 'ARG', 'ASP', 'ASN', 'THR', 'GLU', 'VAL', 'ALA', 'TYR', 'LEU', 'LYS', 'THR', 'LEU', 'THR', 'LYS', 'GLU', 'ASP', 'ILE', 'ILE', 'LYS', 'PHE', 'TYR', 'LYS', 'GLU', 'MET', 'LEU', 'ALA', 'VAL', 'ASP', 'ALA', 'PRO', 'ARG', 'ARG', 'HSD', 'LYS', 'VAL', 'SER', 'VAL', 'HSD', 'VAL', 'LEU', 'ALA', 'ARG', 'GLU', 'MET', 'ASP', 'SER', 'ASN', 'PRO', 'VAL', 'VAL', 'GLY', 'GLU', 'PHE', 'PRO', 'ALA', 'GLN', 'ASN', 'ASP', 'ILE', 'ASN', 'LEU', 'SER', 'GLN', 'ALA', 'PRO', 'ALA', 'LEU', 'PRO', 'GLN', 'PRO', 'GLU', 'VAL', 'ILE', 'GLN', 'ASN', 'MET', 'THR', 'GLU', 'PHE', 'LYS', 'ARG', 'GLY', 'LEU', 'PRO', 'LEU', 'PHE', 'PRO', 'LEU', 'VAL', 'LYS', 'PRO', 'HSD', 'ILE', 'ASN', 'PHE', 'MET', 'ALA', 'ALA', 'LYS', 'LEU')
ABETA_SHORT_STRING = "".join(ABETA_SHORT_SEQUENCE)
ABETA_SHORT_RESIDUE_COUNT = len(ABETA_SHORT_SEQUENCE)
