import pandas as pd
import numpy as np
import os

from .stride import parse_with_stride
from .constants import *

def parse_ss_freq(struct_files: list[str], filetype="csv") -> dict:
    """
    parse_ss_freq: takes a directory of parsed .pdb/.ent files as .csv files
    and returns the frequency of each secondary strcuture in the abeta-chain

    input:
        struct_files: full paths
    
    output:
        A dictionary containing an array of the normalized number of that secondary
        structure at the given index. 
        in the form {"Alpha": alpha, "Beta": beta, "Turn": turn, "Other": other}
    """

    if filetype not in ["csv", "ent", "pdb"]:
        raise BaseException("parse_ss_freq: filetype must be one of csv ent psb")

    total_chain_count = 0
    long_chain_count  = 0

    alpha = np.zeros(ABETA_LONG_RESIDUE_COUNT, dtype=np.float32)
    beta  = np.zeros(ABETA_LONG_RESIDUE_COUNT, dtype=np.float32)
    turn  = np.zeros(ABETA_LONG_RESIDUE_COUNT, dtype=np.float32)
    other = np.zeros(ABETA_LONG_RESIDUE_COUNT, dtype=np.float32)

    for infile in struct_files:
        df = pd.read_csv(infile) if filetype == "csv" else parse_with_stride(infile)

        # iterate though each chain in the complex
        chain_names = set(df["Chain"])
        for cname in chain_names:
            chain = df.loc[df['Chain'] == cname]
            
            # check if the current chain is an ABeta monomer by checking if it is 
            # a substring of the total chain. Sequence the strcuture if it is.
            aa_sequence = list(chain["Residue Name"])
            if "".join(aa_sequence) in ABETA_LONG_STRING:
                
                chain_start = int(min(chain["PDB ResID"]))
                chain_stop  = int(max(chain["PDB ResID"]))

                total_chain_count += 1
                if chain_stop > 40:
                    long_chain_count += 1

                # iterate through each residue in the chain. If it is in the structure,
                # add it to out sequence
                for res_index in range(min(ABETA_LONG_RESIDUE_COUNT, chain_stop)):
                    # select data by residueID. Note that this is 1-indexed
                    resid = res_index + 1
                    residue_data = chain.loc[chain["PDB ResID"] == resid]

                    if not residue_data.empty:
                        structure = list(residue_data["SS Code"])[0]
                        
                        # alpha-helix
                        if structure in ["H", "G", "I"]:
                            alpha[res_index] += 1
                        
                        # beta strand or sheet
                        elif structure in ["B", "b", "E"]:
                            beta[res_index] += 1
                        
                        # turn
                        elif structure in ["T"]:
                            turn[res_index] += 1

                        # other, listed as a coil
                        else:
                            other[res_index] += 1
                    
                    else:
                        other[res_index] += 1
                
                # fill in blank residues at the tail with "other"
                for res_index in range(chain_stop, ABETA_SHORT_RESIDUE_COUNT):
                    other[res_index] += 1

    # normalize by dividing the quantity at each residue by the total number of chains
    normalizer = np.repeat(total_chain_count, ABETA_LONG_RESIDUE_COUNT)
    normalizer[-1] = long_chain_count if long_chain_count != 0 else 1
    normalizer[-2] = long_chain_count if long_chain_count != 0 else 1

    freq_dict = {"Alpha": alpha, "Beta": beta, "Turn": turn, "Other": other}
    for k, v in freq_dict.items():
        freq_dict[k] = np.divide(v, normalizer)

    return freq_dict
