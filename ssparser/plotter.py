import numpy as np
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
import pandas as pd
import os

from .constants import *


def plot_muiltiple(template_dict: dict, output_file=None):

    # remove coils and turns
    template_exclude_unstrcutured = dict(template_dict)
    del template_exclude_unstrcutured["Other"]

    x = np.arange(ABETA_LONG_RESIDUE_COUNT)  
    width = 0.33  # the width of the bars
    multiplier = 3

    fig, ax = plt.subplots(layout='constrained')
    plt.ylim(0.0, 1.1)
    # plt.xticks(np.arange(1, ABETA_LONG_RESIDUE_COUNT+1, dtype=int), ABETA_LONG_SEQUENCE, rotation=-60)

    # add line at 50%
    # plt.axhline(y=0.5, color='r', linestyle='-')

    for attribute, measurement in template_exclude_unstrcutured.items():
        offset = width * multiplier
        rects = ax.bar(x + offset, measurement, width, label=attribute)
        multiplier += 1

    ax.set_ylabel("Frequency of Secondary Structure for")
    ax.set_xlabel("Residue")
    ax.set_title("Frequency of Secondary Structure per Residue")
    ax.legend(loc='upper right')

    if output_file:
        plt.savefig(output_file, dpi=200)
    else:
        plt.show()
    plt.close()

def plot_complex_hits(ss_sequence_list, template_dict, output_file=None):
    beta_color = "orange"
    turn_color = "green"
    helix_color = "blue"
    coil_color = "red"

    chain_count = len(ss_sequence_list)

    x = np.arange(ABETA_LONG_RESIDUE_COUNT, dtype=int)
    fig, ax = plt.subplots(2, layout='constrained')
    
    width = 0.7
    offset = 0

    beta_color_strength = [0 for _ in range(ABETA_LONG_RESIDUE_COUNT)]
    turn_color_strength = [0 for _ in range(ABETA_LONG_RESIDUE_COUNT)]

    for ss_sequence in ss_sequence_list:
        for i in range(0, ABETA_LONG_RESIDUE_COUNT):
            beta_color_strength[i] += (1 / chain_count) if ss_sequence[i] in ["B", "b", "E"] else 0
            turn_color_strength[i] += (1 / chain_count) if ss_sequence[i] in ["T"] else 0

    data_x = x
    fig, ax = plt.subplots(2, layout='constrained')

    b_cmap = plt.cm._colormaps['GnBu'] # GnBu
    t_cmap = plt.cm._colormaps['GnBu']
    b_colors = b_cmap(beta_color_strength)
    t_colors = t_cmap(turn_color_strength)

    ax[0].bar(data_x, template_dict["Beta"], color=b_colors)
    ax[1].bar(data_x, template_dict["Turn"], color=t_colors)

    sm = ScalarMappable(cmap=b_cmap, norm=plt.Normalize(0,1))
    sm.set_array([])

    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Freq of SS in AlphaFold Complex', rotation=270,labelpad=25)

    ax[0].set_ylim([0.0, 1.05])
    ax[1].set_ylim([0.0, 1.05])

    ax[0].set_ylabel("Freq of Beta SS in PDB")
    ax[1].set_ylabel("Freq of Turn SS in PDB")
    ax[1].set_xlabel("Residue")
    ax[0].set_title("Frequency of Secondary Structure per Residue of ABeta-Monomer")

    if output_file:
        plt.savefig(output_file, dpi=400)
        plt.close()
    else:
        plt.show()

    plt.close()
    


def plot_monomer_hits(ss_sequence: list, template_dict: dict, output_file=None):
    beta_color = "orange"
    turn_color = "green"
    helix_color = "blue"
    coil_color = "red"

    x = np.arange(ABETA_LONG_RESIDUE_COUNT, dtype=int)

    fig, ax = plt.subplots(2, layout='constrained')
    ax[0].set_xticks(x, ABETA_LONG_SEQUENCE, rotation=-45, fontsize = 6)
    ax[1].set_xticks(x, ABETA_LONG_SEQUENCE, rotation=-45, fontsize = 6)
    ax[0].set_ylim([0.0, 1.05])
    ax[1].set_ylim([0.0, 1.05])

    #plt.figure(figsize=(10, 5))
    plt.ylim(0.0, 1.0)

    beta_data = template_dict["Beta"]
    turn_data = template_dict["Turn"]
    
    width = 0.7
    offset = 0

    for i, measurement in enumerate(beta_data):
        alpha = 1 
        ax[0].bar(offset, measurement, width, color=beta_color, alpha=alpha)
        offset += 1
    
    offset = 0
    for i, measurement in enumerate(turn_data):
        alpha = 1 
        ax[1].bar(offset, measurement, width, color=turn_color, alpha=alpha)
        offset += 1

    ax[0].set_ylabel("Frequency of Beta SS in RCSB")
    ax[1].set_ylabel("Frequency of Turn SS  in RCSB")
    ax[1].set_xlabel("Residue")
    ax[0].set_title("Frequency of Secondary Structure per Residue of ABeta-Monomer")


    if output_file:
        plt.savefig(output_file, dpi=400)
        plt.close()
    else:
        plt.show()
    plt.close()
    