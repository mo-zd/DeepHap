import os

import matplotlib.pyplot as plt

from read_write import *
from hadi import *

color1 = "#4f80c9"
color2 = "#e36f6f"
color3 = "#f70601"

def plot(denoiser, chr, single_parent=False):
    denoised = denoiser.denoised

    if not single_parent:
        chr_index = denoiser.order_dict[chr]
        fig, axes = plt.subplots(nrows=3, ncols=2)
        fig.tight_layout()
        col_names = dict([(el, denoised[el].columns[-2]) for el in ('M1', 'M2', 'P1', 'P2', 'Logr')])
        denoised_chr = dict([(el, denoiser.split(denoised[el])[chr_index]) for el in ('M1', 'M2', 'P1', 'P2', 'Logr')])
        pos_chr = dict([(el, denoised_chr[el]['Position'].tolist()) for el in ('M1', 'M2', 'P1', 'P2', 'Logr')])
        val_chr = dict([(el, denoised_chr[el][col_names[el]].tolist()) for el in ('M1', 'M2', 'P1', 'P2', 'Logr')])
        val_denoised_chr = dict([(el, denoised_chr[el]["denoised"].tolist()) for el in ('M1', 'M2', 'P1', 'P2', 'Logr')])
        
        baf_ticks_label = ["0.0", "0.33", "0.50", "0.66", "1.0"]
        baf_ticks_pos = [0.0, 1/3., 1/2., 2/3., 1.]
        ax = axes[0, 0]
        pos_1 = pos_chr['P1']
        pos_2 = pos_chr['P2']
        val_color1 = val_denoised_chr["P1"]
        val_color2 = val_denoised_chr["P2"]
        ax.scatter(pos_1, val_color1, c=color1, s=0.2)
        ax.scatter(pos_2, val_color2, c=color2, s=0.2)
        ax.set_ylabel('Paternal')
        ax.set_yticks(baf_ticks_pos)
        ax.set_yticklabels(baf_ticks_label)
        ax.grid(True, alpha=0.3)

        ax = axes[0, 1]
        val_color1 = val_chr["P1"]
        val_color2 = val_chr["P2"]
        ax.scatter(pos_1, val_color1, c=color1, s=0.2)
        ax.scatter(pos_2, val_color2, c=color2, s=0.2)
        ax.set_ylabel('Paternal (noisy)')
        ax.set_yticks(baf_ticks_pos)
        ax.set_yticklabels(baf_ticks_label)
        ax.grid(True, alpha=0.3)

        ax = axes[1, 0]
        pos_1 = pos_chr['M1']
        pos_2 = pos_chr['M2']
        val_color1 = val_denoised_chr["M1"]
        val_color2 = val_denoised_chr["M2"]
        ax.scatter(pos_1, val_color1, c=color1, s=0.2)
        ax.scatter(pos_2, val_color2, c=color2, s=0.2)
        ax.set_ylabel('Maternal')
        ax.set_yticks(baf_ticks_pos)
        ax.set_yticklabels(baf_ticks_label)
        ax.grid(True, alpha=0.3)

        ax = axes[1, 1]
        val_color1 = val_chr["M1"]
        val_color2 = val_chr["M2"]
        ax.scatter(pos_1, val_color1, c=color1, s=0.2)
        ax.scatter(pos_2, val_color2, c=color2, s=0.2)
        ax.set_ylabel('Matrnal (noisy)')
        ax.set_yticks(baf_ticks_pos)
        ax.set_yticklabels(baf_ticks_label)
        ax.grid(True, alpha=0.3)

        ax = axes[2, 0]
        pos = pos_chr['Logr']
        val_color = val_denoised_chr["Logr"]
        ax.scatter(pos, val_color, c=color3, s=0.2)
        ax.set_ylabel('LogR')
        ax.set_yticks([0, 1, 2, 3, 4])
        ax.grid(True, alpha=0.3)

        ax = axes[2, 1]
        val_color = np.array(val_chr["Logr"]) + 2
        ax.scatter(pos, val_color, c=color3, s=0.2)
        ax.set_ylabel('LogR (noisy)')
        ax.set_yticks([0, 1, 2, 3, 4])
        ax.grid(True, alpha=0.3)
