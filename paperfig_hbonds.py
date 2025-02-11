import numpy as np
import pandas as pd
import seaborn as sns
from scipy.ndimage import convolve
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec


#Description
# ----------------------------------------------------------------------------------
# This version adjusts vmin and vmax so it looks more binary (like v2)
# This time, we want to make pretty axes (nice indexing in x-label),
# indicate regions where nucleosomes form and remap the y-axis according to this.

min_val = 4  # Any value below this will be 0 and above 1
#enhance_res = True
enhance_res = False

if enhance_res:
    # Convolution kernel to expand 1s to neighboring rows
    #kernel = np.array([[0, 1, 0], [0, 1, 0], [0, 1, 0]])
    kernel = np.array([[1], [1], [1]])

# System info
nbp = 358  # number of base-pairs in minicircle
n_length = 147  # nucleosome length
l_length = int((nbp - n_length*2)/2)  # linker length (arms)
nucleosome1_position = 17
nucleosome2_position = 196
nucleosome1_apex = 90
nucleosome2_apex = 269

# Correct apexes
n1_apex = nucleosome1_apex - nucleosome1_position
n2_apex = nucleosome2_apex - nucleosome1_position
apex_color = 'red'
apex_ls = '--'
apex_lw = 1

# Workout regions, labels, colors and their alphas
tick_positions = [0, n_length-1, n_length+l_length-1, n_length*2+l_length-1, nbp-1]  # Positions on the y-axis
tick_labels = ["NUC 1", "L1", "NUC 2", "L2"]  # Labels for the regions
colors = ["red", "blue", "green", "blue"]  # Colors for the regions

alpha0 = 0.2
alpha1 = 0.7
alpha_list = [
    [alpha0, alpha0, alpha0, alpha0],
    [alpha0, alpha0, alpha0, alpha0],
    [alpha0, alpha0, alpha0, alpha0],
    [alpha0, alpha0, alpha0, alpha0],
    [alpha0, alpha0, alpha1, alpha0],
#    [alpha0, alpha0, alpha1, alpha0],
    [alpha1, alpha0, alpha1, alpha0],
    [alpha1, alpha0, alpha1, alpha0],
]

#Figure parms
# ----------------------------------------------------------------------------------
width = 6
height = 4#4
cmap = 'plasma'
tick_s = 10
label_s = 14
title_s = 16

#yticks = np.linspace(0, nbp, 10, dtype=int)
yticks = np.arange(0, nbp+1, 50, dtype=int)

# And plot
#fig, axs = plt.subplots(4,2, figsize=(width*2, height*3), tight_layout=True)
fig = plt.figure(figsize=(width*2, height*4), tight_layout=True)
gs = gridspec.GridSpec(4, 2, figure=fig)  # 4 rows, 2 columns

# Create subplots for the first two rows
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])

# Create a single subplot spanning the third row
#ax5 = fig.add_subplot(gs[2, :])  # Span both columns
#ax5 = fig.add_subplot(gs[2, 0])  # Place it in the first column of the 3rd row
ax5 = fig.add_axes([0.25, 0.28, 0.5, 0.2])  # [left, bottom, width, height]

# Create subplots for the fourth row
ax6 = fig.add_subplot(gs[3, 0])
ax7 = fig.add_subplot(gs[3, 1])

# List of dicts with info
info_list = [
    {'name':"no_nuc_CA", "title":"Empty minicircle Ca r1", "ax": ax1},
    {'name': "no_nuc_CA_r1", "title": "Empty minicircle Ca r2", "ax": ax2},
    {'name': "no_nuc_Na", "title": "Empty minicircle Na r1", "ax": ax3},
    {'name': "no_nuc_Na_r1", "title": "Empty minicircle Na r2", "ax": ax4},
    {'name': "one_nuc_Na", "title": "Mononucleosome Ca", "ax": ax5},
    {'name': "two_nuc_CA", "title": "Dinucleosome Ca", "ax": ax6},
    {'name': "two_nuc_Na", "title": "Dinucleosome Na", "ax": ax7},
]

for i, info_dict in enumerate(info_list):
    ax = info_dict['ax']
    title = info_dict['title']
    name = info_dict['name']

    ax.set_title(title,fontsize=title_s)

    # Let's load the hbonds dfs
    hbond_df = pd.read_csv("hbonds_{}.csv".format(name))
    hbond_df = (hbond_df >= min_val).astype(int)

    if enhance_res:
        expanded_data = convolve(hbond_df.values, kernel, mode="constant", cval=0)
        expanded_data[expanded_data > 0] = 1
        hbond_df = pd.DataFrame(expanded_data, columns=hbond_df.columns)

    # Rearrange rows
    hbond_df = pd.concat([hbond_df.iloc[nucleosome1_position-1:], hbond_df.iloc[:nucleosome1_position-1]], ignore_index=True)

    # Plot
    sns.heatmap(hbond_df, ax=ax, cmap=cmap, cbar=False)

    # Draw nucleosome regions
    # -----------------------------------------------------------
    ax.set_yticks(tick_positions)  # Draw ticks for the nucleosome positions

    # Add labels at the middle of each region and draw a rectangle
    for k in range(len(tick_labels)):
        start = tick_positions[k]
        end = tick_positions[k+1]
        label = tick_labels[k]
        midpoint = (start + end) / 2  - 3
        alpha = alpha_list[i][k]
        color = colors[k]
        # Draw region label
        ax.text(
            -0.5, midpoint, label,  # -0.5 positions the label outside the heatmap
            ha="right", va="center", rotation=90, fontsize=tick_s
        )
        # Draw region with rectangle
        rect = patches.Rectangle(
            (-6., start),  # x, y: starting position (-1.5 places it outside the heatmap)
            6.,  # width
            end - start,  # height of the rectangle
            linewidth=0,  # No border
            edgecolor=None,
            facecolor=color,
            clip_on=False,
            alpha=alpha  # Transparency
        )
        ax.add_patch(rect)

    # Finally, draw apexes
    time = len(hbond_df.iloc[0])
    ax.plot([0, time], [n1_apex,n1_apex], lw=apex_lw, ls=apex_ls, color=apex_color, alpha=0.5)
    ax.plot([0, time], [n2_apex, n2_apex], lw=apex_lw, ls=apex_ls, color=apex_color, alpha=0.5)

    # SORT TICKS
    # -----------------------------------------------------------
    xticks = np.arange(0, time+1, 25, dtype=int)
#    ax.set_yticks([])
    #ax.set_yticks(yticks)
    #ax.set_yticklabels(yticks, fontsize=tick_s)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks, fontsize=tick_s)
    ax.set_xlabel('Time (ns)', fontsize=label_s)
    #ax.set_ylabel('DNA Base Pair', fontsize=label_s)
    # ax.set_ylabel(' ')
    ax.grid(True, alpha=0.5)


if enhance_res:
    plt.savefig("hbonds_paper_e.pdf")
    plt.savefig("hbonds_paper_e.png")
else:
    plt.savefig("hbonds_paper.pdf")
    plt.savefig("hbonds_paper.png")
plt.show()



