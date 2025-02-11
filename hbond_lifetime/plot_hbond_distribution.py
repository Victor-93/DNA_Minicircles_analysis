import numpy as np
import pandas as pd
import seaborn as sns
from scipy.ndimage import convolve
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec


#Description
# ----------------------------------------------------------------------------------
# Here, we want to load the hbond dataframes (csv's) and do an analysis of lifetime.
# Should we first start by average lifetime? Maybbe a histogram of lifetime and then just
# count the occurances

min_val = 4  # Any value below this will be 0 and above 1

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

# TODO: I also don't think I need this

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

# Functions
# ----------------------------------------------------------------------------------
# Function to calculate run lengths of 1s in a binary array
def calculate_run_lengths(row):
    runs = []
    current_run = 0
    for val in row:
        if val == 1:
            current_run += 1
        else:
            if current_run > 0:
                runs.append(current_run)
            current_run = 0
    if current_run > 0:
        runs.append(current_run)  # Append the last run if it exists
    return runs

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
    {'name':"no_nuc_CA", "title":"No Nucleosome Ca", "ax": ax1},
    {'name': "no_nuc_CA_r1", "title": "No Nucleosome Ca (r1)", "ax": ax2},
    {'name': "no_nuc_Na", "title": "No Nucleosome Na", "ax": ax3},
    {'name': "no_nuc_Na_r1", "title": "No Nucleosome Na (r1)", "ax": ax4},
    {'name': "one_nuc_Na", "title": "One Nucleosome Ca", "ax": ax5},
    {'name': "two_nuc_CA", "title": "Two Nucleosomes Ca", "ax": ax6},
    {'name': "two_nuc_Na", "title": "Two Nucleosomes Na", "ax": ax7},
]

for i, info_dict in enumerate(info_list):
    ax = info_dict['ax']
    title = info_dict['title']
    name = info_dict['name']

    ax.set_title(title,fontsize=title_s)

    # Let's load the hbonds dfs
    hbond_df = pd.read_csv("../high-res_hbonds_{}.csv".format(name))
    hbond_df = (hbond_df >= min_val).astype(int)

    # Collect run lengths for all rows
    all_run_lengths = []
    for row in hbond_df.values:
        all_run_lengths.extend(calculate_run_lengths(row))

    all_run_lengths = np.array(all_run_lengths)/10.

    # Plot
    print(max(all_run_lengths))
    # print(all_run_lengths)
    # ax.hist(all_run_lengths, bins=range(1, max(all_run_lengths) + 2), align='left', color='blue', edgecolor='black')
    ax.hist(all_run_lengths, bins=100, align='left', color='blue', edgecolor='black',range=(2,150))
    #ax.hist(all_run_lengths, align='left', color='blue', edgecolor='black')

    # SORT TICKS
    # -----------------------------------------------------------
    #xticks = np.arange(0, time+1, 25, dtype=int)
    #ax.set_xticks(xticks)
    #ax.set_xticklabels(xticks, fontsize=tick_s)
    ax.set_xlabel('Time (ns)', fontsize=label_s)
    #ax.set_ylabel('DNA Base Pair', fontsize=label_s)
    # ax.set_ylabel(' ')
    #ax.grid(True, alpha=0.5)


#plt.savefig("hbonds_binary3.pdf")
#plt.savefig("hbonds_binary3.png")
plt.show()



