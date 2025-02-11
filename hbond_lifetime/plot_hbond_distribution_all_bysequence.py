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
# This time just one histogram using all cases and separating by sequence AT or GC

min_val = 4  # Any value below this will be 0 and above 1

# Load sequence
# ----------------------------------------------------------------------------------
# ChatGPT filtered it (don't judge me, it was the fastest way and I was lazy).
file_path = '../../dna_2hbonds_2dist.txt'
with open(file_path, 'r') as file:
    lines = file.readlines()

# Extract unique indices and their corresponding base1
unique_bases = {}
for line in lines:
    parts = line.split()
    index1, base1 = parts[0], parts[1]
    if index1 not in unique_bases:
        unique_bases[index1] = base1

# Convert the dictionary values to a list for the sequence
sequence = list(unique_bases.values())

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
width = 7
height = 4#4
cmap = 'plasma'
tick_s = 10
label_s = 14
title_s = 16

fig, axs = plt.subplots(3, figsize=(width, 3*height), tight_layout=True)

# List of dicts with info
info_list = [
    {'name':"no_nuc_CA", "title":"No Nucleosome Ca"},
    {'name': "no_nuc_CA_r1", "title": "No Nucleosome Ca (r1)"},
    {'name': "no_nuc_Na", "title": "No Nucleosome Na"},
    {'name': "no_nuc_Na_r1", "title": "No Nucleosome Na (r1)"},
    {'name': "one_nuc_Na", "title": "One Nucleosome Ca"},
    {'name': "two_nuc_CA", "title": "Two Nucleosomes Ca"},
    {'name': "two_nuc_Na", "title": "Two Nucleosomes Na"}
]

all_run_lengths = []
AT_runs = []
GC_runs = []
# Collect run lengths for all rows
for i, info_dict in enumerate(info_list):
    title = info_dict['title']
    name = info_dict['name']

    # Let's load the hbonds dfs
    hbond_df = pd.read_csv("../high-res_hbonds_{}.csv".format(name))
    hbond_df = (hbond_df >= min_val).astype(int)

    for s, row in enumerate(hbond_df.values):
        seq = sequence[s]
        if seq == 'DT' or seq == 'DA':
            AT_runs.extend(calculate_run_lengths(row))
        else:
            GC_runs.extend(calculate_run_lengths(row))
        all_run_lengths.extend(calculate_run_lengths(row)) # And all

all_run_lengths = np.array(all_run_lengths)/10.
AT_runs = np.array(AT_runs)/10.
GC_runs = np.array(GC_runs)/10.

# First all
ax = axs[0]
ax.set_title('All Base-Pairs', fontsize=title_s)
ax.hist(all_run_lengths, bins=100, align='left', color='blue', edgecolor='black',range=(1,101))

#AT
ax = axs[1]
ax.set_title('AT Base-Pairs', fontsize=title_s)
ax.hist(AT_runs, bins=100, align='left', color='red', edgecolor='black',range=(1,101))

#AT
ax = axs[2]
ax.set_title('GC Base-Pairs', fontsize=title_s)
ax.hist(GC_runs, bins=100, align='left', color='green', edgecolor='black',range=(1,101))

# SORT TICKS
# -----------------------------------------------------------
xticks = np.arange(0, 101, 10, dtype=int)
yticks = np.arange(0, 19, 3, dtype=int)
for ax in [axs[0], axs[1], axs[2]]:
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xlabel('Life time (ns)', fontsize=label_s)
    ax.set_ylabel('Counts', fontsize=label_s)
    ax.grid(True)

plt.savefig("hbond_lifetime1.pdf")
plt.savefig("hbond_lifetime1.png")
plt.show()



