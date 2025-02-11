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
# This time just one histogram using all cases

min_val = 4  # Any value below this will be 0 and above 1

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

fig, axs = plt.subplots(1, figsize=(width, height), tight_layout=True)
ax = axs

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

ax.set_title('Hbond lifetime', fontsize=title_s)

all_run_lengths = []

# Collect run lengths for all rows
for i, info_dict in enumerate(info_list):
    title = info_dict['title']
    name = info_dict['name']

    # Let's load the hbonds dfs
    hbond_df = pd.read_csv("../high-res_hbonds_{}.csv".format(name))
    hbond_df = (hbond_df >= min_val).astype(int)

    for row in hbond_df.values:
        all_run_lengths.extend(calculate_run_lengths(row))

all_run_lengths = np.array(all_run_lengths)/10.
# Plot
print(max(all_run_lengths))
# print(all_run_lengths)
# ax.hist(all_run_lengths, bins=range(1, max(all_run_lengths) + 2), align='left', color='blue', edgecolor='black')
ax.hist(all_run_lengths, bins=100, align='left', color='blue', edgecolor='black',range=(1,101))
#ax.hist(all_run_lengths, align='left', color='blue', edgecolor='black')

# SORT TICKS
# -----------------------------------------------------------
xticks = np.arange(0, 101, 10, dtype=int)
yticks = np.arange(0, 19, 3, dtype=int)
ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_xlabel('Life time (ns)', fontsize=label_s)
ax.set_ylabel('Counts', fontsize=label_s)
#ax.grid(True)

#plt.savefig("hbonds_binary3.pdf")
#plt.savefig("hbonds_binary3.png")
plt.show()



