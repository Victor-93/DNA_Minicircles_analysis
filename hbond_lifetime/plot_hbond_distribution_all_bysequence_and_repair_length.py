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
# We want to calculate if they repaired themselves as well.

min_val = 4  # Any value below this will be 0 and above 1

repaired_colour = 'gray'
not_repaired_colour = 'yellow'

outfile = "hbond_lifetime_repairs"
showlabels = True # Indicate if you want to print labels

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

# Very similar to calculate_run_lengths, but here we only count events that repaired themselves.
def calculate_run_repairs(row):
    repaired = []  # The ones that repaired themselves
    not_repaired = []  # The ones that didn't repair
    current_run = 0
    for val in row:
        if val == 1:
            current_run += 1
        else:
            if current_run > 0:
                repaired.append(current_run)
            current_run = 0
    if current_run > 0:  # So, basically, only the disruptions at the end of simulations are the ones that didn't repair themselves.
        not_repaired.append(current_run)  # Append the last run if it exists
    return repaired, not_repaired

#Figure parms
# ----------------------------------------------------------------------------------
width = 6.5
height = 3.5
cmap = 'plasma'
tick_s = 12
label_s = 16
title_s = 18
legend_s=9

fig, axs = plt.subplots(2, figsize=(width, 2*height), tight_layout=True)

# List of dicts with info
info_list = [
    {'name':"no_nuc_CA", "title":"No Nucleosome Ca", 'length':1500},
    {'name': "no_nuc_CA_r1", "title": "No Nucleosome Ca (r1)", 'length':1750},
    {'name': "no_nuc_Na", "title": "No Nucleosome Na", 'length':1500},
    {'name': "no_nuc_Na_r1", "title": "No Nucleosome Na (r1)", 'length':1750},
    {'name': "one_nuc_Na", "title": "One Nucleosome Ca", 'length':1000},
    {'name': "two_nuc_CA", "title": "Two Nucleosomes Ca", 'length':1000},
    {'name': "two_nuc_Na", "title": "Two Nucleosomes Na", 'length':1000}
]

all_run_lengths = []
all_repaired = []
all_not_repaired = []
AT_runs = []
AT_repaired = []
AT_not_repaired = []
GC_runs = []
GC_repaired  = []
GC_not_repaired = []
# Collect run lengths for all rows
for i, info_dict in enumerate(info_list):
    title = info_dict['title']
    name = info_dict['name']
    l = info_dict['length']

    # Let's load the hbonds dfs
    hbond_df = pd.read_csv("../high-res_hbonds_{}.csv".format(name),index_col=0)

    # Filter columns: keep only those with names <= l -> Column names are the frames
    hbond_df = hbond_df.loc[:, hbond_df.columns.astype(float) <= l]
    hbond_df = (hbond_df >= min_val).astype(int)

    for s, row in enumerate(hbond_df.values):
        seq = sequence[s]
        if seq == 'DT' or seq == 'DA':
            AT_runs.extend(calculate_run_lengths(row))
            repaired, not_repaired = calculate_run_repairs(row)
            AT_repaired.extend(repaired)
            AT_not_repaired.extend(not_repaired)
            # if 508 in repaired:
            #    print(name)
        else:
            GC_runs.extend(calculate_run_lengths(row))
            repaired, not_repaired = calculate_run_repairs(row)
            GC_repaired.extend(repaired)
            GC_not_repaired.extend(not_repaired)

        all_run_lengths.extend(calculate_run_lengths(row)) # And all
        repaired, not_repaired = calculate_run_repairs(row)
        all_repaired.extend(repaired)
        all_not_repaired.extend(not_repaired)

# This was an idea for mixing arrays but didn't work
#first = all_repaired.copy()
#first.extend(all_not_repaired)
#second = all_repaired
#first = np.array(first)/10.
#second = np.array(second)/10.
#ax.hist(first, bins=100, align='left', color='yellow', edgecolor='black',range=(1,101))
#ax.hist(second, bins=100, align='left', color='blue', edgecolor='black',range=(1,101))

all_run_lengths = np.array(all_run_lengths)/10.
all_repaired = np.array(all_repaired )/10.
all_not_repaired = np.array(all_not_repaired)/10.
AT_runs = np.array(AT_runs)/10.
AT_repaired = np.array(AT_repaired)/10.
AT_not_repaired = np.array(AT_not_repaired)/10.
GC_runs = np.array(GC_runs)/10.
GC_repaired = np.array(GC_repaired)/10.
GC_not_repaired = np.array(GC_not_repaired)/10.

# Print the repaired values
print(np.array(list(set(AT_repaired))))

# And print how many events
print("All events ", len(all_run_lengths))
print("All AT events ", len(AT_runs))
print("All AT repairs ", len(AT_repaired))
print("All AT not repaires ", len(AT_not_repaired))
print("All GC events ", len(GC_runs))
print("All GC repairs ", len(GC_repaired))
print("All GC not repaires ", len(GC_not_repaired))

# First all
#ax = axs[0]
#ax.set_title('All Base-Pairs', fontsize=title_s)
#ax.hist(all_run_lengths, bins=100, align='left', color='blue', edgecolor='black',range=(1,101))
#ax.hist(all_repaired, bins=100, align='left', color='blue', edgecolor='black',range=(1,101))
# Second histogram (stacked on top with transparency)
#ax.hist(all_not_repaired, bins=100, align='left', color='red', edgecolor='black', range=(1, 101), alpha=0.7, label='Not Repaired')
#ax.hist(all_not_repaired, bins=100, align='left', color='blue', edgecolor='black',range=(1,101))

#Titles
if showlabels:
    axs[0].set_title('AT Base-Pairs', fontsize=title_s)
    axs[1].set_title('GC Base-Pairs', fontsize=title_s)

#AT
ax = axs[0]
ax.hist(AT_repaired, bins=100, align='left', color=repaired_colour, edgecolor='black',range=(1,101), label='Repaired')
ax.hist(AT_not_repaired, bins=100, align='left', color=not_repaired_colour, edgecolor='black', range=(1, 101), alpha=0.7, label='Not Repaired')

#AT
ax = axs[1]
ax.hist(GC_repaired, bins=100, align='left', color=repaired_colour, edgecolor='black',range=(1,101), label='Repaired')
ax.hist(GC_not_repaired, bins=100, align='left', color=not_repaired_colour, edgecolor='black', range=(1, 101), alpha=0.7, label='Not Repaired')

# SORT TICKS
# -----------------------------------------------------------
xticks = np.arange(0, 101, 10, dtype=int)
yticks = np.arange(0, 16, 3, dtype=int)
if showlabels:
    axs[0].legend(loc='best', fontsize=legend_s)
for ax in [axs[0], axs[1]]:#, axs[2]]:
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xlabel('Life time (ns)', fontsize=label_s)
    ax.set_ylabel('Counts', fontsize=label_s)
    ax.grid(True)

#plt.savefig(outfile+".eps")
plt.savefig(outfile+".png",dpi=600)
plt.show()



