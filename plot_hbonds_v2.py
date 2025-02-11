import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt

# Inputs
# ----------------------------------------------------------------------------------
Na_list = ["no_nuc_Na", "no_nuc_Na_r1", "two_nuc_Na"]
Ca_list = ["no_nuc_CA", "no_nuc_CA_r1", "two_nuc_CA"]

#Figure parms
# ----------------------------------------------------------------------------------
width = 8
height = 4
cmap = 'plasma'
nbp = 358
tick_s = 10
label_s = 14
title_s = 16

yticks = np.linspace(0, nbp, 6, dtype=int)

# And plot
fig, axs = plt.subplots(3,2, figsize=(width*2, height*3), tight_layout=True)

# First Na
titles = ["No Nucleosome Na", "No Nucleosome Na (r1)", "Two Nucleosomes Na"]
for i, name in enumerate(Na_list):
    ax = axs[i,0]
    title = titles[i]
    ax.set_title(title,fontsize=title_s)

    # Let's load the hbonds dfs
    hbond_df = pd.read_csv("hbonds_{}.csv".format(name))

    time = len(hbond_df.iloc[0])
    xticks = np.linspace(0, time, 10, dtype=int)

    sns.heatmap(hbond_df, ax=ax, cmap=cmap, vmin=3, vmax=10)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks, fontsize=tick_s)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks, fontsize=tick_s)
    ax.set_xlabel('Time (ns)', fontsize=label_s)
    ax.set_ylabel('DNA Base Pair', fontsize=label_s)

# Second Ca
titles = ["No Nucleosome Ca", "No Nucleosome Ca (r1)", "Two Nucleosomes Ca"]
for i, name in enumerate(Ca_list):
    ax = axs[i,1]
    title = titles[i]
    ax.set_title(title,fontsize=title_s)

    # Let's load the hbonds dfs
    hbond_df = pd.read_csv("hbonds_{}.csv".format(name))

    time = len(hbond_df.iloc[0])
    xticks = np.linspace(0, time, 10, dtype=int)

    sns.heatmap(hbond_df, ax=ax, cmap=cmap, vmin=3, vmax=10)

    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks, fontsize=tick_s)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks, fontsize=tick_s)
    ax.set_xlabel('Time (ns)', fontsize=label_s)
    ax.set_ylabel('DNA Base Pair', fontsize=label_s)

#plt.savefig("hbonds.pdf")
#plt.savefig("hbonds.png")
plt.show()



