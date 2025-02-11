import numpy as np
import pandas as pd
import seaborn as sns
from scipy.ndimage import convolve
import matplotlib.pyplot as plt

#Description
# ----------------------------------------------------------------------------------
#This version adjusts vmin and vmax so it looks more binary

min_val = 4  # Any value below this will be 0 and above 1
#enhance_res = True
enhance_res = False

if enhance_res:
    # Convolution kernel to expand 1s to neighboring rows
    #kernel = np.array([[0, 1, 0], [0, 1, 0], [0, 1, 0]])
    kernel = np.array([[1], [1], [1]])

#Figure parms
# ----------------------------------------------------------------------------------
width = 8
height = 4
cmap = 'plasma'
nbp = 358
tick_s = 10
label_s = 14
title_s = 16

#yticks = np.linspace(0, nbp, 10, dtype=int)
yticks = np.arange(0, nbp+1, 50, dtype=int)

# And plot
fig, axs = plt.subplots(3,2, figsize=(width*2, height*3), tight_layout=True)

# List of dicts with info
info_list = [
    {'name':"no_nuc_CA", "title":"No Nucleosome Ca", "ax": axs[0,0]},
    {'name': "no_nuc_CA_r1", "title": "No Nucleosome Ca (r1)", "ax": axs[0, 1]},
    {'name': "no_nuc_Na", "title": "No Nucleosome Na", "ax": axs[1, 0]},
    {'name': "no_nuc_Na_r1", "title": "No Nucleosome Na (r1)", "ax": axs[1, 1]},
    {'name': "two_nuc_CA", "title": "Two Nucleosomes Ca", "ax": axs[2, 0]},
    {'name': "two_nuc_Na", "title": "Two Nucleosomes Na", "ax": axs[2, 1]},
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

    time = len(hbond_df.iloc[0])
    xticks = np.linspace(0, time, 10, dtype=int)

    sns.heatmap(hbond_df, ax=ax, cmap=cmap)#, vmin=3, vmax=4)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks, fontsize=tick_s)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks, fontsize=tick_s)
    ax.set_xlabel('Time (ns)', fontsize=label_s)
    ax.set_ylabel('DNA Base Pair', fontsize=label_s)
    ax.grid(True, alpha=0.2)


if enhance_res:
    plt.savefig("hbonds_binary_e.pdf")
    plt.savefig("hbonds_binary_e.png")
else:
    plt.savefig("hbonds_binary.pdf")
    plt.savefig("hbonds_binary.png")
plt.show()



