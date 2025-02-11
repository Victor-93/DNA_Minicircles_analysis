import numpy as np
import pandas as pd
import seaborn as sns
from scipy.ndimage import convolve
import matplotlib.pyplot as plt

#Description
# ----------------------------------------------------------------------------------
#This version adjusts vmin and vmax so it looks more binary


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
fig, axs = plt.subplots(3, figsize=(width, height*3), tight_layout=True)

# List of dicts with info
info_list = [
    {'name':"Tom_two-nucs", "title":"Tom's simulation: two-nucs", "ax": axs[0]},
    {'name': "Tom_two-nucs-minus", "title": "Tom's simulation: two-nucs-minus", "ax": axs[1]},
    {'name': "Tom_two-nucs-plus", "title": "Tom's simulation: two-nucs-plus", "ax": axs[2]}
]

for i, info_dict in enumerate(info_list):
    ax = info_dict['ax']
    title = info_dict['title']
    name = info_dict['name']

    ax.set_title(title,fontsize=title_s)

    # Let's load the hbonds dfs
    hbond_df = pd.read_csv("hbonds_{}.csv".format(name))

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


plt.savefig("Tom-hbonds.pdf")
plt.savefig("Tom-hbonds.png")
plt.show()



