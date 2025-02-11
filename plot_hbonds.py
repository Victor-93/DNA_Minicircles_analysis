import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt

#Inputs
names_list = [
   # "no_nuc_CA",
   # "no_nuc_CA_r1",
   # "no_nuc_Na",
  #  "no_nuc_Na_r1",
  #  ####"one_nuc_Na",
    "two_nuc_CA",
    "two_nuc_Na"
]
#names_list = [
#    "Tom_two-nucs",
#    "Tom_two-nucs-minus",
#    "Tom_two-nucs-plus"
#]


#Figure parms
width = 8
height = 4
cmap = 'cool'
nbp = 358
yticks = np.linspace(0, nbp, 6, dtype=int)
print(yticks)

for name in names_list:
    # Let's load the hbonds dfs
    hbond_df = pd.read_csv("hbonds_{}.csv".format(name))

    time = len(hbond_df.iloc[0])
    xticks = np.linspace(0, time, 10, dtype=int)

    # And plot
    fig, ax = plt.subplots(figsize=(width, height), tight_layout=True)
    sns.heatmap(hbond_df, ax=ax, xticklabels=10, yticklabels=10, cmap=cmap)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('DNA Base Pair')
    plt.show()



