import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Let's sort the inputs
info_list = [
    {'path': "../../Tom_simulations/dinucsims/two-nucs/Tw-Wr/", "name": "two-nucs", "color": "black"},
    {'path': "../../Tom_simulations/dinucsims/two-nucs-minus/Tw-Wr/", "name": "two-nucs-minus", "color": "red"},
    {'path': "../../Tom_simulations/dinucsims/two-nucs-plus/Tw-Wr/", "name": "two-nucs-plus", "color": "green"}
]

B_DNA_turn = 10.5 # How many base-pairs per B-DNA turn
#Figure parms
# ----------------------------------------------------------------------------------
width = 10
height = 4
cmap = 'plasma'
nbp = 358
tick_s = 10
label_s = 14
title_s = 16

#yticks = np.linspace(0, nbp, 10, dtype=int)
yticks = np.arange(0, nbp+1, 50, dtype=int)

# And plot
fig, axs = plt.subplots(5, figsize=(width*1, height*5), tight_layout=True)

# Let's try to plot while we collect
# ------------------------------------------------------------------------------------------------
ax = axs[0]
ax.set_title("Twist", fontsize=title_s)
for i, info_dict in enumerate(info_list):
    color = info_dict['color']
    label = info_dict['name']
    twist_file = info_dict['path'] + 'tw.ser'
    writhe_file = info_dict['path'] + 'writhe.ser'
    twist_data = np.loadtxt(twist_file)
    writhe_data = np.loadtxt(writhe_file)
    nbp=twist_data.shape[1]  # Number of base-pairs
    n = twist_data.shape[0]

    time  = np.arange(n)/10

    # Calculate total twist, extract writhe and calculate linking number
    total_twist = np.sum(twist_data, axis=1)/360
    writhe = writhe_data[:,1]
    Lk = total_twist + writhe
    Lk0 = nbp/B_DNA_turn
    dLk = Lk - Lk0
    superhelical = dLk/Lk0

    # Let's plot
    axs[0].plot(time,total_twist, color=color, label=label)
    axs[1].plot(time,writhe, color=color, label=label)
    axs[2].plot(time,Lk, color=color, label=label)
    axs[3].plot(time,dLk, color=color, label=label)
    axs[4].plot(time,superhelical, color=color, label=label)

for i in range(5):
    axs[i].grid(True, alpha=0.2)
    axs[i].set_xlabel('Time (ns)', fontsize=label_s)

axs[0].set_ylabel('Twist', fontsize=label_s)
axs[1].set_ylabel('Writhe', fontsize=label_s)
axs[2].set_ylabel(r'$Lk$', fontsize=label_s)
axs[3].set_ylabel(r'$\Delta Lk$', fontsize=label_s)
axs[4].set_ylabel('Superhelical density', fontsize=label_s)
axs[0].legend(loc='best')
plt.savefig("Tom_tw_wr_Lk.pdf")
plt.savefig("Tom_tw_wr_Lk.png")

#axs[4].set_xlabel('Time (ns)', fontsize=label_s)
plt.show()