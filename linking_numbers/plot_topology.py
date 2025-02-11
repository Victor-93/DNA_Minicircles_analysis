import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Let's sort the inputs
info_list = [
    {'path': "../../no_nuc_CA/writhe/no_nuc/", "name": "no_nuc_CA", "color": "blue"},
    {'path': "../../no_nuc_CA_r1/writhe/no_nuc_CA_r1/", "name": "no_nuc_CA_r1", "color": "red"},
    {'path': "../../no_nuc_Na/writhe/no_nuc_NA/", "name": "no_nuc_Na", "color": "green"},
    {'path': "../../no_nuc_Na_r1/writhe/no_nuc_NA_r1/", "name": "no_nuc_Na_r1", "color": "yellow"},
    {'path': "../../one_nuc/writhe/one_nuc/", "name": "one_nuc_Na", "color": "purple"},
    {'path': "../../two_nuc_CA/writhe/two_nuc/", "name": "two_nuc_CA", "color": "black"},
    {'path': "../../two_nuc_Na/writhe/two_nuc_Na/", "name": "two_nuc_Na", "color": "orange"}
]

B_DNA_turn = 10.5 # How many base-pairs per B-DNA turn
#Figure parms
# ----------------------------------------------------------------------------------
width = 6
height = 3
cmap = 'plasma'
nbp = 358
tick_s = 10
label_s = 14
title_s = 16

#yticks = np.linspace(0, nbp, 10, dtype=int)
yticks = np.arange(0, nbp+1, 50, dtype=int)

# And plot
fig, axs = plt.subplots(2, 2, figsize=(width*2, height*2), tight_layout=True, sharex=True)

# Let's try to plot while we collect
# ------------------------------------------------------------------------------------------------
fig.suptitle("DNA Topology", fontsize=title_s)
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
    axs[0,0].plot(time,total_twist, color=color, label=label)
    axs[0,1].plot(time,writhe, color=color, label=label)
#    axs[2].plot(time,Lk, color=color, label=label)
    axs[1,0].plot(time,dLk, color=color, label=label)
    axs[1,1].plot(time,superhelical, color=color, label=label)

axs[0,0].grid(True, alpha=0.2)
axs[0,1].grid(True, alpha=0.2)
axs[1,0].grid(True, alpha=0.2)
axs[1,1].grid(True, alpha=0.2)
axs[1,0].set_xlabel('Time (ns)', fontsize=label_s)
axs[1,1].set_xlabel('Time (ns)', fontsize=label_s)

axs[0,0].set_ylabel('Twist', fontsize=label_s)
axs[0,1].set_ylabel('Writhe', fontsize=label_s)
#axs[2].set_ylabel(r'$Lk$', fontsize=label_s)
axs[1,0].set_ylabel(r'$\Delta Lk$', fontsize=label_s)
axs[1,1].set_ylabel('Superhelical density', fontsize=label_s)
axs[0,0].legend(loc='best')
plt.savefig("topology_calcs.pdf")
plt.savefig("topology_calcs.png")

#axs[4].set_xlabel('Time (ns)', fontsize=label_s)
plt.show()