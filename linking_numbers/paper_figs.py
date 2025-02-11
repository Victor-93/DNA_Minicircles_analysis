import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Let's sort the inputs
info_list = [
    {'path': "../../no_nuc_CA/writhe/no_nuc/", "name": "Empty Ca r1", "color": "blue"},
    {'path': "../../no_nuc_CA_r1/writhe/no_nuc_CA_r1/", "name": "Empty Ca r2", "color": "red"},
    {'path': "../../no_nuc_Na/writhe/no_nuc_NA/", "name": "Empty Na r1", "color": "green"},
    {'path': "../../no_nuc_Na_r1/writhe/no_nuc_NA_r1/", "name": "Empty Na r2", "color": "yellow"},
    {'path': "../../one_nuc/writhe/one_nuc/", "name": "Mononuc Na", "color": "purple"},
    {'path': "../../two_nuc_CA/writhe/two_nuc/", "name": "Dinuc Ca", "color": "black"},
    {'path': "../../two_nuc_Na/writhe/two_nuc_Na/", "name": "Dinuc Na", "color": "orange"}
]

B_DNA_turn = 10.5 # How many base-pairs per B-DNA turn
t_name = "Dinuc CA"  # Name that we want to skip first measurement
#Figure parms
# ----------------------------------------------------------------------------------
width = 6.5
height = 3.5
cmap = 'plasma'
nbp = 358
tick_s = 12
label_s = 16
title_s = 18
legend_s=9

#yticks = np.linspace(0, nbp, 10, dtype=int)
yticks = np.arange(0, nbp+1, 50, dtype=int)

# Average dLK and average superhelical
# ***************************************
dLK_array = []  # These will store the average linking difference and superhelical density of all simulations
                # within the first 100ns
superhelical_array = []

# Figure 1: Twist and Writhe
# ------------------------------------------------------------------------------------------------
# And plot
fig, axs = plt.subplots(1, 2, figsize=(width*2, height*1), tight_layout=True, sharex=True)

outfile = 'Tw-Wr'
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

    if label==t_name:
        total_twist = total_twist[1:]
        superhelical = superhelical[1:]
        writhe = writhe[1:]
        dLk = dLk[1:]
        time =time[1:]

    # Let's plot
    axs[0].plot(time,total_twist, color=color, label=label)
    axs[1].plot(time,writhe, color=color, label=label)

    # Collect measurements for overalls
    dLK_array.append(dLk[:1000])
    superhelical_array.append(superhelical[:1000])

axs[0].grid(True, alpha=0.2)
axs[1].grid(True, alpha=0.2)
axs[0].set_xlabel('Time (ns)', fontsize=label_s)
axs[1].set_xlabel('Time (ns)', fontsize=label_s)

axs[0].set_ylabel('Twist', fontsize=label_s)
axs[1].set_ylabel('Writhe', fontsize=label_s)
axs[0].legend(loc='best',fontsize=legend_s)
plt.savefig(outfile+".pdf")
plt.savefig(outfile+".png")
plt.show()


# ********* ARRAY PART
dLK_array = np.concatenate(dLK_array)
superhelical_array = np.concatenate(superhelical_array)

print("Average Linking Difference", np.mean(dLK_array), np.std(dLK_array))
print("Average Superhelical", np.mean(superhelical_array), np.std(superhelical_array))

# Supfigure: dLK and Superhelical
# ------------------------------------------------------------------------------------------------
# And plot
fig, axs = plt.subplots(1, 2, figsize=(width*2, height*1), tight_layout=True, sharex=True)

outfile = 'dLK-sigma'
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

    if label==t_name:
        total_twist = total_twist[1:]
        superhelical = superhelical[1:]
        writhe = writhe[1:]
        dLk = dLk[1:]
        time =time[1:]

    # Let's plot
    axs[0].plot(time,dLk, color=color, label=label)
    axs[1].plot(time,superhelical, color=color, label=label)

axs[0].grid(True, alpha=0.2)
axs[1].grid(True, alpha=0.2)
axs[0].set_xlabel('Time (ns)', fontsize=label_s)
axs[1].set_xlabel('Time (ns)', fontsize=label_s)

axs[0].set_ylabel(r'$\Delta Lk$', fontsize=label_s)
axs[1].set_ylabel('Superhelical density', fontsize=label_s)
axs[0].legend(loc='best',fontsize=legend_s)
plt.savefig(outfile+".pdf")
plt.savefig(outfile+".png")
plt.show()
