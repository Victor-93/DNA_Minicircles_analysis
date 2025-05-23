import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import distances

# Let's sort the inputs
info_list = [
    {'parm': "../no_nuc_CA/writhe/no_nuc_nobox.parm7", "coords": "../no_nuc_CA/writhe/no_nuc_100ps.nc", "name": "no_nuc_CA"},
    {'parm': "../no_nuc_CA_r1/writhe/no_nuc_nobox.parm7", "coords": "../no_nuc_CA_r1/writhe/nonuc_CA_r1_100ps.nc", "name": "no_nuc_CA_r1"},
    {'parm': "../no_nuc_Na/writhe/no_nuc_nobox.parm7", "coords": "../no_nuc_Na/writhe/no_nuc_Na_100ps.nc", "name": "no_nuc_Na"},
    {'parm': "../no_nuc_Na/writhe/no_nuc_nobox.parm7", "coords": "../no_nuc_Na_r1/writhe/no_nuc_Na_r1_100ps.nc", "name": "no_nuc_Na_r1"},
    {'parm': "../no_nuc_Na/writhe/no_nuc_nobox.parm7", "coords": "../one_nuc/writhe/one_nuc_100ps.nc", "name": "one_nuc_Na"},
    {'parm': "../no_nuc_Na/writhe/no_nuc_nobox.parm7", "coords": "../two_nuc_CA/writhe/two_nuc_100ps.nc", "name": "two_nuc_CA"},
    {'parm': "../no_nuc_Na/writhe/no_nuc_nobox.parm7", "coords": "../two_nuc_Na/writhe/two_nuc_Na_100ps.nc", "name": "two_nuc_Na"}
]

# Process
for info_dict in info_list:
    parmfile = info_dict['parm']
    coordsfile = info_dict['coords']
    name = info_dict['name']

    print("Doing " + name)


    coords = mda.Universe(parmfile, coordsfile)

    dist = []
    time = []
    base = []
    atom = []

    n_frames = coords.trajectory.n_frames
    f = open("../dna_2hbonds_2dist.txt", "r")
    for i, line in enumerate(f):
        res1 = int(line.split()[0])
        atom1 = line.split()[2]
        res2 = line.split()[3]
        atom2 = line.split()[5]

        for ts in coords.trajectory[1:n_frames:10]:
        #for ts in coords.trajectory[1:1500:10]:
            select1 = coords.select_atoms('resid {} and name {}'.format(res1, atom1))
            select2 = coords.select_atoms('resid {} and name {}'.format(res2, atom2))
            resids1, resids2, dist1 = distances.dist(select1, select2)
            base_n = str(res1)
            base_z = base_n.zfill(3)
            time.append(int((ts.frame - 1) / 10))
            atom.append(str(atom1))
            base.append(str(base_z))
            dist.append(float(dist1))

    dist_df = pd.DataFrame(list(zip(base, atom, time, dist)), columns=['base', 'atom', 'time', 'distance'])

    pivot_dist_df = dist_df.pivot(index=['base', 'atom'], columns='time', values='distance')

    hbonds_df = pivot_dist_df.groupby('base').mean()
    hbonds_df.index = np.arange(0, len(hbonds_df))
    hbonds_df.to_csv("hbonds_{}.csv".format(name), index=False)

