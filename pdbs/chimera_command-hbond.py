import pandas as pd
import numpy as np

#Description
# ----------------------------------------------------------------------------------
#This version adjusts vmin and vmax so it looks more binary

min_val = 4  # Any value below this will be 0 and above 1
# Description
# ----------------------------------------------------------------------------------------------------------------------
# This script will create a chimera script that will load the pdb structure (last frame) and will color the base-pairs
# that would denature

# Let's sort the inputs
info_list = [
    # Sarah's
    {'parm': "../../no_nuc_CA/writhe/no_nuc_nobox.parm7", "coords": "../../no_nuc_CA/writhe/no_nuc_100ps.nc", "name": "no_nuc_CA"},
    {'parm': "../../no_nuc_CA_r1/writhe/no_nuc_nobox.parm7", "coords": "../../no_nuc_CA_r1/writhe/nonuc_CA_r1_100ps.nc", "name": "no_nuc_CA_r1"},
    {'parm': "../../no_nuc_Na/writhe/no_nuc_nobox.parm7", "coords": "../../no_nuc_Na/writhe/no_nuc_Na_100ps.nc", "name": "no_nuc_Na"},
    {'parm': "../../no_nuc_Na/writhe/no_nuc_nobox.parm7", "coords": "../../no_nuc_Na_r1/writhe/no_nuc_Na_r1_100ps.nc", "name": "no_nuc_Na_r1"},
    {'parm': "../../no_nuc_Na/writhe/no_nuc_nobox.parm7", "coords": "../../one_nuc/writhe/one_nuc_100ps.nc", "name": "one_nuc_Na"},
    {'parm': "../../no_nuc_Na/writhe/no_nuc_nobox.parm7", "coords": "../../two_nuc_CA/writhe/two_nuc_100ps.nc", "name": "two_nuc_CA"},
    {'parm': "../../no_nuc_Na/writhe/no_nuc_nobox.parm7", "coords": "../../two_nuc_Na/writhe/two_nuc_Na_100ps.nc", "name": "two_nuc_Na"},

    # Toms
    {'parm': "../../Tom_simulations/dinucsims/two-nucs/sys.parm7",
     "coords": "../../Tom_simulations/dinucsims/two-nucs/summary-0.1ns-per-frame.nc",
     "name": "Tom_two-nucs"},
    {'parm': "../../Tom_simulations/dinucsims/two-nucs-minus/sys.parm7",
     "coords": "../../Tom_simulations/dinucsims/two-nucs-minus/summary-0.1ns-per-frame.nc",
     "name": "Tom_two-nucs-minus"},
    {'parm': "../../Tom_simulations/dinucsims/two-nucs-plus/sys.parm7",
     "coords": "../../Tom_simulations/dinucsims/two-nucs-plus/summary-0.1ns-per-frame.nc",
     "name": "Tom_two-nucs-plus"}
]

for info_dict in info_list:
    topology = info_dict['parm']
    coordinates = info_dict['coords']
    name = info_dict['name']
    script_filename = name+"_command.csc"

    # Process
    # ------------------------------------------------------------------------------------------------------------------
    # Let's load the hbonds dfs
    hbond_df = pd.read_csv("../hbonds_{}.csv".format(name))
    hbond_df = (hbond_df >= min_val).astype(int) #Filter

    hbonds = hbond_df.to_numpy()
    x = hbonds[:,-1]
    nbp = len(x)

    # Find indices where x == 1
    indices = np.where(x == 1)[0]

    # Create the lines for the script
    script_lines = [
        "background solid white",
        "open last-frame_" + name + '.pdb',
        "display",
        "nucleotide sidechain fill/fill",
        "color blue #0"
    ]

    # Add the line with the indices
    index_line = f"color yellow :{','.join(map(str, indices))}"
    script_lines.append(index_line)

    # Get complementary bases:
    if len(indices) >0:
        #cindices = 2*nbp -1  - indices
        cindices = 2*nbp +1 - indices
        index_line = f"color yellow :{','.join(map(str, cindices))}"
        script_lines.append(index_line)

    # Save the script to a file
    with open(script_filename, "w") as script_file:
        script_file.write("\n".join(script_lines))


