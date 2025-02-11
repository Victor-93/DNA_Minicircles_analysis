import pandas as pd
import numpy as np

# ----------------------------------------------------------------------------------------------------------------------
# This script will create a chimera script that will load the pdb structure and will color the base-pairs
# that would denature.
# The purpose of this script is to be modified, so you can adjust the scripts to produce pdbs that tell the
# evolution of hbonds with pictures.
# So, just tell it which frame you want to filter and it'll highlight the base-pairs that are denaturated in that frame

min_val = 4  # Any value below this will be 0 and above 1

filter_frame = 80 # Which frame to filter

# Description

# Let's sort the inputs
info_list = [
    # Sarah's
    {'parm': "../../../no_nuc_Na/writhe/no_nuc_nobox.parm7", "coords": "../../../no_nuc_Na/writhe/no_nuc_Na_100ps.nc", "name": "no_nuc_Na"}
]

for info_dict in info_list:
    topology = info_dict['parm']
    coordinates = info_dict['coords']
    name = info_dict['name']
    script_filename = name+"_frame-"+str(filter_frame)+"_command.csc"
    pdb_name = name+"_frame-"+str(filter_frame)+"_command.pdb"

    # Process
    # ------------------------------------------------------------------------------------------------------------------
    # Let's load the hbonds dfs
    hbond_df = pd.read_csv("../../hbonds_{}.csv".format(name))
    hbond_df = (hbond_df >= min_val).astype(int) #Filter

    hbonds = hbond_df.to_numpy()

    x = hbonds[:,filter_frame-1]
    nbp = len(x)

    # Find indices where x == 1
    indices = np.where(x == 1)[0]

    # Create the lines for the script
    script_lines = [
        "background solid white",
        "open " + pdb_name,
        "display",
        "nucleotide sidechain fill/fill",
        "color blue #0"
    ]

    # Add the line with the indices
    if len(indices) > 0:
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


