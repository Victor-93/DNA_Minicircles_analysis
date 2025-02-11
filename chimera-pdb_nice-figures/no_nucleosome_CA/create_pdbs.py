import MDAnalysis as mda

# Description
# ----------------------------------------------------------------------------------------------------------------------
# We want to make a script that creates a pdb of the last frame of the trajectory

# Let's sort the inputs
info_list = [
    # Sarah's
    {'parm': "../../../no_nuc_CA/writhe/no_nuc_nobox.parm7", "coords": "../../../no_nuc_CA/writhe/no_nuc_100ps.nc", "name": "no_nuc_CA"},
]

filter_frame = 130  # Which frame to filter

for info_dict in info_list:
    topology = info_dict['parm']
    coordinates = info_dict['coords']
    name = info_dict['name']
    output_pdb = name+"_frame-"+str(filter_frame)+"_command.pdb"


    # PDB part
    # ------------------------------------------------------------------------------------------------------------------
    # Load the trajectory
    u = mda.Universe(topology, coordinates)

    # Select the last frame
    u.trajectory[filter_frame-1]  # Move to the last frame

    # Save the last frame as a PDB file
    with mda.Writer(output_pdb, multiframe=False) as pdb_writer:
        pdb_writer.write(u.atoms)

    print(f"Last frame saved as: {output_pdb}")