import MDAnalysis as mda

# Description
# ----------------------------------------------------------------------------------------------------------------------
# We want to make a script that creates a pdb of the last frame of the trajectory

# Let's sort the inputs
info_list = [
    # Sarah's
    {'parm': "../../../no_nuc_CA/writhe/no_nuc_nobox.parm7", "coords": "../../../no_nuc_CA/writhe/no_nuc_100ps.nc", "name": "no_nuc_CA"},
]

filter_frame = 10  # Which frame to filter

for info_dict in info_list:
    topology = info_dict['parm']
    coordinates = info_dict['coords']
    name = info_dict['name']
    output_pdb = name+"_frame-"+str(filter_frame)+"_command.pdb"


    # PDB part
    # ------------------------------------------------------------------------------------------------------------------
    # Load the trajectory
    u = mda.Universe(topology, coordinates)

    # Save the trajectory as a .mdcrd file
    #with mda.Writer("output.crd", u.atoms.n_atoms) as writer:
    #    for ts in u.trajectory:
    #        writer.write(u.atoms)

    print("Trajectory saved as 'output.mdcrd'")

    # Write to Amber .prmtop and .inpcrd formats
    u.save("output.prmtop", format="amber")
