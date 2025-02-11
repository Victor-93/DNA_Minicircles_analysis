import MDAnalysis as mda

# Description
# ----------------------------------------------------------------------------------------------------------------------
# We want to make a script that creates a pdb of the last frame of the trajectory

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
    output_pdb = 'last-frame_'+ name + '.pdb'

    # PDB part
    # ------------------------------------------------------------------------------------------------------------------
    # Load the trajectory
    u = mda.Universe(topology, coordinates)

    # Select the last frame
    u.trajectory[-1]  # Move to the last frame

    # Save the last frame as a PDB file
    with mda.Writer(output_pdb, multiframe=False) as pdb_writer:
        pdb_writer.write(u.atoms)

    print(f"Last frame saved as: {output_pdb}")