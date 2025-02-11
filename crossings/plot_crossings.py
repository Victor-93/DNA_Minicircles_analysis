import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
import matplotlib.patches as patches

def add_custom_labels(ax, labels, colors, rect_width=0.05, fontsize=16):
    """
    Adds colored rectangles with labels outside the plot.

    Parameters:
    - ax: The axis to which the labels will be added.
    - labels: List of label names.
    - colors: List of colors corresponding to each label.
    - position: 'right' or 'left' (where to place the labels).
    - spacing: Spacing between labels.
    - rect_width: Width of the rectangle.

    Returns:
    - None (modifies the plot).
    """
    # Get axis limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()


    x_size = abs(xlim[1] - xlim[0])
    y_size = abs(ylim[1] - ylim[0])
    x_start = x_size*.05
    #x_spacing = x_size*0.25
    x_spacing = x_size*0.4
    y_start = -y_size*0.05

    for i, (label, color) in enumerate(zip(labels, colors)):
        y_pos = y_start
        x_rect = x_start + x_spacing * i
        x_text = x_rect + rect_width*x_size*1.2
        y_text = y_pos - y_size*0.02
        # y_pos = y_start - i * y_step * (ylim[1] - ylim[0])

        # Add rectangle
        rect = patches.Rectangle(
            (x_rect, y_pos),  # (x, y) bottom-left corner
            rect_width * (xlim[1] - xlim[0]),  # width
            0.05 * (ylim[1] - ylim[0]),  # height
            linewidth=1,
            edgecolor='black',
            facecolor=color,
            transform=ax.transData,
            clip_on=False  # This ensures it appears outside the plot
        )
        ax.add_patch(rect)

        # Add text next to rectangle
        ax.text(x_text, y_text, label, va='center', ha='left', fontsize=fontsize, transform=ax.transData)

# Description
# ----------------------------------------------------------
# We want to collect crossings from different cases and join them into one graph

# Inputs
# ----------------------------------------------------------
paths = [
    'data/no_nuc_CA-',
    'data/no_nuc_CA_r1-'  #,
    # 'data/two_nuc_CA-',
    # 'data/one_nuc-'
]

title = 'DNA Crossings'

infile = 'color_op1.txt'
outfile = 'crossings'

cmap = ListedColormap(['white', 'red', 'blue'])  #,'green','purple'])#, 'blue'])#, 'green', 'yellow'])

# System info
nbp = 358  # number of base-pairs in minicircle
n_length = 147  # nucleosome length
l_length = int((nbp - n_length * 2) / 2)  # linker length (arms)
nucleosome1_position = 17
nucleosome2_position = 196
nucleosome1_apex = 90
nucleosome2_apex = 269

# Correct apexes
n1_apex = nucleosome1_apex - nucleosome1_position
n2_apex = nucleosome2_apex - nucleosome1_position
apex_color = 'gray'
apex_ls = '--'
apex_lw = 1

# Workout regions, labels, colors and their alphas
tick_positions = [0, n_length - 1, n_length + l_length - 1, n_length * 2 + l_length - 1,
                  nbp - 1]  # Positions on the y-axis
tick_labels = ["NUC 1", "L1", "NUC 2", "L2"]  # Labels for the regions
colors = ["red", "blue", "green", "blue"]  # Colors for the regions

alpha0 = 0.2

# Plotting variables
# -------------------------------------------------------------------
width = 6.5
height = 3.5

tick_s = 12
label_s = 16
title_s = 18
legend_s = 12
subtitle_s = 14
bar_tick = 10
s_shrink = .5
tmin = 1.75
groove_s = 8  # Font for groove
colormap = 'inferno'

# Process
# ----------------------------------------------------------
fig, axs = plt.subplots(1, figsize=(width, height), dpi=300, tight_layout=True)

nuc_color = 'silver'
# axs.plot([0,2103], [17,17], '--', lw=0.5,color=nuc_color,alpha=0.5)
# axs.plot([0,2103], [164,164], '--', lw=0.5,color=nuc_color,alpha=0.5)
# axs.plot([0,2103], [196,196], '--', lw=0.5,color=nuc_color,alpha=0.5)
# axs.plot([0,2103], [343,343], '--', lw=0.5,color=nuc_color,alpha=0.5)

# Initialise
cross_matrix = np.loadtxt(paths[1] + infile)  # Careful, load the biggest one
cross_matrix = np.zeros_like(cross_matrix, dtype=int)

for path in paths:
    icross = np.loadtxt(path + infile)
    a, b = np.shape(icross)
    print(a, b)
    cross_matrix[:, :b] = cross_matrix[:, :b] + icross[:, :b]

fig.suptitle(title, fontsize=title_s)

hbond_matrix = cross_matrix
# Rearrange rows
hbond_matrix = np.vstack((hbond_matrix[nucleosome1_position - 1:], hbond_matrix[:nucleosome1_position - 1]))
cross_matrix = hbond_matrix
im = axs.imshow(cross_matrix, cmap=cmap, aspect=2)

# Draw nucleosome regions
# -----------------------------------------------------------
ax = axs
ax.set_yticks(tick_positions)  # Draw ticks for the nucleosome positions

# Add labels at the middle of each region and draw a rectangle
for k in range(len(tick_labels)):
    start = tick_positions[k]
    end = tick_positions[k + 1]
    label = tick_labels[k]
    midpoint = (start + end) / 2 - 3
    alpha = alpha0
    color = colors[k]
    # Draw region label
    ax.text(-0.5, midpoint, label,  # -0.5 positions the label outside the heatmap
            ha="right", va="center", rotation=90, fontsize=tick_s)

    # Draw region with rectangle
    rect = patches.Rectangle(
        (-60., start),  # x, y: starting position (-1.5 places it outside the heatmap)
        60.,  # width
        end - start,  # height of the rectangle
        linewidth=0,  # No border
        edgecolor=None,
        facecolor=color,
        clip_on=False,
        alpha=alpha  # Transparency
    )
    ax.add_patch(rect)
    # Finally, draw apexes
    time = 2103
    ax.plot([0, time], [n1_apex, n1_apex], lw=apex_lw, ls=apex_ls, color=apex_color, alpha=0.5)
    ax.plot([0, time], [n2_apex, n2_apex], lw=apex_lw, ls=apex_ls, color=apex_color, alpha=0.5)

# Add labels (Rectangles)
# -----------------------------------------------------------

# Call the function to add labels
labels = ["Empty minicircle Ca r1", "Empty minicircle Ca r2"]
colors = ["red", "blue"]
add_custom_labels(ax, labels, colors, fontsize=legend_s)


# SORT TICKS
# -----------------------------------------------------------
xticks = np.arange(0, time + 1, 250, dtype=int)
xtick_labels = np.arange(0, time / 10 + 1, 25, dtype=int)
axs.set_xticks(xticks)
axs.set_xticklabels(xtick_labels, fontsize=tick_s)
axs.set_xlabel('Time (ns)', fontsize=label_s)
axs.set_xlim(0, 2103)
axs.grid(True, alpha=0.5)
axs.tick_params(labelleft=False)

fig.savefig(outfile + '.pdf')
fig.savefig(outfile + '.png')

plt.show()