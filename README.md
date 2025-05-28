# DNA Minicircles Analysis

This repository contains Python scripts for analysing molecular dynamics (MD) simulations of DNA minicircles, used in the manuscript *"Nucleosome-mediated conformational switches in ultra-small eccDNAs."*

The toolkit includes:

* **Hydrogen bond analysis** in terms of inter-base-pair distances.
* **Topology calculations**, including twist, writhe, and superhelical density.
* **Crossing analyses**, identifying regions where two DNA segments come into close contact (writhed structures).

If you find the scripts, data, or methodology useful, please cite our work.
Feel free to use and modify the scripts provided here.

---

## Requirements

The main Python libraries used are:

* `MDAnalysis`
* `Pandas`
* `NumPy`
* `Matplotlib`

---

## Scripts and Methodology

The analysis scripts can be run independently, although some rely on preprocessed data or outputs from others.

---

### Hydrogen Bond Calculations

Hydrogen bonds are assessed by calculating average distances between heavy atoms:

* **AT base-pairs**: N–O and N–N distances.
* **GC base-pairs**: O–N and N–O distances (excluding the middle hydrogen bond).

If the average distance exceeds **4 Å**, the base-pair is classified as a **DNA defect**.

**Scripts:**

* `calculate_hbonds.py`: Computes inter-base-pair distances every 10 frames. Outputs CSV files with averaged distances.
* `calculate_hbonds_high_res.py`: Similar to the above but computes distances for **every frame** (higher temporal resolution). Intended for analysing hydrogen bond lifetimes.
* `plot_hbond_matrix.py`: Converts inter-base-pair distances to binary matrices (1 = defect, 0 = intact) and plots heatmaps. The threshold is defined by `min_val=4`.

---

### Hydrogen Bond Lifetime

Quantifies how long defects persist across simulations.

**Script:**

* `hbond_lifetime/plot_hbond_distribution_all_bysequence_and_repair_length.py`: Takes high-res CSV files from `calculate_hbonds_high_res.py`, classifies defects by base-pair type (AT or GC) and by repair status:

  * **Repaired**: Defects that resolve during the simulation.
  * **Not\_repaired**: Defects that persist.

---

### Topology Calculations

Includes calculation of:

* **Twist**
* **Writhe**
* **Superhelical density**
* **Linking difference**

These calculations rely on preprocessed outputs from the [WrLINE](https://github.com/agnesnoy/reWrLINE) software, which must be run separately to obtain twist and writhe.

**Note:** Total twist is calculated as the sum of base-step twist values.

**Script:**

* `linking_numbers/paper_figs_topology.py`: Processes WrLINE outputs and computes/plots the topological quantities.

---

### DNA Crossings

Defines DNA crossings based on proximity:

* If two base-pairs are **≤ 10 Å apart** and **≥ 90 bases apart** along the contour (to avoid short-range contacts), the region is flagged as a **crossing**.

The contour is based on base-pair positions derived from WrLINE outputs.

**Script:**

* `crossings/plot_crossings.py`: Detects and visualises crossing regions from MD trajectories.

