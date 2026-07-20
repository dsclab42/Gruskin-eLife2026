# Gruskin-eLife2026

Analysis code for:

> Gruskin DC, Vieira DJ, Lee JK, Patel GH (2026).
> **Heritability of movie-evoked brain activity and connectivity.**
> *eLife* 14. https://doi.org/10.7554/eLife.106081

This repository contains the MATLAB, R, and Python code used to quantify the
heritability of movie-evoked brain activity (intersubject correlation of BOLD
time courses) and functional connectivity, and to decompose that heritability
into genetic control over cortical topography and neural timescale, using 7T
fMRI data from the Human Connectome Project (HCP) Young Adult twin sample.

## Repository structure

- `matlab_code/` — ISC/FC computation, parcellation, hyperalignment, and neural timescale analyses
- `r_code/` — multidimensional heritability (h²) estimation
- `python_code/` — downstream analysis and figure generation (Jupyter notebooks)

## Software

Analyses were performed in MATLAB (R2023b), Python, and R. Cortical surface
visualizations were produced with Connectome Workbench.

## Data availability

The raw neuroimaging data come from the Human Connectome Project (HCP) Young
Adult 7T release and can be downloaded from ConnectomeDB
(https://db.humanconnectome.org) under the HCP data use terms. This repository
contains analysis code only and does not redistribute HCP data.

## Archived version

A permanent snapshot of this code is archived at Software Heritage.

## Citation

If you use this code, please cite the paper above.

## License

This project is released under the MIT License. See the [`LICENSE`](LICENSE) file for details.
