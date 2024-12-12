# Optimization of Embedded Element Patterns of Reactively Loaded Antenna Arrays

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14412066.svg)](https://doi.org/10.5281/zenodo.14412066)

This repository contains MATLAB codes and data for reproducing the results of the research paper 

> A. Salmi, M. Capek, L. Jelinek, A. Lehtovuori, and V. Viikari, "Optimization of Embedded Element Patterns of Reactively Loaded Antenna Arrays," in *Transactions on Antennas and Propagation*, 2024 [Submitted].

The codes optimize embedded element patterns of reactively loaded antenna arrays and compare multiple optimization methods.

## Requirements
The codes require the following MATLAB toolboxes and packages to be installed in the system and added to the MATLAB path:
- [RF Toolbox](https://se.mathworks.com/products/rftoolbox.html)
- [Manopt 7.1](https://www.manopt.org/downloads.html)
- [CVX](https://cvxr.com/cvx/download/)
- [matlab2tikz](https://se.mathworks.com/matlabcentral/fileexchange/22022-matlab2tikz-matlab2tikz) (Only required for exporting the figures in Tikz format from MATLAB.)

## Usage
The reactive terminations of the passive antenna elements are optimized in script *compute_terminations.m*. At the beginning of the script, you must specify the folders where to save the results, and where the simulation data can be found from. In addition, you should define the target sector and gain. 

The optimization results are saved as a MATLAB workspace after running *compute_terminations.m*. The figures of the paper are plotted using the scripts in the folder *plot_scripts*. Each plot script requires a path to the corresponding Matlab workspace file. The plotting scripts generate tsv data files and a tex file, which are processed using Latex.

There are two additional scripts: *compare_time.m* and *tune_manopt.m*. The first one is used to generate the computational complexity analysis of the studied optimization methods. The latter one is for improving the MO algorithm's convergence speed by tuning the hyperparameters.
