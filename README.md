# Optimization of Embedded Element Patterns of Reactively Loaded Antenna Arrays

This repositery contains MATLAB codes and data for reproducing the results of the research paper 

> A. Salmi, M. Capek, L. Jelinek, A. Lehtovuori, and V. Viikari, "Optimization of Embedded Element Patterns of Reactively Loaded Antenna Arrays," in *Transactions on Antennas and Propagation*, 2024 [Submitted].

The codes optimize embedded element patterns of reactively loaded antenna arrays, and compare multiple optimization methods.

## Requirements
The codes require the following MATLAB toolboxes and packages to be installed in the system and added to the MATLAB path:
- [RF Toolbox](https://se.mathworks.com/products/rftoolbox.html)
- [Manopt 7.1](https://www.manopt.org/downloads.html)
- [CVX](https://cvxr.com/cvx/download/)
- [matlab2tikz](https://se.mathworks.com/matlabcentral/fileexchange/22022-matlab2tikz-matlab2tikz) (Only required for exporting the figures in Tikz format from MATLAB.)

## Usage
The reactive terminations of the passive antenna elements are optimized in script *compute_terminations.m*. In the beginning of the script, you must specify the folders where to save the results, and where the simulation data can be found from. The results are saved as a MATLAB workspace.

The figures are plotted in the scripts in the folder *plot_scripts*. Each plot script requires path to the corresponding matlab workspace file.