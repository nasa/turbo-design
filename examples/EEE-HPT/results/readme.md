# EEE-2 Stage HPT Simulations with ADS

## Background
There are two EEE HPT Designs: A single stage, developed by Pratt and Whitney, and a dual stage, developed by GE. The single stage is not used although the report was very nicely written. At NASA, we are using the dual stage EEE. The files below are all for the dual stage. The purpose of this page is to supply the reader with all the files needed to set up a CFD Simulation so geometry and boundary conditions along with a few results from the experiments to compare with. 

## The code
`main.py` loads up the iges files which were provided in the [EEE report](https://ntrs.nasa.gov/citations/20150003286). This report presents a portion of the Computational Fluid Dynamics (CFD) results for the High Pressure Compressor (HPC). However, CFD data for the High Pressure Turbine (HPT) is currently unavailable. The code processes IGES files by extracting curves and smoothing the geometry, which initially appears jagged with repeated points near the trailing edge. After smoothing, the geometry is converted into an AGF file format, which is compatible with ADS Wand for mesh generation and subsequent simulation using LEO.


## Prerequisites
>   `pip install -r pyturbo-aero numpy scipy matplotlib`


## How to run

`python main.py` 
main.py
- Reads iges files
- Extracts curves
- creates agf
- runs wand
- runs leo
- plots convergence