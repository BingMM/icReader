# IMAGE conductance (Builder)

icBuilder is a code to process WIC, SI12, and SI13 data from IMAGE and estimated ionospheric condutance with propagated uncertainties.

## Project description

This code was created to robustly estimate ionospheric conductances. All data between ? and ? have been processed and is available at. The intension with the code is not for others to run it (however, they can), but simply documentation of the data processing procedure.

## Dependencies

The main dependency is fuvpy (https://github.com/aohma/fuvpy).
The only other dependency is tqdm (https://github.com/tqdm/tqdm), used for creating progres bars, and can be omitted with minor changes to the code.

## Install

mamba activate your_environment
git clone icBuilder
cd icBulder
pip install -e .

## Step-by-step guide

The code was written to be run in the following sequence.
### make_orbit_h5_files.py
Looks through available WIC, SI12, SI13 data and assing each data file to and orbit based on orbitsdates.csv
### make_orbit_nc_files.py
Reads the WIC, SI12, and SI13 data, applies a background removal algorithm, saves the results into a series of netcdf files.
### make_background_removal_figures.py (optional)
Plots the data in the netcdf files made above
### grid_resolution_determination.py (optional)
Analyzes a series of the netcdf files, generated above, to determine optimal grid resolution for binned statistics of WIC, SI12, and SI13 data based on the trade-off between standard error and the amount of bins with 30 or more measurements (lower limit of the central limit theorem). Suggested grid resolution is 225 km for WIC and 450 for SI12 and SI13.
### make_conductance_orbit_files.py
Ingests the WIC, SI12, and SI13 netcdf files generates above and generates estimates of conductance
### make_spline_model.py
Not available yet

## Literature

Ohma, A., Madelaire, M., Laundal, K. M., Reistad, J. P., Hatch, S. M., Gasparini, S., and Walker, S. J. (2024). Background removal from global auroral images: Data-driven dayglow modeling. Earth Planet. Phys., 8(1), 247–257. DOI:  10.26464/epp2023051
Anders
