#%% Import

import os
import numpy as np
import matplotlib.pyplot as plt
import icreader
from icreader import SplineImage
from polplot import Polarplot

#%% Paths

# Get path to root of icreader repo
base = os.path.dirname(os.path.abspath(icreader.__file__))
base = os.path.abspath(os.path.join(base, '..'))  # Move up to repo root

# Paths
orbit = 150
path_in = os.path.join(base, 'example_data', 'spline', f'or_{str(orbit).zfill(4)}.nc')
fig1_out = os.path.join(base, 'figures', 'test_plot_sI.png')
fig2_out = os.path.join(base, 'figures', 'test_plot_sI_fun.png')

#%% Load conductance Image

sI = SplineImage(path_in)

#%% define grid

coordinate_system = 'geo'

lat = np.arange(40, 90)
lon = np.arange(0, 360)
lat, lon = np.meshgrid(lat, lon)

#%% Set space

sI.set_space(lon=lon, lat=lat, coord_sys=coordinate_system)

#%% Set time

time = sI.time_[101]
sI.set_time(time)

#%% Plot by calling variables directly

plt.ioff()
fig, axs = plt.subplots(2, 2, figsize=(15, 15))
paxs = [Polarplot(ax, minlat=40) for ax in axs.flatten()]

im0 = paxs[0].plotimg(lat, lon/15%24, sI.H)
im1 = paxs[1].plotimg(lat, lon/15%24, sI.P)
im2 = paxs[2].plotimg(lat, lon/15%24, sI.dH)
im3 = paxs[3].plotimg(lat, lon/15%24, sI.dP)

for ax, im in zip(axs.flatten(), [im0, im1, im2, im3]):
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_axis_off()

plt.savefig(fig1_out, bbox_inches='tight')
plt.close('all')
plt.ion()

#%% Create Lompe conductance function

H_fun  = sI.get_H_fun(coord_sys=coordinate_system)
P_fun  = sI.get_P_fun(coord_sys=coordinate_system)
dH_fun = sI.get_dH_fun(coord_sys=coordinate_system)
dP_fun = sI.get_dP_fun(coord_sys=coordinate_system)

#%% Plot by calling Lompe conductance functions

# Time can be changed, the conductance functions with update automatically
sI.set_time(sI.time_[100])

plt.ioff()
fig, axs = plt.subplots(2, 2, figsize=(15, 15))
paxs = [Polarplot(ax, minlat=40) for ax in axs.flatten()]

im0 = paxs[0].plotimg(lat, lon/15%24, H_fun(lon, lat))
im1 = paxs[1].plotimg(lat, lon/15%24, P_fun(lon, lat))
im2 = paxs[2].plotimg(lat, lon/15%24, dH_fun(lon, lat))
im3 = paxs[3].plotimg(lat, lon/15%24, dP_fun(lon, lat))

for ax, im in zip(axs.flatten(), [im0, im1, im2, im3]):
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_axis_off()

plt.savefig(fig2_out, bbox_inches='tight')
plt.close('all')
plt.ion()