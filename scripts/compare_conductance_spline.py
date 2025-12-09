#%% Import

import os
import numpy as np
from datetime import timedelta
import matplotlib.pyplot as plt
import icreader
from icreader import ConductanceImage, SplineImage
from secsy import CSgrid, CSprojection
from polplot import Polarplot

#%% Paths

# Get path to root of icreader repo
base = os.path.dirname(os.path.abspath(icreader.__file__))
base = os.path.abspath(os.path.join(base, '..'))  # Move up to repo root

# Paths
orbit = 150
path_conductance = os.path.join(base, 'example_data', 'conductance', f'or_{str(orbit).zfill(4)}.nc')
path_spline = os.path.join(base, 'example_data', 'spline', f'or_{str(orbit).zfill(4)}.nc')
path_out = '/home/bing/Dropbox/work/temp_storage/spline_comp_150/'

#%% Load images

cI = ConductanceImage(path_conductance)
sI = SplineImage(path_spline)

#%% Fine grid

position = (0, 90) # lon, lat
orientation = (0, 1) # east, north
L, Lres = 20000e3, 50e3
grid = CSgrid(CSprojection(position, orientation), L, L, Lres, Lres, R = 6481.2e3)

t_res = 30 # seconds
total_time = (sI.time_[-1] - sI.ref).total_seconds()
time = np.array([sI.ref + timedelta(seconds=int(s)) for s in np.arange(0, total_time, t_res)])

#%% 

sI.set_space(x=grid.xi, y=grid.eta)

for i, time_ in enumerate(time):
    
    ii = np.argmin(abs(cI.time - time_))
    sI.set_time(time_)
    fig, axs = plt.subplots(2, 2, figsize=(15,15))
    paxs = [Polarplot(ax) for ax in axs.flatten()]
    f = ~np.isnan(cI.H[ii])
    cc0 = paxs[0].tricontourf(cI.grid.lat[f], cI.grid.lon[f]/15%24, cI.H[ii][f], levels=np.linspace(0, 70, 100))
    f = ~np.isnan(cI.dH[ii])
    cc2 = paxs[2].tricontourf(cI.grid.lat[f], cI.grid.lon[f]/15%24, cI.dH[ii][f], levels=np.linspace(0, 70, 100))
    cc1 = paxs[1].tricontourf(grid.lat.flatten(), grid.lon.flatten()/15%24, sI.H.flatten(), levels=np.linspace(0, 70, 100))
    cc3 = paxs[3].tricontourf(grid.lat.flatten(), grid.lon.flatten()/15%24, sI.dH.flatten(), levels=np.linspace(0, 5, 100))
    
    cbar = fig.colorbar(cc3, ax=axs.ravel().tolist(),
                        orientation='horizontal',
                        fraction=0.03, pad=0.05)
    cbar.set_ticks([0, 2.5, 5])
    cbar.set_ticklabels(['0', '1/2 max', 'max'], fontsize=15)
    cbar.set_label('[mho]', fontsize=20)
    
    for ax, tit, s in zip(axs.flatten(), 
                          ['Hall', 'Hall model', 'Hall uncertainty', 'Hall uncertainty model'],
                          [70, 70, 70, 5]):
        ax.set_axis_off()
        ax.set_title(tit, fontsize=25)
        ax.text(0.1, .9, f'max={s}', fontsize=15, ha='left', va='center', transform=ax.transAxes)
    
    for ax in axs[1, :]:
        ax.text(0.5, 0, 'mlt=0', fontsize=15, ha='center', va='center', transform=ax.transAxes)
    
    axs[0,0].text(0.85, .15, 'mlat=50', fontsize=15, ha='center', va='center', transform=axs[0,0].transAxes)
    
    plt.suptitle(time_, fontsize=35)
    
    plt.savefig(f'/home/bing/Dropbox/work/temp_storage/spline_comp_150/{str(i).zfill(5)}.png', bbox_inches='tight', dpi=150)
    plt.close('all')
                       
