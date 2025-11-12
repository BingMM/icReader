#%% Import

import os
import numpy as np
import matplotlib.pyplot as plt
import icreader
from icreader import SplineImage

#%% Paths

# Get path to root of icreader repo
base = os.path.dirname(os.path.abspath(icreader.__file__))
base = os.path.abspath(os.path.join(base, '..'))  # Move up to repo root

# Paths
path_in = os.path.join(base, 'example_data', 'splineimages', 'or_0099.nc')
fig_out = os.path.join(base, 'figures', 'test_plot_sI.png')

#%% Load conductance Image

sI = SplineImage(path_in)

#%% Plot something

var = sI.P[51]
#var[var < 0] = np.nan

plt.ioff()
plt.figure(figsize=(10,10))
plt.imshow(var)
plt.colorbar()
plt.savefig(fig_out, bbox_inches='tight')
plt.close('all')
plt.ion()


#%%




