#%% Import

import os
import matplotlib.pyplot as plt
import icreader
from icreader import ConductanceImage

#%% Paths

# Get path to root of icreader repo
base = os.path.dirname(os.path.abspath(icreader.__file__))
base = os.path.abspath(os.path.join(base, '..'))  # Move up to repo root

# Paths
path_in = os.path.join(base, 'example_data', 'conductanceimages', 'or_0099.nc')
fig_out = os.path.join(base, 'figures', 'test_plot_cI.png')

#%% Load conductance Image

cI = ConductanceImage(path_in)

#%% Plot something

plt.ioff()
plt.figure(figsize=(10,10))
plt.imshow(cI.P[100])
plt.colorbar()
plt.savefig(fig_out, bbox_inches='tight')
plt.close('all')
plt.ion()
