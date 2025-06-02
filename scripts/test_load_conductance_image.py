#%% Import

import matplotlib.pyplot as plt
from icreader import ConductanceImage

#%% Paths

base = '/Home/siv32/mih008/repos/icReader/'
path_in = base + 'example_data/or_0085.nc'

#%% Load conductance Image

cI = ConductanceImage(path_in)

#%% Plot something

plt.ioff()
plt.figure(figsize=(10,10))
plt.imshow(cI.H[70])
plt.colorbar()
plt.savefig(base + 'figures/test_plot.png', bbox_inches='tight')
plt.close('all')
plt.ion()
