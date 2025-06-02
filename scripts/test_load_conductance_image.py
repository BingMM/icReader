#%% Import

import matplotlib.pyplot as plt
from icreader import conductanceimage

#%% Paths

base = '/Home/siv32/mih008/repos/icReader/'
path_in = base + 'example_data/or_0085.nc'

#%% Load conductance Image

cI = conductanceimage(path_in)

#%% Plot something

plt.figure(figsize=(10,10))
plt.imshow(cI.H[50])
plt.colorbar()
