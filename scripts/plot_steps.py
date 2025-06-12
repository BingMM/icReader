#%% Import

import os
import matplotlib.pyplot as plt
import icreader
from icreader import ConductanceImage
from polplot import pp

#%% Paths

# Get path to root of icreader repo
base = os.path.dirname(os.path.abspath(icreader.__file__))
base = os.path.abspath(os.path.join(base, '..'))  # Move up to repo root

# Paths
path_in = os.path.join(base, 'example_data', 'or_0085.nc')
fig_out = os.path.join(base, 'figures/')

#%% Load conductance Image

cI = ConductanceImage(path_in)

#%% Step 1 : Binning

lt = (cI.grid.lon/15)%24
lat = cI.grid.lat

vmax = [4000, 12, 25,
         800, 12, 25]

plt.ioff()
fig, axs = plt.subplots(2, 3, figsize=(15, 10))
plt.subplots_adjust(hspace=0, wspace=0)
paxs = [pp(ax) for ax in axs.flatten()]

im = paxs[0].plotimg(lat, lt, cI.wic_avg[0], crange=(0,vmax[0]))
paxs[1].plotimg(lat, lt, cI.s13_avg[0], crange=(0,vmax[1]))
paxs[2].plotimg(lat, lt, cI.s12_avg[0], crange=(0,vmax[2]))

paxs[3].plotimg(lat, lt, cI.wic_std[0], crange=(0,vmax[3]))
paxs[4].plotimg(lat, lt, cI.s13_std[0], crange=(0,vmax[4]))
paxs[5].plotimg(lat, lt, cI.s12_std[0], crange=(0,vmax[5]))

for ax, title in zip(axs[0, :], ['WIC', 'SI13', 'SI12']):
    ax.text(.5, 1, title, ha='center', va='center', fontsize=16, transform=ax.transAxes)

for ax, label in zip(axs[:, 0], ['Median', 'Std']):
    ax.text(0, .5, label, ha='right', va='center', fontsize=16, transform=ax.transAxes, rotation='vertical')

for ax in axs[-1, :]:
    ax.text(.5, 0, 'MLT=00', ha='center', va='center', fontsize=12, transform=ax.transAxes)

for ax, vmax_ in zip(axs.flatten(), vmax):
    ax.text(.9, .85, f'max={vmax_}', ha='center', va='center', fontsize=12, transform=ax.transAxes)

ax = axs[1,2]
ax.text(.85, .15, '50$^{\circ}$', ha='center', va='center', fontsize=12, transform=ax.transAxes)

cbar = fig.colorbar(im, ax=axs, orientation='horizontal', fraction=0.04, pad=0.04)
cbar.set_ticks([0, .5*vmax[0], vmax[0]])
cbar.set_ticklabels(['0', '.5', 'max'], fontsize=16)

plt.suptitle('Step 1: Binned statistics', fontsize=25, y=.95)

plt.savefig(fig_out + 'step1.png', bbox_inches='tight')
plt.close('all')
plt.ion()

#%% Step 2 : Frey

lt = (cI.grid.lon/15)%24
lat = cI.grid.lat

vmax = [150, 25, 50,
        5000, 30, 10]

plt.ioff()
fig, axs = plt.subplots(2, 3, figsize=(15, 10))
plt.subplots_adjust(hspace=0, wspace=0)
paxs = [pp(ax) for ax in axs.flatten()]

im = paxs[0].plotimg(lat, lt, cI.R[0], crange=(0,vmax[0]))
paxs[1].plotimg(lat, lt, cI.E0[0], crange=(0,vmax[1]))
paxs[2].plotimg(lat, lt, cI.Fe[0], crange=(0,vmax[2]))

paxs[3].plotimg(lat, lt, cI.dR[0], crange=(0,vmax[3]))
paxs[4].plotimg(lat, lt, cI.dE0[0], crange=(0,vmax[4]))
paxs[5].plotimg(lat, lt, cI.dFe[0], crange=(0,vmax[5]))

for ax, title in zip(axs[0, :], ['R (WIC*/SI13*)', 'E0', 'Fe']):
    ax.text(.5, 1, title, ha='center', va='center', fontsize=16, transform=ax.transAxes)

for ax, label in zip(axs[:, 0], ['Mean', 'Std']):
    ax.text(0, .5, label, ha='right', va='center', fontsize=16, transform=ax.transAxes, rotation='vertical')

for ax in axs[-1, :]:
    ax.text(.5, 0, 'MLT=00', ha='center', va='center', fontsize=12, transform=ax.transAxes)

for ax, vmax_ in zip(axs.flatten(), vmax):
    ax.text(.9, .85, f'max={vmax_}', ha='center', va='center', fontsize=12, transform=ax.transAxes)

ax = axs[1,2]
ax.text(.85, .15, '50$^{\circ}$', ha='center', va='center', fontsize=12, transform=ax.transAxes)

cbar = fig.colorbar(im, ax=axs, orientation='horizontal', fraction=0.04, pad=0.04)
cbar.set_ticks([0, .5*vmax[0], vmax[0]])
cbar.set_ticklabels(['0', '.5', 'max'], fontsize=16)

plt.suptitle('Step 2: Frey et al. 2003', fontsize=25, y=.95)

plt.savefig(fig_out + 'step2.png', bbox_inches='tight')
plt.close('all')
plt.ion()

#%% Step 3 : Robinson

lt = (cI.grid.lon/15)%24
lat = cI.grid.lat

vmax = [ 70, 10, 1,
        150, 150]

plt.ioff()
fig, axs = plt.subplots(2, 3, figsize=(15, 10))
plt.subplots_adjust(hspace=0, wspace=0)
paxs = [pp(ax) for ax in axs.flatten()[:-1]]

im = paxs[0].plotimg(lat, lt, cI.H[0], crange=(0,vmax[0]))
paxs[1].plotimg(lat, lt, cI.P[0], crange=(0,vmax[1]))
paxs[2].plotimg(lat, lt, 1-cI.w[0], crange=(0,vmax[2]))

paxs[3].plotimg(lat, lt, cI.dH[0], crange=(0,vmax[3]))
paxs[4].plotimg(lat, lt, cI.dP[0], crange=(0,vmax[4]))
axs[1,-1].axis('off')

for ax, title in zip(axs[0, :], ['Hall', 'Pedersen', 'Weight']):
    ax.text(.5, 1, title, ha='center', va='center', fontsize=16, transform=ax.transAxes)

for ax, label in zip(axs[:, 0], ['Mean', 'Std']):
    ax.text(0, .5, label, ha='right', va='center', fontsize=16, transform=ax.transAxes, rotation='vertical')

for ax in [axs[1,0], axs[1,1], axs[0,2]]:
    ax.text(.5, 0, 'MLT=00', ha='center', va='center', fontsize=12, transform=ax.transAxes)

for ax, vmax_ in zip(axs.flatten()[:-1], vmax):
    ax.text(.9, .85, f'max={vmax_}', ha='center', va='center', fontsize=12, transform=ax.transAxes)

ax = axs[0,2]
ax.text(.85, .15, '50$^{\circ}$', ha='center', va='center', fontsize=12, transform=ax.transAxes)
ax = axs[1,2]
ax.text(.15, .6, 'Weights determined from\nIRLS fit of dayglow\nsubtraction and spherical\nharmonic correction.', ha='left', va='center', fontsize=20, transform=ax.transAxes)

cbar = fig.colorbar(im, ax=axs, orientation='horizontal', fraction=0.04, pad=0.04)
cbar.set_ticks([0, .5*vmax[0], vmax[0]])
cbar.set_ticklabels(['0', '.5', 'max'], fontsize=16)

plt.suptitle('Step 3: Robinson et al. 1987', fontsize=25, y=.95)

plt.savefig(fig_out + 'step3.png', bbox_inches='tight')
plt.close('all')
plt.ion()

