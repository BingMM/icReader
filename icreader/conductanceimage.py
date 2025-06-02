#%% Import

import numpy as np
from secsy import CSgrid, CSprojection
from scipy.io import netcdf_file
from datetime import datetime, timedelta

#%% Conductance Image class

class ConductanceImage:
    """
    Container for loading and accessing conductance data from a NetCDF file.

    This class handles precomputed quantities such as Hall and Pedersen conductances, 
    and associated uncertainties.

    Attributes
    ----------
    E0 : np.ndarray
        Characteristic energy [keV] for each pixel.
    dE0 : np.ndarray
        Uncertainty in E0.
    Fe : np.ndarray
        Electron energy flux [erg/cm^2/s].
    dFe : np.ndarray
        Uncertainty in Fe.
    R : np.ndarray
        WIC and SI13 ratio.
    dR : np.ndarray
        Uncertainty in R.
    P : np.ndarray
        Pedersen conductance [mho].
    H : np.ndarray
        Hall conductance [mho].
    dP : np.ndarray
        Uncertainty in P.
    dH : np.ndarray
        Uncertainty in H.
    dP2 : np.ndarray
        Secondary uncertainty or alternate error metric for P.
    dH2 : np.ndarray
        Secondary uncertainty or alternate error metric for H.
    Ep : float
        Proton characteristic energy [keV]
    dEp : float
        Uncertainty in Ep.
    shape : tuple
        Shape of the conductance arrays (typically [y, x]).
    time : np.ndarray, optional
        Array of datetimes corresponding to the image snapshots, if available.
    grid : CSgrid
        Cubbed sphere grid on which the data is projected.
    """
    def __init__(self, filename: str):
        """
        Load conductance image data from a NetCDF file.

        Parameters
        ----------
        filename : str
            Path to the NetCDF file containing the conductance and grid data.
        """
        with netcdf_file(filename, 'r') as nc:
            def load_var(name):
                return np.copy(nc.variables[name][:])

            for attr in ['wic_avg', 'wic_std', 's12_avg', 's12_std', 's13_avg', 
                         's13_std', 'E0', 'dE0', 'Fe', 'dFe', 'R', 'dR', 'P', 'H', 
                         'dP', 'dH', 'dP2', 'dH2']:
                setattr(self, attr, load_var(attr))

            self.Ep = float(nc.Ep)
            self.dEp = float(nc.dEp)
            self.shape = self.E0.shape

            if "time" in nc.variables:
                ref = datetime.strptime(nc.reference_time.decode(), "%Y-%m-%dT%H:%M:%S")
                self.time = np.array([ref + timedelta(seconds=int(s)) for s in nc.variables["time"][:]])

            self.grid = CSgrid(
                CSprojection(nc.position, nc.orientation),
                L=nc.L, W=nc.W, Lres=nc.Lres, Wres=nc.Wres, R=nc.gridR
                )
