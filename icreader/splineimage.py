#%% Import

import numpy as np
from secsy import CSgrid, CSprojection
from scipy.io import netcdf_file
from datetime import datetime, timedelta
import scipy.sparse as sp
from sksparse.cholmod import cholesky

#%% Conductance Image class

class SplineImage:
    """
    Container for loading and accessing spline (conductance) data from a NetCDF file.

    This class handles precomputed quantities such as Hall and Pedersen conductances, 
    and associated uncertainties.

    Attributes
    ----------
    P : np.ndarray
        Pedersen conductance [mho].
    H : np.ndarray
        Hall conductance [mho].
    dP : np.ndarray
        Uncertainty in P.
    dH : np.ndarray
        Uncertainty in H.
    time : np.ndarray, optional
        Array of datetimes corresponding to the image snapshots, if available.
    grid : CSgrid
        Cubbed sphere grid on which the data is projected.
    mH : np.ndarray
        Hall spline model coefficients.
    mP : np.ndarray
        Pedersen spline model coefficients.
    G : np.ndarray
        2D spatial spline basis functions
    kt : int
        something
    nkt : int
        something
    LH : csc matrix
        Lower triangular matrix from Cholesky factorization of the Hall sensitivity matrix
    LP : csc matrix
        Lower triangular matrix from Cholesky factorization of the Pedersen sensitivity matrix
    PH : np.ndarray
        Permuation of LH
    PP : np.ndarray
        Permutation of LP
    """
    def __init__(self, filename: str):
        """
        Load spline image data from a NetCDF file.

        Parameters
        ----------
        filename : str
            Path to the NetCDF file containing the conductance and grid data.
        """
        self.load_file(filename)
        self.reset_time()
        self.reset_C()

#%% Load file

    def load_file(self, filename: str):
        with netcdf_file(filename, 'r') as nc:
            def load_var(name):
                return np.copy(nc.variables[name][:])
        
            for attr in ['H', 'P', 'dH', 'dP', 'mH', 'mP', 'G', 'kt', 'nkt', 'PH', 'PP']:
                setattr(self, attr, load_var(attr))
            
            data    = nc.variables["LH_data"][:].copy()
            indices = nc.variables["LH_indices"][:].copy()
            indptr  = nc.variables["LH_indptr"][:].copy()
            shape   = tuple(nc.LH_shape)  # stored as attribute (tuple of ints)
            self.LH = sp.csc_matrix((data, indices, indptr), shape=shape)
        
            data    = nc.variables["LP_data"][:].copy()
            indices = nc.variables["LP_indices"][:].copy()
            indptr  = nc.variables["LP_indptr"][:].copy()
            shape   = tuple(nc.LP_shape)  # stored as attribute (tuple of ints)
            self.LP = sp.csc_matrix((data, indices, indptr), shape=shape)
        
            if "time" in nc.variables:
                ref = datetime.strptime(nc.reference_time.decode(), "%Y-%m-%dT%H:%M:%S")
                self.time = np.array([ref + timedelta(seconds=int(s)) for s in nc.variables["time"][:]])
        
            self.grid = CSgrid(
                CSprojection(nc.position, nc.orientation),
                L=nc.L, W=nc.W, Lres=nc.Lres, Wres=nc.Wres, R=nc.gridR
                )

#%% Define internal time
    
    def reset_time(self):
        self._t = None
    
    @property
    def t(self):
        if self._t is None:
            self._t = np.array([self.get_t(t) for t in self.time])
        return self._t

    def get_t(self, dt: datetime):
        if (dt < self.time[0]) or (dt > self.time[-1]):
            raise ValueError('The requested date is outside the models temporal constraints.')
        
        return (dt - self.time[0]).total_seconds()

#%% Define posterior model covariance matrix
        
    def reset_C(self):
        self._CH = None
        self._CP = None
    
    @property
    def CH(self):
        if self._CH is None:
            self._CH = self.compute_covariance_cholmod(comp='hall')
        return self._CH
    
    @property
    def CP(self):
        if self._CP is None:
            self._CP = self.compute_covariance_cholmod(comp='pedersen')
        return self._CP
    
    def compute_covariance_cholmod(self, comp='hall'):
        """
        Recreate cholmod Factor object and use its inv() method
        """
        if comp == 'hall':
            L, P = self.LH, self.PH
        else:
            L, P = self.LP, self.PP
        
        # Create a minimal Factor-like object
        # Note: This is a bit hacky but works
        
        # First create a dummy factor to get the structure
        
        # Create A from L and P (we need this just once)
        P_inv = np.empty_like(P)
        P_inv[P] = np.arange(len(P))
        
        A_perm = L @ L.T  # This is P*A*P^T
        A = A_perm[P_inv, :][:, P_inv]  # Get back A
        
        # Now factor it to get a proper Factor object
        factor = cholesky(A)
        C = factor.inv()
        
        return C






















