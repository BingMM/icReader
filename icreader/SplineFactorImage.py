#%% Import

import numpy as np
from secsy import CSgrid, CSprojection
from scipy.io import netcdf_file
from datetime import datetime, timedelta
import scipy.sparse as sp
from sksparse.cholmod import cholesky
from scipy.interpolate import BSpline
from scipy.sparse import kron, vstack, csc_matrix

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
        self.reset_spatial_evaluation()

#%% Load file

    def load_file(self, filename: str):
        with netcdf_file(filename, 'r') as nc:
            
            self._H = nc.variables['H'][:].copy()
            self._P = nc.variables['P'][:].copy()
            self._dH = nc.variables['dH'][:].copy()
            self._dP = nc.variables['dP'][:].copy()
            
            self.mH = nc.variables['mH'][:].copy()
            self.mP = nc.variables['mP'][:].copy()
            
            self.PH = nc.variables['PH'][:].copy()
            self.PP = nc.variables['PP'][:].copy()
            
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

            self.k, self.nk = nc.k, nc.nk
            self.kt, self.nkt = nc.kt, nc.nkt

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

#%% Splines
    @property
    def knots(self):
        return np.r_[[self.x.min()]*self.k, np.linspace(self.x.min(), self.x.max(), self.nk-2*self.k), [self.x.max()]*self.k]

    @property
    def tknots(self):
        return np.r_[[self.t.min()]*self.kt, np.linspace(self.t.min(), self.t.max(), self.nkt-2*self.kt), [self.t.max()]*self.kt]

#%% Define posterior model covariance matrix
        
    def reset_C(self):
        self._CH = None
        self._CP = None
    
    @property
    def CH(self):
        if self._CH is None:
            print('Forming Hall posterior model covariance matrix. This can take a while.')
            self._CH = self.compute_covariance_cholmod(comp='hall')
        return self._CH
    
    @property
    def CP(self):
        if self._CP is None:
            print('Forming Pedersen posterior model covariance matrix. This can take a while.')
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

#%% Spatial grid

    def reset_spatial_evaluation(self):
        self._x = None
        self._y = None
        self._G = None
    
    def set_spatial_evaluation(self, lon, lat):
        self.reset_spatial_evaluation()
        self._x, self._y = self.grid.projection.geo2cube(lon, lat)        

    @property
    def x(self):
        if self._x is None:
            return self.grid.xi
        else:
            return self._x
    
    @property
    def y(self):
        if self._y is None:
            return self.grid.eta
        else:
            return self._y

    @property
    def G(self):
        if self._G is None:
            self._G = self.generate_G_2d(return_sparse=True)
        return self._G

    def generate_G_2d(self, return_sparse=False):
        G = np.zeros((self.x.size, self.ncp**2))
        for i, (xi, yi) in enumerate(zip(self.x.flatten(), self.y.flatten())):
            Gx = BSpline.design_matrix(xi, self.knots, self.k).todense()
            Gy = BSpline.design_matrix(yi, self.knots, self.k).todense()
            Gy = np.kron(np.eye(self.ncp), Gy)
            G[i, :] = Gx.dot(Gy)
        if return_sparse:
            return csc_matrix(G)
        return G
        
#%% Evaluate

    @property
    def Gt(self):
        if self._Gt is None:
            self._Gt = self.generate_G_3d()
        return self._Gt
        
    def generate_G_3d(self):
        Gt = []
        for i, ti in enumerate(self.t):
            Gt_ = BSpline.design_matrix(ti, self.tknots, self.kt)  # Already sparse
            Gt_ = kron(np.eye(self.G.shape[1]), Gt_, format='csr')  # Ensure sparse output
            Gt.append(self.G @ Gt_)  # Matrix multiplication remains sparse
        return vstack(Gt, format='csr')  # Efficient sparse stacking






















