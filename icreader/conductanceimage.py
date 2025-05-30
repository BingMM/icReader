#%% Import

import numpy as np
from secsy import CSgrid, CSprojection
from scipy.io import netcdf_file
from datetime import datetime, timedelta

#%% Conductance Image class

class ConductanceImage:
    def __init__(self, filename: str):
        with netcdf_file(filename, 'r') as nc:
            def load_var(name):
                return np.copy(nc.variables[name][:])

            for attr in ["E0", "dE0", "Fe", "dFe", "R", "dR", "P", "H", "dP", "dH", "dP2", "dH2"]:
                setattr(self, attr, load_var(attr))

            self.Ep = float(nc.Ep)
            self.dEp = float(nc.dEp)
            self.shape = self.E0.shape

            if "time" in nc.variables:
                ref = datetime.strptime(nc.reference_time.decode(), "%Y-%m-%dT%H:%M:%S")
                self.time = np.array([ref + timedelta(seconds=int(s)) for s in nc.variables["time"][:]])

            self.grid = CSgrid(
                CSprojection(nc.position, nc.orientation),
                nc.L, nc.W, nc.Lres, nc.Wres
                )
