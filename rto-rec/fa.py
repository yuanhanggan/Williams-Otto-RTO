# rto-rec/fa.py
# Defining and discretizing arbitrary function over t_f - 0
import os, numpy as np
from pyomo.dae import *
from pyomo.environ import *

# Data
dt, dT, t_f = int(os.environ.get('dt')), int(os.environ.get('dT')), int(os.environ.get('t_f'))
da = DataPortal()
da.load(filename='wo-rec.dat')
fa_amb = da['fa_init']

# Function 
t_infl = np.array([50, 150], dtype=float)
fa_sp = np.array([5, -5], dtype=float)
fa_profile = np.insert(np.cumsum(fa_sp).astype(float) + fa_amb, 0, fa_amb)[np.searchsorted(t_infl, np.arange(0, t_f + dT, dT), side='right')]
def fa(t):
    return dict(zip(np.arange(0, t_f + dT, dT), fa_profile)).get(t)





