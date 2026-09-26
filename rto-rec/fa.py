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
t_i = np.arange(0, t_f + dT, dT)

# Function 
t_infl = [50, 150]
fa_sp = [5, -5]
fa_profile = np.insert(fa_sp, 0, fa_amb)[np.searchsorted(t_infl, t_i, side='right')]
def fa(t):
    return dict(zip(t_i, fa_profile))





