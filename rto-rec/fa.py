# rto-rec/fa.py
# Defining and discretizing arbitrary function over t_f - 0
import os, numpy as np
from pyomo.dae import *
from pyomo.environ import *

# Data
dt, dT, t_f = int(os.environ.get('dt')), int(os.environ.get('dT')), int(os.environ.get('t_f'))
da = DataPortal()
da.load(filename=os.path.join(os.getcwd(), 'data', 'wo-rec.dat'))
fa_amb = da['fa_init']

# Arbitrary step function 
# t_infl = np.array([100, 200, 400, 500, 700, 800], dtype=float)
# fa_sp = np.array([-1.2, 1.2, -1.2, 1.2, -1.2, 1.2], dtype=float)
# t_fa = np.arange(0, t_f + dT, dT)
# fa_prof = np.insert(np.cumsum(fa_sp).astype(float) + fa_amb, 0, fa_amb)[np.searchsorted(t_infl, t_fa, side='right')]

# Arbitrary linear function 
t_fa = np.arange(0, t_f + dT, dT)
t_infl = [0, t_f + dT]
v_bound = [2, 10]
fa_prof =np.interp(t_fa, t_infl, v_bound)
