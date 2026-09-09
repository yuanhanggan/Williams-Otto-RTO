from pyomo.environ import * 
from pyomo.dae import * 

# 
# Model 
# 

mo = AbstractModel()
mo.n = Param(within=NonNegativeIntegers) # 1x1
mo.m = Param(within=NonNegativeIntegers) # 1x1
mo.I = RangeSet(1, mo.m) # mx1 
mo.J = RangeSet(1, mo.n) # nx1
mo.A = Param(mo.I) # mx1 
mo.Ea = Param(mo.I) # mx1
mo.W = Param(within=PositiveReals) # 1x1 
mo.fa = Param(within=PositiveReals) # 1x1 may change this, this is a disturbance variable 
mo.fb_cons = Param(RangeSet(1, 2), within=PositiveReals) # 2x1
mo.Tr_cons = Param(RangeSet(1, 2), within=PositiveReals) # 2x1 

mo.t = ContinuousSet() # (tf - t0)/h x 1
mo.x = Var(mo.J, mo.t)
mo.fb = Var(mo.t)
mo.Tr = Var(mo.t)

# Rate equation 
def k(am, I):
    return am.A[I] * exp(am.Ea[I] / am.Tr)

# Objective expression 
def obj_rule(am):
    return -1 * ((5554.1 * (am.fa + am.fb) * am.x[6]) \
    + (125.91 * (am.fa + am.fb) * am.x[4]) \
    - (370.3 * am.fa) \
    - (555.42 * am.fb))

# Derivative expressions
mo.dxa_dt = DerivativeVar(mo.x[1], wrt=mo.t)
mo.dxb_dt = DerivativeVar(mo.x[2], wrt=mo.t)
mo.dxc_dt = DerivativeVar(mo.x[3], wrt=mo.t)
mo.dxe_dt = DerivativeVar(mo.x[4], wrt=mo.t)
mo.dxg_dt = DerivativeVar(mo.x[5], wrt=mo.t)
mo.dxp_dt = DerivativeVar(mo.x[6], wrt=mo.t)



