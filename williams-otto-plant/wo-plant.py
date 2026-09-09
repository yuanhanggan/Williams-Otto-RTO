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
mo.fa = Var(mo.t)
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

# Differential mass balances
def _xa_rule(am, t):
    if t == 0: 
        return Constraint.Skip
    return am.dxa_dt[t] == (1 / am.w) * (am.fa - ((am.fa + am.fb) * am.x[1]) \
    - (k(am, 1) * am.x[1] * am.x[2] * am.W))
def _xb_rule(am, t):
    if t == 0: 
        return Constraint.Skip
    return am.dxb_dt[t] == (1 / am.w) * (am.fb - ((am.fa + am.fb) * am.x[2]) \
    - (k(am, 1) * am.x[1] * am.x[2] * am.W) \
    - (k(am, 2) * am.x[2] * am.x[3] * am.W))
def _xc_rule(am, t): 
    if t == 0: 
        return Constraint.Skip
    return am.dxc_dt[t] == (1 / am.w) * ((-1 * (am.fa + am.fb) * am.x[3]) \
    + 2 * (k(am, 1) * am.x[1] * am.x[2] * am.W) \
    - 2 * (k(am, 2) * am.x[2] * am.x[3] * am.W) \
    - (k(am, 3) * am.x[3] * am.x[6] * am.W))
def _xe_rule(am, t): 
    if t == 0: 
        return Constraint.Skip
    return am.dxe_dt[t] == (1 / am.w) * ((-1 * (am.fa + am.fb) * am.x[4]) \
    + 2 * (k(am, 2) * am.x[2] * am.x[3] * am.W))
def _xg_rule(am, t): 
    if t == 0: 
        return Constraint.Skip
    return am.dxg_dt[t] == (1 / am.w) * ((-1 * (am.fa + am.fb) * am.x[5]) \
    + 1.5 * (k(am, 3) * am.x[3] * am.x[6] * am.W))
def _xp_rule(am, t): 
    if t == 0: 
        return Constraint.Skip
    return am.dxp_dt[t] == (1 / am.w) * ((-1 * (am.fa + am.fb) * am.x[6]) \
    + (k(am, 2) * am.x[2] * am.x[3] * am.W) \
    - 0.5 * (k(am, 3) * am.x[3] * am.x[6] * am.W))

# Inequality constraints
def fb_rule(am):
    return inequality(am.fb_cons[1], am.fb, am.fb_cons[2])
def tr_rule(am):
    return inequality(am.Tr_cons[1], am.Tr, am.Tr_cons[2])

# Setting initial conditions 
mo.x[1, 0] = mo.x_init[1]
mo.x[2, 0] = mo.x_init[2]
mo.x[3, 0] = mo.x_init[3]
mo.x[4, 0] = mo.x_init[4]
mo.x[5, 0] = mo.x_init[5]
mo.x[6, 0] = mo.x_init[6]
mo.Tr[0] = mo.Tr_init
mo.fa[0] = mo.fa_init
mo.fb[0] = mo.fb_init