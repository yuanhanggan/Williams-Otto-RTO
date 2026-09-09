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
# mo.fa = Param(within=PositiveReals) # 1x1 may change this, this is a disturbance variable 
mo.fb_cons = Param(RangeSet(1, 2), within=PositiveReals) # 2x1 
mo.Tr_cons = Param(RangeSet(1, 2), within=PositiveReals) # 2x1 

mo.t = ContinuousSet(bounds = (0, 50)) # (tf - t0)/h x 1
mo.x = Var(mo.J, mo.t)
mo.fb = Var(mo.t)
mo.fa = Var(mo.t)
mo.Tr = Var(mo.t)

x_init = Param(mo.J, initialize=(0.09576092651018556, 0.38800254944550255, 0.015975230722923693, 0.28492243697452413, 0.10931645678940315, 0.10602239955746103))
Tr_init = Param(initialize=365.9045240211792)
fa_init = Param(initialize=2.4)
fb_init = Param(initialize=6.091109369944852)

# Rate equation 
def k(am, I):
    return am.A[I] * exp(am.Ea[I] / am.Tr)

# Objective expression 
# def obj_rule(am):
#     return -1 * ((5554.1 * (am.fa + am.fb) * am.x[6]) \
#     + (125.91 * (am.fa + am.fb) * am.x[4]) \
#     - (370.3 * am.fa) \
#     - (555.42 * am.fb))

# Derivative expressions
mo.dxa_dt = DerivativeVar(mo.x[1], wrt=mo.t)
mo.dxb_dt = DerivativeVar(mo.x[2], wrt=mo.t)
mo.dxc_dt = DerivativeVar(mo.x[3], wrt=mo.t)
mo.dxe_dt = DerivativeVar(mo.x[4], wrt=mo.t)
mo.dxg_dt = DerivativeVar(mo.x[5], wrt=mo.t)
mo.dxp_dt = DerivativeVar(mo.x[6], wrt=mo.t)

discretizer = TransformationFactory('dae.finite_difference')
discretizer.apply(mo, nfe=60, wrt = mo.t, scheme='FORWARD')



# Differential mass balances
def _xa_rule(am, t):
    if t == 0: 
        return Constraint.Skip
    return am.dxa_dt[t] == (1 / am.w) * (am.fa - ((am.fa + am.fb) * am.x[1]) \
    - (k(am, 1) * am.x[1] * am.x[2] * am.W))
mo.xa_rule = Constraint(mo.t, rule=_xa_rule)
def _xb_rule(am, t):
    if t == 0: 
        return Constraint.Skip
    return am.dxb_dt[t] == (1 / am.w) * (am.fb - ((am.fa + am.fb) * am.x[2]) \
    - (k(am, 1) * am.x[1] * am.x[2] * am.W) \
    - (k(am, 2) * am.x[2] * am.x[3] * am.W))
mo.xb_rule = Constraint(mo.t, rule=_xb_rule)
def _xc_rule(am, t): 
    if t == 0: 
        return Constraint.Skip
    return am.dxc_dt[t] == (1 / am.w) * ((-1 * (am.fa + am.fb) * am.x[3]) \
    + 2 * (k(am, 1) * am.x[1] * am.x[2] * am.W) \
    - 2 * (k(am, 2) * am.x[2] * am.x[3] * am.W) \
    - (k(am, 3) * am.x[3] * am.x[6] * am.W))
mo.xc_rule = Constraint(mo.t, rule=_xc_rule)
def _xe_rule(am, t): 
    if t == 0: 
        return Constraint.Skip
    return am.dxe_dt[t] == (1 / am.w) * ((-1 * (am.fa + am.fb) * am.x[4]) \
    + 2 * (k(am, 2) * am.x[2] * am.x[3] * am.W))
mo.xe_rule = Constraint(mo.t, rule=_xe_rule)
def _xg_rule(am, t): 
    if t == 0: 
        return Constraint.Skip
    return am.dxg_dt[t] == (1 / am.w) * ((-1 * (am.fa + am.fb) * am.x[5]) \
    + 1.5 * (k(am, 3) * am.x[3] * am.x[6] * am.W))
mo.xg_rule = Constraint(mo.t, rule=_xg_rule)
def _xp_rule(am, t): 
    if t == 0: 
        return Constraint.Skip
    return am.dxp_dt[t] == (1 / am.w) * ((-1 * (am.fa + am.fb) * am.x[6]) \
    + (k(am, 2) * am.x[2] * am.x[3] * am.W) \
    - 0.5 * (k(am, 3) * am.x[3] * am.x[6] * am.W))
mo.xp_rule = Constraint(mo.t, rule=_xp_rule)
# Inequality constraints
def fb_rule(am):
    return inequality(am.fb_cons[1], am.fb, am.fb_cons[2])
def tr_rule(am):
    return inequality(am.Tr_cons[1], am.Tr, am.Tr_cons[2])
mo.fb_rule = Constraint(rule=fb_rule)
mo.tr_rule = Constraint(rule=tr_rule)
# Setting initial conditions 
def _initxi(am, i): 
    return am.x[i, 0] == am.x_init[i, 0]
mo.x_con = Constraint(mo.n, rule=_initxi)
def _inittr(am):
    return am.Tr[am.T].fix(am.Trinit)
mo.tr_con = Constraint(mo.T, rule=_inittr)
def _initfa(am):
    return am.fa[am.T].fix(am.fa_init)
mo.fa_con = Constraint(mo.T, rule=_initfa)
def _initfb(am):
    return am.fb[am.T].fix(am.fb_init) 
mo.fb_con = Constraint(mo.T, rule=_initfb)

# # Piecewise constant inputs 
# def _piecewisefa(am, t):
