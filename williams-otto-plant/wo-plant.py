from pyomo.environ import * 
from pyomo.dae import * 

# 
# Model 
# 
mo = AbstractModel()

# Parameter 
mo.n = Param() # 1x1
mo.m = Param() # 1x1
mo.I = RangeSet(1, mo.m) # mx1 
mo.J = RangeSet(1, mo.n) # nx1
mo.A = Param(mo.I) # mx1 
mo.Ea = Param(mo.I) # mx1
mo.W = Param() # 1x1  
mo.fb_cons = Param(RangeSet(1, 2)) # 2x1 
mo.Tr_cons = Param(RangeSet(1, 2)) # 2x1 
mo.t_range = Param(RangeSet(1, 2)) # 2x1 
mo.x_init = Param(mo.J) # nx1
mo.Tr_init = Param() # 1x1 
mo.fa_init = Param() # 1x1 
mo.fb_init = Param() # 1x1 

# Variable 
mo.t = ContinuousSet(bounds = (0, 50)) # tx1
mo.x = Var(mo.J, mo.t) # nxt 
mo.fb = Var(mo.t)
mo.fa = Var(mo.t)
mo.Tr = Var(mo.t)

# Note that we need to make sure fa, fb, Tr are functions of time
# Rate  
def k(am, I):
    return am.A[I] * exp(am.Ea[I] / am.Tr)

# Derivative 
mo.dx_dt = DerivativeVar(mo.x, wrt=mo.t) # nxt 

# Differential mass balances
def _xa_rule(am, t):
    if t == 0: 
        return Constraint.Skip
    return am.dx_dt[1, t] == (1 / am.W) * (am.fa - ((am.fa + am.fb) * am.x[1]) \
    - (k(am, 1) * am.x[1] * am.x[2] * am.W))
mo.xa_rule = Constraint(mo.t, rule=_xa_rule)
def _xb_rule(am, t):
    if t == 0: 
        return Constraint.Skip
    return am.dx_dt[2, t] == (1 / am.W) * (am.fb - ((am.fa + am.fb) * am.x[2]) \
    - (k(am, 1) * am.x[1] * am.x[2] * am.W) \
    - (k(am, 2) * am.x[2] * am.x[3] * am.W))
mo.xb_rule = Constraint(mo.t, rule=_xb_rule)
def _xc_rule(am, t): 
    if t == 0: 
        return Constraint.Skip
    return am.dx_dt[3, t] == (1 / am.W) * ((-1 * (am.fa + am.fb) * am.x[3]) \
    + 2 * (k(am, 1) * am.x[1] * am.x[2] * am.W) \
    - 2 * (k(am, 2) * am.x[2] * am.x[3] * am.W) \
    - (k(am, 3) * am.x[3] * am.x[6] * am.W))
mo.xc_rule = Constraint(mo.t, rule=_xc_rule)
def _xe_rule(am, t): 
    if t == 0: 
        return Constraint.Skip
    return am.dx_dt[4, t] == (1 / am.W) * ((-1 * (am.fa + am.fb) * am.x[4]) \
    + 2 * (k(am, 2) * am.x[2] * am.x[3] * am.W))
mo.xe_rule = Constraint(mo.t, rule=_xe_rule)
def _xg_rule(am, t): 
    if t == 0: 
        return Constraint.Skip
    return am.dx_dt[5, t] == (1 / am.W) * ((-1 * (am.fa + am.fb) * am.x[5]) \
    + 1.5 * (k(am, 3) * am.x[3] * am.x[6] * am.W))
mo.xg_rule = Constraint(mo.t, rule=_xg_rule)
def _xp_rule(am, t): 
    if t == 0: 
        return Constraint.Skip
    return am.dxp_dt[6, t] == (1 / am.W) * ((-1 * (am.fa + am.fb) * am.x[6]) \
    + (k(am, 2) * am.x[2] * am.x[3] * am.W) \
    - 0.5 * (k(am, 3) * am.x[3] * am.x[6] * am.W))
mo.xp_rule = Constraint(mo.t, rule=_xp_rule)


# Inequality constraints
def fb_rule(am):
    return inequality(am.fb_cons[1], am.fb, am.fb_cons[2])
mo.fb_rule = Constraint(rule=fb_rule)
def tr_rule(am):
    return inequality(am.Tr_cons[1], am.Tr, am.Tr_cons[2])
mo.tr_rule = Constraint(rule=tr_rule)

# Setting initial conditions 
def _initxi(am, i): 
    return (am.x[i, 0] == am.x_init[i])

mo.x_con = Constraint(mo.J, rule=_initxi)
# print(f"test{mo.fb_cons[1]}") # code breaks here
def _inittr(am):
    return am.Tr[am.T].fix(am.Tr_init)
mo.tr_con = Constraint(rule=_inittr)
def _initfa(am):
    return am.fa[am.T].fix(am.fa_init)
mo.fa_con = Constraint(rule=_initfa)
def _initfb(am):
    return am.fb[am.T].fix(am.fb_init) 
mo.fb_con = Constraint(rule=_initfb)


# discretizer = TransformationFactory('dae.finite_difference')
# discretizer.apply_to(mo, nfe=100, wrt=mo.t, scheme='FORWARD')
# print(f"test{mo.fb_cons[1]}") # code breaks here