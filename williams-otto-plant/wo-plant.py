from pyomo.environ import * 
from pyomo.dae import * 
import pandas as pd 
import os, datetime

# Model 
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
mo.t = ContinuousSet(bounds = (0, 3000)) # tx1
mo.x = Var(mo.J, mo.t, bounds = (0, 1)) # nxt 
mo.fb = Var(mo.t) # 1xt 
mo.fa = Var(mo.t) # 1xt
mo.Tr = Var(mo.t, initialize=300.0) # 1xt

# Rate  
def k(am, I, t):
    return am.A[I] * exp(am.Ea[I] / am.Tr[t])

# Derivative 
mo.dx_dt = DerivativeVar(mo.x, wrt=mo.t) # nxt 

# Differential mass balances 
def _xa_rule(am, t):
    if t == 0: 
        return Constraint.Skip
    return am.dx_dt[1, t] == (1 / am.W) * (am.fa[t] - ((am.fa[t] + am.fb[t]) * am.x[1, t]) \
    - (k(am, 1, t) * am.x[1, t] * am.x[2, t] * am.W))
mo.xa_rule = Constraint(mo.t, rule=_xa_rule)
def _xb_rule(am, t):
    if t == 0: 
        return Constraint.Skip
    return am.dx_dt[2, t] == (1 / am.W) * (am.fb[t] - ((am.fa[t] + am.fb[t]) * am.x[2, t]) \
    - (k(am, 1, t) * am.x[1, t] * am.x[2, t] * am.W) \
    - (k(am, 2, t) * am.x[2, t] * am.x[3, t] * am.W))
mo.xb_rule = Constraint(mo.t, rule=_xb_rule)
def _xc_rule(am, t): 
    if t == 0: 
        return Constraint.Skip
    return am.dx_dt[3, t] == (1 / am.W) * ((-1 * (am.fa[t] + am.fb[t]) * am.x[3, t]) \
    + 2 * (k(am, 1, t) * am.x[1, t] * am.x[2, t] * am.W) \
    - 2 * (k(am, 2, t) * am.x[2, t] * am.x[3, t] * am.W) \
    - (k(am, 3, t) * am.x[3, t] * am.x[6, t] * am.W))
mo.xc_rule = Constraint(mo.t, rule=_xc_rule)
def _xe_rule(am, t): 
    if t == 0: 
        return Constraint.Skip
    return am.dx_dt[4, t] == (1 / am.W) * ((-1 * (am.fa[t] + am.fb[t]) * am.x[4, t]) \
    + 2 * (k(am, 2, t) * am.x[2, t] * am.x[3, t] * am.W))
mo.xe_rule = Constraint(mo.t, rule=_xe_rule)
def _xg_rule(am, t): 
    if t == 0: 
        return Constraint.Skip
    return am.dx_dt[5, t] == (1 / am.W) * ((-1 * (am.fa[t] + am.fb[t]) * am.x[5, t]) \
    + 1.5 * (k(am, 3, t) * am.x[3, t] * am.x[6, t] * am.W))
mo.xg_rule = Constraint(mo.t, rule=_xg_rule)
def _xp_rule(am, t): 
    if t == 0: 
        return Constraint.Skip
    return am.dx_dt[6, t] == (1 / am.W) * ((-1 * (am.fa[t] + am.fb[t]) * am.x[6, t]) \
    + (k(am, 2, t) * am.x[2, t] * am.x[3, t] * am.W) \
    - 0.5 * (k(am, 3, t) * am.x[3, t] * am.x[6, t] * am.W))
mo.xp_rule = Constraint(mo.t, rule=_xp_rule)

# Inequality constraints
def fb_rule(am, t):
    return inequality(am.fb_cons[1], am.fb[t], am.fb_cons[2])
mo.fb_rule = Constraint(mo.t, rule=fb_rule)
def tr_rule(am, t):
    return inequality(am.Tr_cons[1], am.Tr[t], am.Tr_cons[2])
mo.tr_rule = Constraint(mo.t, rule=tr_rule)

# Setting initial conditions 
def _initxi(am, i): 
    return (am.x[i, 0] == am.x_init[i])
mo.x_con = Constraint(mo.J, rule=_initxi)
# def _initxf(am, i):
    # return (am.x[i, am.t.last()] == am.x_init[i])
# mo.x_conf = Constraint(mo.J, rule=_initxf)
def _inittr(am):
    return am.Tr[0] == am.Tr_init
mo.tr_con = Constraint(rule=_inittr)
def _initfa(am, i):
    return am.fa[i] == am.fa_init
mo.fa_con = Constraint(mo.t, rule=_initfa)
def _initfb(am, i):
    return am.fb[i] == am.fb_init
mo.fb_con = Constraint(mo.t, rule=_initfb)

# Run
sv_dir = os.path.join(os.getcwd(), "sims", datetime.datetime.now().strftime("%y%m%d%H%M"))
os.makedirs(sv_dir)
ist = mo.create_instance('wo-plant.dat')
discretizer = TransformationFactory('dae.finite_difference')
discretizer.apply_to(ist, nfe=100, wrt=ist.t, scheme='FORWARD')
ist.obj = Objective(expr=1)
solver = SolverFactory('ipopt')
solver.options['halt_on_ampl_error'] = 'yes'
results = solver.solve(ist, tee=True, keepfiles=True, logfile=os.path.join(sv_dir, "wo.log"))

# Saving results 

res = pd.DataFrame(index=list(ist.t))
res.index.name = 't'
for v in ist.component_objects(Var, active=True):
    subsets = list(v.index_set().subsets())
    if v.dim() == 1: 
        res[v.name] = [value(v[i]) for i in ist.t]
    elif v.dim() == 2: 
        J, _ = subsets
        for j in J:   
            res[f"{v.name}_{j}"] = [value(v[j, i]) for i in ist.t]
res.to_csv(os.path.join(sv_dir, "wo.csv"))
with open(os.path.join(sv_dir, "wo.txt") , 'w') as file:
    ist.pprint(ostream=file)


