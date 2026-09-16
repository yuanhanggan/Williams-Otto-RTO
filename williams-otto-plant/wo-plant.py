from pyomo.environ import * 
from pyomo.dae import * 
import pandas as pd, os, datetime

# Model 
mo = ConcreteModel()
da = DataPortal()
da.load(filename='wo-plant.dat')

# Parameters 
mo.I = RangeSet(1, 3) # mx1 
mo.J = RangeSet(1, 6) # nx1
mo.A = Param(mo.I, initialize=da['A']) # mx1 
mo.Ea = Param(mo.I, initialize=da['Ea']) # mx1
mo.W = Param(initialize=da['W']) # 1x1  
mo.x_init = Param(mo.J, initialize=da['x_init']) # nx1
mo.Tr_init = Param(initialize=da['Tr_init']) # 1x1 
mo.fa_init = Param(initialize=da['fa_init']) # 1x1 
mo.fb_init = Param(initialize=da['fb_init']) # 1ddx1 
mo.q = Param(mo.J, initialize=0.01)

# Variables 
mo.t = ContinuousSet(bounds=(0, 1200)) # tx1
def init_x(am, j, t):
    return am.x_init[j]
mo.x = Var(mo.J, mo.t, bounds=(0, 1), initialize=init_x) # nxt 
mo.fb = Var(mo.t, bounds=(2, 10)) # 1xt 
mo.fa = Var(mo.t, initialize=mo.fa_init) # 1xt
mo.Tr = Var(mo.t, bounds=(323.15, 423.15)) # 1xt

# Rate  
def k(am, I, t):
    return am.A[I] * exp(am.Ea[I] / am.Tr[t])

# Objective 
def l2_rule(am):
    return sum(sum((am.q[i] * (am.x[i, t]  - am.x_init[i])) ** 2 for i in am.x_init) for t in am.t) # 1x1

# Derivative 
mo.dx_dt = DerivativeVar(mo.x, wrt=mo.t) # nxt 
# mo.dfa_dt = DerivativeVar(mo.fa, wrt=mo.t, initialize=0)

# Differential mass balances 
def _xa_rule(am, t):
    if t == am.t.first(): 
        return Constraint.Skip
    return am.dx_dt[1, t] == (1 / am.W) * (am.fa[t] - ((am.fa[t] + am.fb[t]) * am.x[1, t]) \
    - (k(am, 1, t) * am.x[1, t] * am.x[2, t] * am.W))
mo.xa_rule = Constraint(mo.t, rule=_xa_rule)
def _xb_rule(am, t):
    if t == am.t.first(): 
        return Constraint.Skip
    return am.dx_dt[2, t] == (1 / am.W) * (am.fb[t] - ((am.fa[t] + am.fb[t]) * am.x[2, t]) \
    - (k(am, 1, t) * am.x[1, t] * am.x[2, t] * am.W) \
    - (k(am, 2, t) * am.x[2, t] * am.x[3, t] * am.W))
mo.xb_rule = Constraint(mo.t, rule=_xb_rule)
def _xc_rule(am, t): 
    if t == am.t.first(): 
        return Constraint.Skip
    return am.dx_dt[3, t] == (1 / am.W) * ((-1 * (am.fa[t] + am.fb[t]) * am.x[3, t]) \
    + 2 * (k(am, 1, t) * am.x[1, t] * am.x[2, t] * am.W) \
    - 2 * (k(am, 2, t) * am.x[2, t] * am.x[3, t] * am.W) \
    - (k(am, 3, t) * am.x[3, t] * am.x[6, t] * am.W))
mo.xc_rule = Constraint(mo.t, rule=_xc_rule)
def _xe_rule(am, t): 
    if t == am.t.first(): 
        return Constraint.Skip
    return am.dx_dt[4, t] == (1 / am.W) * ((-1 * (am.fa[t] + am.fb[t]) * am.x[4, t]) \
    + 2 * (k(am, 2, t) * am.x[2, t] * am.x[3, t] * am.W))
mo.xe_rule = Constraint(mo.t, rule=_xe_rule)
def _xg_rule(am, t): 
    if t == am.t.first(): 
        return Constraint.Skip
    return am.dx_dt[5, t] == (1 / am.W) * ((-1 * (am.fa[t] + am.fb[t]) * am.x[5, t]) \
    + 1.5 * (k(am, 3, t) * am.x[3, t] * am.x[6, t] * am.W))
mo.xg_rule = Constraint(mo.t, rule=_xg_rule)
def _xp_rule(am, t): 
    if t == am.t.first(): 
        return Constraint.Skip
    return am.dx_dt[6, t] == (1 / am.W) * ((-1 * (am.fa[t] + am.fb[t]) * am.x[6, t]) \
    + (k(am, 2, t) * am.x[2, t] * am.x[3, t] * am.W) \
    - 0.5 * (k(am, 3, t) * am.x[3, t] * am.x[6, t] * am.W))
mo.xp_rule = Constraint(mo.t, rule=_xp_rule)

# Boundary conditions 
for j in mo.J:
    mo.x[j, 0].fix(mo.x_init[j])
    mo.dx_dt[j, 0].fix(0)
    # mo.x[j, mo.t.last()].fix(mo.x_init[j])
mo.Tr[0].fix(mo.Tr_init)
def _initfa(am, i):
    if i >= 120:
        return am.fa[i]==1.2
    else:
        return am.fa[i]==am.fa_init
mo.fa_con = Constraint(mo.t, rule=_initfa)
# def _initfb(am, i):
    # return am.fb[i]==am.fb_init
# mo.fb_con = Constraint(mo.t, rule=_initfb)
mo.fb[0].fix(mo.fb_init)
# Run
sv_dir = os.path.join(os.getcwd(), 'sims', datetime.datetime.now().strftime('%y%m%d%H%M'))
os.makedirs(sv_dir)
# dis = TransformationFactory('dae.collocation')
# dis.apply_to(mo, nfe=20, ncp=2, scheme='LAGRANGE-RADAU')
dis = TransformationFactory('dae.finite_difference')
dis.apply_to(mo, nfe=100, wrt=mo.t, scheme='BACKWARD')
mo.obj = Objective(rule=l2_rule, sense=minimize)
solver = SolverFactory('ipopt')
solver.options['halt_on_ampl_error'] = 'yes'
results = solver.solve(mo, tee=True, keepfiles=True, logfile = os.path.join(sv_dir, 'wo.log'))

# Saving results 
res = pd.DataFrame(index=list(mo.t))
bounds = {}
res.index.name = 't'
for v in mo.component_objects(Var, active = True):
    subsets = list(v.index_set().subsets())
    if v.dim() == 1: 
        res[v.name] = [value(v[i]) for i in mo.t]
        var = v[mo.t.first()]
        bounds[v.name] = (var.lb, var.ub)
    elif v.dim() == 2: 
        J, _ = subsets
        for j in J:   
            res[f'{v.name}_{j}'] = [value(v[j, i]) for i in mo.t]
            var = v[j, mo.t.first()]
            bounds[f'{v.name}_{j}'] = (var.lb, var.ub) 
res.to_csv(os.path.join(sv_dir, 'wo.csv'))
with open(os.path.join(sv_dir, 'wo.txt') , 'w') as file:
    mo.pprint(ostream = file)
pd.DataFrame.from_dict(bounds, orient='index', columns=['lower', 'upper']).to_csv(os.path.join(sv_dir, 'bounds.csv'))

