from pyomo.environ import * 
from pyomo.dae import * 
import pandas as pd, os

# Model 
mo = ConcreteModel()
da = DataPortal()
da.load(filename='wo-rec.dat')
sv_dir = os.path.join(os.getcwd(), 'sims', 'rec_10_fwd')
da_r = pd.read_csv(os.path.join(sv_dir, 'wo.csv')).set_index('t').dropna().to_dict(orient='index')
t_is = list(da_r.keys())
x_is = ['x_1','x_2','x_3','x_4','x_5','x_6']

# Parameters 
mo.dt = Param(initialize=da['dt']) # 1x1 
mo.x_is = Set(initialize=x_is) # nx1
mo.A = Param(RangeSet(1, 3), initialize=da['A']) # mx1 
mo.Ea = Param(RangeSet(1, 3), initialize=da['Ea']) # mx1
mo.W = Param(initialize=da['W']) # 1x1  
mo.q = Param(mo.x_is, initialize=0.0035) # 1x1 
mo.t = ContinuousSet(bounds=(t_is[-1], t_is[-1]+mo.dt))
mo.x_i0 = Param(mo.x_is, initialize={k: da_r[t_is[-1]][k] for k in da_r[t_is[-1]] & mo.x_is}) # nx1 
mo.Tr_0 = Param(initialize=da_r[t_is[-1]]['Tr']) # 1x1 
mo.fa_0 = Param(initialize=da_r[t_is[-1]]['fa']) # 1x1 
mo.fb_0 = Param(initialize=da_r[t_is[-1]]['fb']) # 1x1  

# Variables 
def init_x(am, j, t):
    return am.x_i0[j]
def init_fb(am, t):
    return am.fb_0
def init_Tr(am, t):
    return am.Tr_0
mo.x = Var(mo.x_is, mo.t, bounds=(0, 1), initialize=init_x) # nxt  
mo.fb = Var(mo.t, bounds=(2, 10), initialize=init_fb) # 1xt 
mo.Tr = Var(mo.t, bounds=(323.15, 423.15), initialize=init_Tr) # 1xt 
mo.fa = Var(mo.t, initialize=mo.fa_0) # 1xt 
mo.dx_dt = DerivativeVar(mo.x, wrt=mo.t, initialize=0) # nxt 

# Rate  
def k(am, I, t):
    return am.A[I] * exp(am.Ea[I] / am.Tr[t])

# Differential mass balances 
def _xa_rule(am, t):
    if t == am.t.first(): 
        return Constraint.Skip
    return am.dx_dt['x_1', t] == (1 / am.W) * (am.fa[t] - ((am.fa[t] + am.fb[t]) * am.x['x_1', t]) \
    - (k(am, 1, t) * am.x['x_1', t] * am.x['x_2', t] * am.W))
mo.xa_rule = Constraint(mo.t, rule=_xa_rule)
def _xb_rule(am, t):
    if t == am.t.first(): 
        return Constraint.Skip
    return am.dx_dt['x_2', t] == (1 / am.W) * (am.fb[t] - ((am.fa[t] + am.fb[t]) * am.x['x_2', t]) \
    - (k(am, 1, t) * am.x['x_1', t] * am.x['x_2', t] * am.W) \
    - (k(am, 2, t) * am.x['x_2', t] * am.x['x_3', t] * am.W))
mo.xb_rule = Constraint(mo.t, rule=_xb_rule)
def _xc_rule(am, t): 
    if t == am.t.first(): 
        return Constraint.Skip
    return am.dx_dt['x_3', t] == (1 / am.W) * ((-1 * (am.fa[t] + am.fb[t]) * am.x['x_3', t]) \
    + 2 * (k(am, 1, t) * am.x['x_1', t] * am.x['x_2', t] * am.W) \
    - 2 * (k(am, 2, t) * am.x['x_2', t] * am.x['x_3', t] * am.W) \
    - (k(am, 3, t) * am.x['x_3', t] * am.x['x_6', t] * am.W))
mo.xc_rule = Constraint(mo.t, rule=_xc_rule)
def _xe_rule(am, t): 
    if t == am.t.first(): 
        return Constraint.Skip
    return am.dx_dt['x_4', t] == (1 / am.W) * ((-1 * (am.fa[t] + am.fb[t]) * am.x['x_4', t]) \
    + 2 * (k(am, 2, t) * am.x['x_2', t] * am.x['x_3', t] * am.W))
mo.xe_rule = Constraint(mo.t, rule=_xe_rule)
def _xg_rule(am, t): 
    if t == am.t.first(): 
        return Constraint.Skip
    return am.dx_dt['x_5', t] == (1 / am.W) * ((-1 * (am.fa[t] + am.fb[t]) * am.x['x_5', t]) \
    + 1.5 * (k(am, 3, t) * am.x['x_3', t] * am.x['x_6', t] * am.W))
mo.xg_rule = Constraint(mo.t, rule=_xg_rule)
def _xp_rule(am, t): 
    if t == am.t.first(): 
        return Constraint.Skip
    return am.dx_dt['x_6', t] == (1 / am.W) * ((-1 * (am.fa[t] + am.fb[t]) * am.x['x_6', t]) \
    + (k(am, 2, t) * am.x['x_2', t] * am.x['x_3', t] * am.W) \
    - 0.5 * (k(am, 3, t) * am.x['x_3', t] * am.x['x_6', t] * am.W))
mo.xp_rule = Constraint(mo.t, rule=_xp_rule)

# Boundary conditions 
for j in mo.x_is:
    mo.x[j, mo.t.first()].fix(mo.x_i0[j])
mo.Tr[mo.t.first()].fix(mo.Tr_0)
def _initfa(am, i):
    return am.fa[i]==am.fa_0
mo.fa_con = Constraint(mo.t, rule=_initfa)
mo.fb[mo.t.first()].fix(mo.fb_0)

# Run
dis = TransformationFactory('dae.finite_difference')
dis.apply_to(mo, nfe=1, wrt=mo.t, scheme='BACKWARD')
solver = SolverFactory('ipopt')
solver.options['halt_on_ampl_error'] = 'yes'
results = solver.solve(mo, tee=True, keepfiles=True, logfile = os.path.join(sv_dir, f'wo{t_is[-1]}.log'))

# Saving results 
res = pd.DataFrame(index=pd.Index([value(mo.t.last())], name='t'))
old = pd.read_csv(os.path.join(sv_dir, 'wo.csv')).set_index('t').dropna()
bounds = {}
res.index.name = 't'
for v in mo.component_objects(Var, active = True):
    subsets = list(v.index_set().subsets())
    if v.dim() == 1: 
        res[v.name] = [value(v[mo.t.last()])]
        bounds[v.name] = (v[mo.t.first()].lb, v[mo.t.first()].ub)
    elif v.dim() == 2: 
        J, _ = subsets
        for j_idx, j in enumerate(J, start=1):   
            res[f'{v.name}_{j_idx}'] = [value(v[j, mo.t.last()])]
            bounds[f'{v.name}_{j_idx}'] = (v[j, mo.t.first()].lb, v[j, mo.t.first()].ub) 
res = pd.concat([old, res]).to_csv(os.path.join(sv_dir, 'wo.csv'), index=True)
with open(os.path.join(sv_dir, f'wo{t_is[-1]}.txt') , 'w') as file:
    mo.pprint(ostream = file)
pd.DataFrame.from_dict(bounds, orient='index', columns=['lower', 'upper']).to_csv(os.path.join(sv_dir, 'bounds.csv'))

