# Optimizer for discrete time point
from pyomo.environ import *
import pandas as pd, os

# Model 
mo = ConcreteModel()
da = DataPortal()
da.load(filename=os.path.join(os.getcwd(), 'data', 'wo-rec.dat'))
sv_dir = os.path.join(os.getcwd(), 'sims', os.environ.get('SV_DIR'))
da_r = pd.read_csv(os.path.join(sv_dir, 'wo.csv')).set_index('t').dropna().to_dict(orient='index')
t_is = list(da_r.keys())
x_is = ['x_1','x_2','x_3','x_4','x_5','x_6']
dx_is = ['dx_dt_1','dx_dt_2','dx_dt_3','dx_dt_4','dx_dt_5','dx_dt_6']

# Parameters 
mo.x_is = Set(initialize=x_is) # nx1
mo.dx_is = Set(initialize=dx_is) # nx1
mo.A = Param(RangeSet(1, 3), initialize=da['A']) # mx1 
mo.Ea = Param(RangeSet(1, 3), initialize=da['Ea']) # mx1
mo.W = Param(initialize=da['W']) # 1x1  
mo.fa = Param(initialize=da_r[t_is[-1]]['fa'])
mo.x_i0 = Param(mo.x_is, initialize={k: da_r[t_is[-1]][k] for k in da_r[t_is[-1]] & mo.x_is}) # nx1
mo.dx_i0 = Param(mo.dx_is, initialize={k: da_r[t_is[-1]][k] for k in da_r[t_is[-1]] & mo.dx_is}) # nx1

# Variables 
def init_x(am, j):
    return am.x_i0[j]
def init_dx(am, j):
    return am.dx_i0[j]
mo.x = Var(mo.x_is, bounds=(0, 1), initialize=init_x) # nx1
# mo.dx = Var(mo.dx_is, within=Reals, initialize=init_dx)
mo.fb = Var(bounds=(2, 10), initialize=da_r[t_is[-1]]['fb']) # 1x1 
mo.Tr = Var(bounds=(323.15, 423.15), initialize=da_r[t_is[-1]]['Tr']) # 1x1

# Rate equation
def k(am, I):
    return am.A[I] * exp(am.Ea[I] / am.Tr)

# Objective expression 
def obj_rule(am):
    return -1 * ((5554.1 * (am.fa + am.fb) * am.x['x_6']) \
    + (125.91 * (am.fa + am.fb) * am.x['x_4']) \
    - (370.3 * am.fa) \
    - (555.42 * am.fb))
mo.obj = Objective(rule=obj_rule, sense=minimize)

# Equality constraint expressions
def xa_bal(am):
    return am.dx_i0['dx_dt_1'] == am.fa - ((am.fa + am.fb) * am.x['x_1']) \
    - (k(am, 1) * am.x['x_1'] * am.x['x_2'] * am.W)
mo.xa_bal = Constraint(rule=xa_bal)
def xb_bal(am):
    return am.dx_i0['dx_dt_2'] == am.fb - ((am.fa + am.fb) * am.x['x_2']) \
    - (k(am, 1) * am.x['x_1'] * am.x['x_2'] * am.W) \
    - (k(am, 2) * am.x['x_2'] * am.x['x_3'] * am.W)
mo.xb_bal = Constraint(rule=xb_bal)
def xc_bal(am):
    return am.dx_i0['dx_dt_3'] == (-1 * (am.fa + am.fb) * am.x['x_3']) \
    + 2 * (k(am, 1) * am.x['x_1'] * am.x['x_2'] * am.W) \
    - 2 * (k(am, 2) * am.x['x_2'] * am.x['x_3'] * am.W) \
    - (k(am, 3) * am.x['x_3'] * am.x['x_6'] * am.W)
mo.xc_bal = Constraint(rule=xc_bal)
def xe_bal(am):
    return am.dx_i0['dx_dt_4'] == (-1 * (am.fa + am.fb) * am.x['x_4']) \
    + 2 * (k(am, 2) * am.x['x_2'] * am.x['x_3'] * am.W)
mo.xe_bal = Constraint(rule=xe_bal)
def xg_bal(am):
    return am.dx_i0['dx_dt_5'] == (-1 * (am.fa + am.fb) * am.x['x_5']) \
    + 1.5 * (k(am, 3) * am.x['x_3'] * am.x['x_6'] * am.W)
mo.xg_bal = Constraint(rule=xg_bal)
def xp_bal(am):
    return am.dx_i0['dx_dt_6']== (-1 * (am.fa + am.fb) * am.x['x_6']) \
    + (k(am, 2) * am.x['x_2'] * am.x['x_3'] * am.W) \
    - 0.5 * (k(am, 3) * am.x['x_3'] * am.x['x_6'] * am.W)
mo.xp_bal = Constraint(rule=xp_bal)

# Solve 
solver = SolverFactory('ipopt')
solver.options['halt_on_ampl_error'] = 'yes'
results = solver.solve(mo, tee=True, keepfiles=True, logfile = os.path.join(sv_dir, 'logs', f'rto{t_is[-1]}.log'))
with open(os.path.join(sv_dir, 'logs', f'rto{t_is[-1]}.txt') , 'w') as file:
    mo.pprint(ostream = file)