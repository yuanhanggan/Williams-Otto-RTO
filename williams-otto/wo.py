import pyomo.environ as pyo 

#
# Model 
# 

model = pyo.AbstractModel()
model.n = pyo.Param(within=pyo.NonNegativeIntegers) # 1x1
model.m = pyo.Param(within=pyo.NonNegativeIntegers) # 1x1
model.I = pyo.RangeSet(1, model.m) # mx1 
model.J = pyo.RangeSet(1, model.n) # nx1
model.A = pyo.Param(model.I) # mx1 
model.Ea = pyo.Param(model.I) # mx1 
model.W = pyo.Param(within=pyo.NonNegativeReals) # 1x1
model.fa = pyo.Param(within=pyo.NonNegativeReals) # 1x1 
model.fb_cons = pyo.Param(pyo.RangeSet(1, 2), within=pyo.NonNegativeReals) # 2x1 
model.Tr_cons = pyo.Param(pyo.RangeSet(1, 2), within=pyo.NonNegativeReals) # 2x1

model.x = pyo.Var(model.J, within=pyo.UnitInterval, initialize=0.1) # mx1
model.fb = pyo.Var(initialize=6, within=pyo.NonNegativeReals) # 1x1 
model.Tr = pyo.Var(initialize=100, within=pyo.NonNegativeReals) # 1x1

# Rate equation
def k(am, I):
    return am.A[I] * pyo.exp(am.Ea[I] / am.Tr)

# Objective expression 
def obj_rule(am):
    return 5554.1 * (am.fa + am.fb) * am.x[6] \
    + 125.91 * (am.fa + am.fb) * am.x[4] \
    - 370.3 * am.fa \
    - 555.42 * am.fb

# Equality constraint expressions
def xa_bal(am):
    return 0 == am.fa - (am.fa + am.fb) * am.x[1] \
    - k(am, 1) * am.x[1] * am.x[2] * am.W
def xb_bal(am):
    return 0 == am.fb - (am.fa + am.fb) * am.x[2] \
    - (k(am, 1) * am.x[1] * am.x[2] * am.W) \
    - (k(am, 2) * am[2] * am[3] * am.W)
def xc_bal(am):
    return 0 == -1 * (am.fa + am.fb) * am.x[3] \
    + 2 * (k(am, 1) * am.x[1] * am.x[2] * am.W) \
    - 2 * (k(am, 2) * am[2] * am[3] * am.W) \
    - k(am, 3) * am[3] * am[6] * am.W
def xe_bal(am):
    return 0 == -1 * (am.fa + am.fb) * am.x[4] \
    + (k(am, 2) * am[2] * am[3] * am.W)
def xg_bal(am):
    return 0 == -1 * (am.fa + am.fb) * am.x[5] \
    + 1.5 * (k(am, 3) * am[3] * am[6] * am.W)
def xp_bal(am):
    return 0 == -1 * (am.fa + am.fb) * am.x[6] \
    + (k(am, 2) * am[2] * am[3] * am.W) \
    - 0.5 * (k(am, 3) * am[3] * am[6] * am.W)
def x_bal(am):
    return 1 == sum(am.x[i] for i in am.I)

# Inequality constraints
def fb_rule(am):
    return pyo.inequality(am.fb_cons[1], am.fb, am.fb_cons[2])
def tr_rule(am):
    return pyo.inequality(am.Tr_cons[1], am.Tr, am.Tr_cons[2])

# Instantiate 
model.obj = pyo.Objective(rule=obj_rule)
model.xa_bal = pyo.Constraint(rule=xa_bal)
model.xb_bal = pyo.Constraint(rule=xb_bal)
model.xc_bal = pyo.Constraint(rule=xc_bal)
model.xe_bal = pyo.Constraint(rule=xe_bal)
model.xg_bal = pyo.Constraint(rule=xg_bal)
model.xp_bal = pyo.Constraint(rule=xp_bal)
model.x_bal = pyo.Constraint(rule=x_bal)
model.fb_rule = pyo.Constraint(rule=fb_rule)
model.tr_rule = pyo.Constraint(rule=tr_rule)

# pyomo solve wo.py wo.dat --solver=glpk