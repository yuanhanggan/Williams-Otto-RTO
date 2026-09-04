import pyomo.environ as pyo 

#
# Model 
# 

model = pyo.AbstractModel()

model.m = pyo.Param(within=pyo.NonNegativeIntegers) # 1x1
model.n = pyo.Param(within=pyo.NonNegativeIntegers) # 1x1
model.J = pyo.RangeSet(1, model.n) # nx1 
model.I = pyo.RangeSet(1, model.m) # mx1 
model.A = pyo.Param(model.I) # mx1 
model.Ea = pyo.Param(model.I) # mx1 
model.W = pyo.Param(within=pyo.NonNegativeReals) # 1x1
model.x = pyo.Param(model.I, within=pyo.UnitInterval) # nx1
model.fa = pyo.Param(within=pyo.NonNegativeReals) # 1x1 

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