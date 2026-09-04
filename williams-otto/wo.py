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
def k(amodel, I):
    return amodel.A[I] * pyo.exp(amodel.Ea[I] / amodel.Tr)

# Objective expression 
def obj_rule(amodel):
    return 5554.1 * (amodel.fa + amodel.fb) * amodel.x[6] \
    + 125.91 * (amodel.fa + amodel.fb) * amodel.x[4] \
    - 370.3 * amodel.fa \
    - 555.42 * amodel.fb

# Equality constraint expressions
def xa_bal(amodel):
    return amodel.fa - (amodel.fa + amodel.fb) * amodel.x[1] \
    - k(amodel, 1) * amodel.x[1] * amodel.x[2] * amodel.W