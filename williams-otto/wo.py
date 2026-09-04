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

model.fb = pyo.Var(initialize=5)

