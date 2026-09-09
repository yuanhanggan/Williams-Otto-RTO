# Taken from Nicholson et al. 2018 
from pyomo.environ import * 
from pyomo.dae import * 


#
# Model 
# 

m = ConcreteModel()
m.pi = Param(initialize=3.1416)
m.t = ContinuousSet(bounds=(0, 2))
m.x = ContinuousSet(bounds=(0, 1))
m.u = Var(m.x, m.t)

# Declare derivatives in the model 
m.dudx = DerivativeVar(m.u, wrt=m.x)
m.dudx2 = DerivativeVar(m.u, wrt=(m.x, m.x))
m.dudt = DerivativeVar(m.u, wrt=m.t)

# Declare the PDE 
def _pde(m, i, j):
    if i == 0 or i == 1 or j == 0: 
        return Constraint.Skip
    return m.pi**2*m.dudt[i, j] == m.dudx2[i, j]

def _initcon(m, i): 
    if i == 0 or i == 1:
        return Constraint.Skip
    return m.u[i, 0] == sin(m.pi*i)
    
