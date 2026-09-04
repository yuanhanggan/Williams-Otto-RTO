# Note that for each example in any of these pyomo things, we need to follow these steps 
# Create model -> instantiate the model -> apply the solver -> interrogate the solver results 

import pyomo.environ as pyo
# Declare model 
model = pyo.AbstractModel()

# First we declare the number of decision variables (n) and the number of constraints (m) as parameters
# Note that these are non-negative integers, for very obvious reasons, they generalize to be indexes 
model.n = pyo.Param(within=pyo.NonNegativeIntegers)
model.m = pyo.Param(within=pyo.NonNegativeIntegers)
# Conveniently, the indexs which we would like to use are not dependent, and are known a priori 
# It becomes a good habit to declare them explicitly 
model.J = pyo.RangeSet(1, model.n)
model.I = pyo.RangeSet(1, model.m)
# Next, given that our indexes have been declare, we can move on to declare our indexed variables, which rely on these indexes
model.a = pyo.Param
