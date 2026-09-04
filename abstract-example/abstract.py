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
model.a = pyo.Param(model.I, model.J) # to give a mxn matrix 
model.b = pyo.Param(model.I) # to give a mx1 vector 
model.c = pyo.Param(model.J) # to give a nx1 vector 

# Next, we are ready to declare our decision variable, note that similar to our parameters, we want to index our decision variable, and that the constraint x_j >= 0 can be conveniently defined with as a parameter 
model.x = pyo.Var(model.J, domain=pyo.NonNegativeReals)

# Implementation of objective and constraint expressions are given as Python function statements 
# Note that the model is always passed as the first argument 
# Note that pyomo has its own inner product function
# Can I just use numpy? here
def obj_expression(m):
    return pyo.summation(m.c, m.x)

# Next, we tell the model that this is specifically the objective expression, which is later minimized 
model.OBJ = pyo.Objective(rule=obj_expression)

# Next, we declare the constraint expression
# Note that constraints are declared by the index i, and to parameterize the expression by i we can include it as a formal parameter to the function
def ax_constraint_rule(m, i): 
    # return the expression for the constraint for a particular i
    return sum(m.a[i, j]* m.x[j] for j in m.J) >= m.b[i]

# Note that we need to declare that we have i number of constraints 
model.AxbConstraint = pyo.Constraint(model.I, rule=ax_constraint_rule)


