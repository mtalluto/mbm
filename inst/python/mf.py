import GPy
import numpy as np

def set_mean_function_py(dim, useMeanFunction = True):
    if useMeanFunction:
        mf = GPy.mappings.linear.Linear(dim, 1)
        mf = GPy.mappings.additive.Additive(GPy.mappings.constant.Constant(dim, 1), \
                GPy.mappings.linear.Linear(dim, 1))
        nm = mf.parameter_names()[1]
        mf[nm][0].constrain_positive()
        # fix other slopes to 0
        for i in range(1, dim):
            mf[nm][i].fix(0)
    else:
        mf = None
    return mf
