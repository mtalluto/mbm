# set kernel constraints
import GPy
import numpy as np


def make_kernel(dim, sparse = False):
    """
    Make a kernel for mbm; presently supports no options
    """
    k = GPy.kern.RBF(input_dim=dim, ARD=True)
    if sparse:
        k.add(GPy.kern.White(dim))
    return k


def set_kernel_constraints_py(kern, lengthscale, which, pr = GPy.priors.Gamma.from_EV(1.,3.)):
    if which == 'all' or which == 'variance':
        kern.variance.set_prior(pr)
    if which == 'all' or which == 'lengthscale':
        kern.lengthscale.set_prior(pr)
    if lengthscale is not None:
        for i in range(len(lengthscale)):
            if not np.isnan(lengthscale[i]) and lengthscale[i] is not None:
                kern.lengthscale[i] = lengthscale[i]
                kern.lengthscale[[i]].fix()
    return(kern)
