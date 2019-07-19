import warnings
import GPy
import numpy as np
# import re
# import sys
# import os





def main():

    # suffix = get_arg('out')[0]
    # link = get_link(get_arg('link')[0])
    # parFile = get_arg('par')[0]
    # n_samples = get_arg('sample')
    # if n_samples is not None:
    #     n_samples = int(n_samples[0])


    # look for fixed lengthscales
    ls = get_arg('ls')
    if ls is not None:
        ls = map(float, ls[0].split(','))

    # do we have a mean function?
    mean_func = '--mf' in sys.argv

    # ues the svgp model?
    if '--svgp' in sys.argv:
        svgp = True
        batchsize = int(get_arg('bs')[0])
        zsize = int(get_arg('inducing')[0])
        svgp_iter = int(get_arg('svgp_iter')[0])
    else:
        svgp = False
        batchsize = 10
        zsize = 10
        svgp_iter = 1000

    # are we resuming a model to predict after the fact?
    pars = None
    if '--resume' in sys.argv:
        pars = read_mbm_data(parFile, reshape = False)

    model = MBM(xDat, yDat, link = link, samples = n_samples, lengthscale = ls, \
            mean_function = mean_func, params=pars, svgp = svgp, batchsize = batchsize, \
            zsize = zsize, svgp_maxiter = svgp_iter)
    fits = model.predict()


# def get_link(linkname):
#     if linkname == 'probit':
#         return GPy.likelihoods.link_functions.Probit() 
#     elif linkname == 'log':
#         return GPy.likelihoods.link_functions.Log()
#     else:
#         return GPy.likelihoods.link_functions.Identity() 


class MBM(object):
    """
    Create an MBM model object

    x: numpy array containing covariates for the model; we assume the first column is 
        distances and others are midpoints y: single-column 2D numpy array containing 
        response data; should be the same number of rows as x
    link: A GPy link function object; see get_link()
    samples: the number of samples to take
    lengthscale: fixed lengthscales to use; if None, all will be optimized; if not None, 
        nan or None elements will be optimized
    mean_function: boolean; should we use a linear increasing mean function for the first 
        x-variable?
    params: an array of parameters; if none, a new model will be created and fit
    svgp: True/False, should the model be fit using stochastic variational (sparse) GP
    z: starting inducing inputs; if none, random numbers will be selected (of size 
        zsize)
    batchsize: batchsize to use with svgp
    zsize: number of inducing inputs (per X dimension)
    svgp_maxiter: maximum number of opimiser iterations for SVGP

    value: Object of class MBM
    """
    def __init__(self, x, y, params = None):
    # def __init__(self, x, y, link, samples, lengthscale, mean_function, params = None, \
    #             svgp = False, z = None, batchsize = 20, zsize = 10, svgp_maxiter = 10000):
        self.init_gp_params(x, y, samples, svgp, lengthscale, link, mean_function)
        # set up the model; either we do it from scratch or we re-initialize if we were 
        # passed a parameter array
        if svgp:
            self.init_svgp(batchsize, z, params, zsize, svgp_maxiter)
        else:
            self.init_gp(params)

        # if params is None:
            self.optimize()
        # else:
            self.model.update_model(False)
            self.model.initialize_parameter()
            self.model[:] = params
            self.model.update_model(True)    

    def init_gp(self, params):
        # initialize = params is None
        # self.model = GPy.core.GP(X=self.X, Y=self.Y, kernel = self.kernel, \
        #         likelihood = self.likelihood, inference_method = self.inference, \
        #         initialize = initialize)
    #     self.model = GPy.core.GP(X=self.X, Y=self.Y, kernel = self.kernel, \
    #             likelihood = self.likelihood, inference_method = self.inference, \
    #             mean_function = self.mean_function, initialize = initialize)

    def init_svgp(self, batchsize, z, params, zsize, svgp_maxiter):
        if z is None:
            self.Z = np.zeros((zsize,np.shape(self.X)[1]))
            # generate inducing inputs along the range of x
            for xind in range(np.shape(self.X)[1]):
                mn = np.amin(self.X[:,xind])
                mx = np.amax(self.X[:,xind])
                self.Z[:,xind] = np.random.rand(zsize)*(mx-mn) + mn
        else:
            self.Z = z
        self.link = GPy.likelihoods.link_functions.Identity()
        self.likelihood = GPy.likelihoods.Gaussian(gp_link = self.link)
        initialize = params is None
        self.batchsize = batchsize
        self.svgp_maxiter = svgp_maxiter
        self.model = GPy.core.SVGP(X=self.X, Y=self.Y, Z = self.Z, kernel = self.kernel, \
                likelihood = self.likelihood,  mean_function = self.mean_function, \
                batchsize = self.batchsize, initialize = initialize)
        if initialize:
            self.model.randomize()
            self.model.Z.unconstrain()

    def optimize(self):
        if isinstance(self.model, GPy.core.svgp.SVGP):
            # import climin
            # opt = climin.Adadelta(self.model.optimizer_array, self.model.stochastic_grad, \
                # step_rate=0.2, momentum=0.9)
            def max_iter(i):
                return i['n_iter'] > self.svgp_maxiter
            opt.minimize_until([max_iter])
        # else:
            # self.model.optimize()

    def predict(self, newX = None):
        if self.samples is None:
            ##### moved to R
        else:
            preds = self.model.posterior_samples_f(newX, self.samples)
        return preds






    # def params(self):
    #     return self.model.param_array

    def param_names(self):
        res = []
        # inducing inputs
        try:
            for i in range(np.shape(self.Z)[0]):
                for j in range(np.shape(self.Z)[1]):
                    res.append("inducing_inputs." + str(i) + "." + str(j))
        except AttributeError:
            pass
        # if self.mean_function is not None:
        #     res.append("mf.intercept")
        #     res = res + ["mf.slope." + str(i) for i in range(np.shape(self.X)[1])]
        # res.append('rbf.variance')
        # res = res + ["rbf.lengthscale." + str(i) for i in range(np.shape(self.X)[1])]
        # add white noise here
        if isinstance(self.kernel, GPy.kern.src.add.Add):
            res.append('White_noise.variance')
        # res.append('Gaussian_noise.variance')
        if isinstance(self.kernel, GPy.kern.src.add.Add):
            nz = np.shape(self.Z)[0]
            nch = (nz * (nz+1))/2
            res = res + ["u_cholesky." + str(i) for i in range(nch)]
            res = res + ["u_mean." + str(i) for i in range(nz)]
        return res






