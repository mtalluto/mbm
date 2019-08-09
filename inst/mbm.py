

def main():
    # are we resuming a model to predict after the fact?
    pars = None
    # if '--resume' in sys.argv:
        # pars = read_mbm_data(parFile, reshape = False)

    # model = MBM(xDat, yDat, link = link, samples = n_samples, lengthscale = ls, \
    #         mean_function = mean_func, params=pars, svgp = svgp, batchsize = batchsize, \
    #         zsize = zsize, svgp_maxiter = svgp_iter)
    # fits = model.predict()




class MBM(object):
    def __init__(self, x, y, params = None):
    # def __init__(self, x, y, link, samples, lengthscale, mean_function, params = None, \
    #             svgp = False, z = None, batchsize = 20, zsize = 10, svgp_maxiter = 10000):
        # set up the model; either we do it from scratch or we re-initialize if we were 
        # passed a parameter array
        if svgp:
        # if params is None:
            self.optimize()
        # else:  #### this is code to use if RE-INITIALIZING model (ie., params was not None)
            self.model.update_model(False)
            self.model.initialize_parameter()
            self.model[:] = params
            self.model.update_model(True)    


    def init_svgp(self, batchsize, z, params, zsize, svgp_maxiter):
        # if z is None:
        #     # self.Z = np.zeros((zsize,np.shape(self.X)[1]))
        #     # generate inducing inputs along the range of x
        #     for xind in range(np.shape(self.X)[1]):
        #         mn = np.amin(self.X[:,xind])
        #         mx = np.amax(self.X[:,xind])
        #         self.Z[:,xind] = np.random.rand(zsize)*(mx-mn) + mn
        # else:
        #     self.Z = z
        # self.link = GPy.likelihoods.link_functions.Identity()
        # self.likelihood = GPy.likelihoods.Gaussian(gp_link = self.link)
        initialize = params is None
        # self.batchsize = batchsize
        # self.svgp_maxiter = svgp_maxiter
        # self.model = GPy.core.SVGP(X=self.X, Y=self.Y, Z = self.Z, kernel = self.kernel, \
        #         likelihood = self.likelihood,  mean_function = self.mean_function, \
        #         batchsize = self.batchsize, initialize = initialize)
        # if initialize:
        #     self.model.randomize()
        #     self.model.Z.unconstrain()


    def predict(self, newX = None):
        if self.samples is None:
            ##### moved to R
        else:
            preds = self.model.posterior_samples_f(newX, self.samples)
        return preds






    # def params(self):
    #     return self.model.param_array

    # def param_names(self):
        # res = []
        # inducing inputs
        # try:
        #     for i in range(np.shape(self.Z)[0]):
        #         for j in range(np.shape(self.Z)[1]):
        #             res.append("inducing_inputs." + str(i) + "." + str(j))
        # except AttributeError:
        #     pass
        # if self.mean_function is not None:
        #     res.append("mf.intercept")
        #     res = res + ["mf.slope." + str(i) for i in range(np.shape(self.X)[1])]
        # res.append('rbf.variance')
        # res = res + ["rbf.lengthscale." + str(i) for i in range(np.shape(self.X)[1])]
        # add white noise here
        # if isinstance(self.kernel, GPy.kern.src.add.Add):
        #     res.append('White_noise.variance')
        # res.append('Gaussian_noise.variance')
        # if isinstance(self.kernel, GPy.kern.src.add.Add):
        #     # nz = np.shape(self.Z)[0]
        #     # nch = (nz * (nz+1))/2
        #     res = res + ["u_cholesky." + str(i) for i in range(nch)]
        #     res = res + ["u_mean." + str(i) for i in range(nz)]
        # return res






