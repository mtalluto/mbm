import warnings
import re
import numpy as np
import sys

## setup to get GPy up and running
def get_arg(arg):
    pat = re.compile('--' + arg + '=(.+)')
    result = filter(None, [pat.match(x) for x in sys.argv])
    if len(result) == 0:
        return None
    else:
        return [r.group(1) for r in result]

# gpy prints lots of warnings during optimization; normally it is safe to ignore these
if '--warn' in sys.argv:
    warn = True
    warnings.simplefilter('default')
else:
    warn = False
    warnings.filterwarnings("ignore")

gpyLoc = get_arg('gpy')
if gpyLoc is not None:
    sys.path.append(gpyLoc)
import GPy

def main():

    suffix = get_arg('out')[0]
    link = get_link(get_arg('link')[0])

    # read x and y data files, prediction files
    xFile = get_arg('x')[0]
    xDat = read_mbm_data(xFile)
    yDat = read_mbm_data(get_arg('y')[0])
    prFiles = get_arg('pr')
    if prFiles is not None:
        prDat = [read_mbm_data(prf) for prf in prFiles]

    model = MBM(xDat, yDat, link = link)
    fits = model.predict()
    np.savetxt(xFile + suffix, fits, delimiter=',')
    if prFiles is not None:
        for prd, prf in zip(prDat, prFiles):
            prFit = model.predict(prd)
            np.savetxt(prf + suffix, prFit, delimiter=',')


def read_mbm_data(fname):
    dat = np.genfromtxt(fname, delimiter=',', skip_header=1, names=None, dtype=float)
    if len(np.shape(dat)) == 1:
        dat = np.expand_dims(dat, 1)
    return dat


def get_link(linkname):
    if linkname == 'probit':
        return GPy.likelihoods.link_functions.Probit() 
    elif linkname == 'log':
        return GPy.likelihoods.link_functions.Log()
    else:
        return GPy.likelihoods.link_functions.Identity() 


class MBM(object):
    """
    Create an MBM model object

    x: numpy array containing covariates for the model; we assume the first column is distances and others are midpoints
    y: single-column 2D numpy array containing response data; should be the same number of rows as x

    value: Object of class MBM
    """
    def __init__(self, x, y, link):
        self.X = x
        self.Y = y
        self.kernel = GPy.kern.RBF(input_dim=np.shape(self.X)[1], ARD=True)
        self.link = link
        self.likelihood = GPy.likelihoods.Gaussian(gp_link = self.link)
        if isinstance(self.likelihood, GPy.likelihoods.Gaussian) and isinstance(self.link, GPy.likelihoods.link_functions.Identity):
            self.inference = GPy.inference.latent_function_inference.ExactGaussianInference()
        else:
            self.inference = GPy.inference.latent_function_inference.Laplace()
        self.model = GPy.core.GP(X=self.X, Y=self.Y, kernel = self.kernel, likelihood = self.likelihood, inference_method = self.inference)
        self.model.optimize()

    def predict(self, newX = None):
        """
        Predict an mbm model

        newX: new dataset, with same number of columns as the original X data; if None, predicts to input data

        value: numpy array of predictions, with same number of rows as newX
        """
        if newX == None:
            newX = self.X
        elif len(np.shape(newX)) == 1:
            newX = np.expand_dims(newX, 1)
        mean, variance = self.model.predict_noiseless(newX)
        return mean

main()