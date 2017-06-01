import warnings
import re
import numpy as np
import sys

## setup to get GPy up and running
def get_arg(arg):
    pat = re.compile('--' + arg + '=(.+)')
    result = filter(None, [pat.match(x) for x in sys.argv])
    if len(result) > 0:
        return result[0].group(1)
    else:
        return None

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

    suffix = get_arg('out')

    # read x and y data files
    xFile = get_arg('x')
    xDat = read_mbm_data(xFile)
    yDat = read_mbm_data(get_arg('y'))

    model = MBM(xDat, yDat)
    fits = model.predict()
    np.savetxt(xFile + suffix, fits, delimiter=',')




def read_mbm_data(fname):
    dat = np.genfromtxt(fname, delimiter=',', skip_header=1, names=None, dtype=float)
    if len(np.shape(dat)) == 1:
        dat = np.expand_dims(dat, 1)
    return dat


class MBM(object):
    """
    Create an MBM model object

    x: numpy array containing covariates for the model; we assume the first column is distances and others are midpoints
    y: single-column 2D numpy array containing response data; should be the same number of rows as x

    value: Object of class MBM
    """
    def __init__(self, x, y):
        self.X = x
        self.Y = y
        self.kernel = GPy.kern.RBF(input_dim=np.shape(self.X)[1], ARD=True)
        self.link = GPy.likelihoods.link_functions.Identity() 
        self.likelihood = GPy.likelihoods.Gaussian(gp_link = self.link)
        self.inference = GPy.inference.latent_function_inference.ExactGaussianInference()
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