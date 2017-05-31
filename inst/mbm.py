import warnings
import re
import numpy as np
import sys



def main():
    # gpy prints lots of warnings during optimization; normally it is safe to ignore these
    if '--warn' in sys.argv:
        warn = True
        warnings.simplefilter('default')
    else:
        warn = False
        warnings.filterwarnings("ignore")

    suffix = get_arg('out')
    gpyLoc = get_arg('gpy')
    if gpyLoc is not None:
        sys.path.append(gpyLoc)
    print sys.path
    import GPy

    # read x and y data files
    xFile = get_arg('x')
    xDat = read_mbm_data(xFile)
    yDat = read_mbm_data(get_arg('y'))

    model = MBM(xDat, yDat)
    fits = model.predict()
    np.savetxt(xFile + suffix, fits, delimiter=',')


def get_arg(arg):
    pat = re.compile('--' + arg + '=(.+)')
    return filter(None, [pat.match(x) for x in sys.argv])[0].group(1)


def read_mbm_data(fname):
    dat = np.genfromtxt(fname, delimiter=',', skip_header=1, names=False, dtype=float)
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
        self.model = GPy.core.GP(X=self.X, Y=self.Y)

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