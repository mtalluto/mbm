### https://climin.readthedocs.io/en/latest/tutorial.html#

def run_svgp(X, Y, Z, kernel, likelihood, mf, bs, init, maxiter, verbose):
    import GPy
    import climin
    import sys
    # import copy

    # some fixed settings
    miniter = 100               # minimum number of iterations before checking for stopping
    lhood_threshold = 0.001     # the smallest proportion change that will trigger stopping
    lhood_run_threshold = 100   # how long should we be below threshold before stopping


    model = GPy.core.SVGP(X=X, Y=Y, Z = Z, kernel = kernel, likelihood = likelihood, \
            mean_function = mf, batchsize = bs, initialize = init)
    if init:
        model.randomize()
        model.Z.unconstrain()
    opt = climin.Adadelta(model.optimizer_array, model.stochastic_grad, step_rate=0.2, \
            momentum=0.9)
    # pars = copy.deepcopy(model.optimizer_array)
    lhood = []
    lhood_run = 0
    for info in opt:
        lhood.append(model.objective_function())
        if info['n_iter'] > 1:
            lhood_change = (abs(lhood[-1] - lhood[-2]) / ((lhood[-1] + lhood[-2])/2))
            if lhood_change <= lhood_threshold:
                lhood_run += 1
            else:
                lhood_run = 0

        if verbose and info['n_iter'] % 500 == 0:
            print("iteration " + str(info['n_iter']) + "   Objective: " + str(lhood[-1]) + \
                    "   delta: " + str(lhood_change))
            sys.stdout.flush()

        # stopping criteria
        if info['n_iter'] >= maxiter or lhood_run >= lhood_run_threshold:
            break
    return [model, len(lhood)]

