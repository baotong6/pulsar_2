import numpy as np
import matplotlib.pyplot as plt
import emcee

file='irs13E2.txt'
date=np.loadtxt(file)[:,0]
RA=np.loadtxt(file)[:,1]
DEC=np.loadtxt(file)[:,2]
RA_err=np.loadtxt(file)[:,3]
DEC_err=np.loadtxt(file)[:,4]
np.random.seed(45)
xdata=date
ydata=RA
yerr=RA_err

# theta_true = (25, 0.5)
# xdata = 100 * np.random.random(20)
# ydata = theta_true[0] + theta_true[1] * xdata
#
# # add scatter to points
# xdata = np.random.normal(xdata, 10)
# ydata = np.random.normal(ydata, 10)


# Create some convenience routines for plotting
# 让我们做一些辅助工作来可视化数据

def compute_sigma_level(trace1, trace2, nbins = 20):
    """From a set of traces, bin by number of standard deviations"""
    L, xbins, ybins = np.histogram2d(trace1, trace2, nbins)
    L[L == 0] = 1E-16
    logL = np.log(L)

    shape = L.shape
    L = L.ravel()

    # obtain the indices to sort and unsort the flattened array
    i_sort = np.argsort(L)[::-1]
    i_unsort = np.argsort(i_sort)

    L_cumsum = L[i_sort].cumsum()
    L_cumsum /= L_cumsum[-1]

    xbins = 0.5 * (xbins[1:] + xbins[:-1])
    ybins = 0.5 * (ybins[1:] + ybins[:-1])

    return xbins, ybins, L_cumsum[i_unsort].reshape(shape)


def plot_MCMC_trace(ax, xdata, ydata, trace, scatter = False, **kwargs):
    """Plot traces and contours"""
    xbins, ybins, sigma = compute_sigma_level(trace[0], trace[1])
    ax.contour(xbins, ybins, sigma.T, levels = [0.683, 0.955,0.9973], **kwargs)
    if scatter:
        ax.plot(trace[0], trace[1], ',k', alpha = 0.1)
    ax.set_xlabel('a0')
    ax.set_ylabel('a1')

def plot_MCMC_trace2(ax, xdata, ydata, trace, scatter = False, **kwargs):
    xbins2, ybins2, sigma2 = compute_sigma_level(trace[1], trace[2])
    ax.contour(xbins2, ybins2, sigma2.T, levels = [0.683, 0.955, 0.9973], **kwargs)
    if scatter:
        ax.plot(trace[1], trace[2], ',k', alpha = 0.1)
    ax.set_xlabel('a1')
    ax.set_ylabel('a2')



def plot_MCMC_model(ax, xdata, ydata, trace):
    """Plot the linear model and 2sigma contours"""
    ax.plot(xdata, ydata, 'ok')

    a0, a1, a2= trace[:3]
    xfit = xdata
    yfit = a0[:, None] + a1[:, None] * (xfit - 2012.0)+a2[:, None] * (xfit - 2012.0) ** 2
    mu = yfit.mean(0)
    sig = 2 * yfit.std(0)

    ax.plot(xfit, mu, '-k')
    ax.fill_between(xfit, mu - sig, mu + sig, color = 'lightgray')
    ax.plot(np.zeros(len(yfit))+2012.0,yfit,'--',color='yellow')
    ax.set_xlabel('date')
    ax.set_ylabel('DEC')


def plot_MCMC_results(xdata, ydata, trace, colors = 'k'):
    """Plot both the trace and the model together"""
    fig, ax = plt.subplots(1, 2, figsize = (10, 5))
    # plot_MCMC_trace(ax[0], xdata, ydata, trace, True, colors = colors)
    plot_MCMC_trace2(ax[0], xdata, ydata, trace, True, colors = colors)
    plot_MCMC_model(ax[1], xdata, ydata, trace)
    # plt.savefig('DEC_3b.eps')


# Define our posterior using Python functions
# for clarity, I've separated-out the prior and likelihood
# but this is not necessary. Note that emcee requires log-posterior

# 用Python函数定义后验，我把先验和似然估计分开写了，其实没必要，主要是显得更简洁。
# 注意emcee需要对数后验证

def log_prior(theta):
    a0, a1, a2= theta
    if a0 > 0:
        return 0.0  # log(0)
    else:
        #return -1.5 * np.log(1 + a1 ** 2) - np.log(sigma)
        return 0.99


def log_likelihood(theta, x, y,yerr):
    a0, a1, a2 = theta
    sigma = yerr
    y_model = a0 + a1 * (x-2012.0)+a2*(x-2012.0)**2
    # sigma = yerr ** 2 + y_model** 2 * np.exp(2 * log_f)
    # return -0.5 * np.sum((y - y_model) ** 2 / sigma2 + np.log(sigma2))

    return -0.5 * np.sum(np.log(2 * np.pi * sigma ** 2) + (y - y_model) ** 2 / sigma ** 2)


def log_posterior(theta, x, y,yerr):
    return log_prior(theta) + log_likelihood(theta, x, y,yerr)


# Here we'll set up the computation. emcee combines multiple "walkers",
# each of which is its own MCMC chain. The number of trace results will
# be nwalkers * nsteps

# 设置计算参数。emcee组合了多个"walkers"，每个都有自己的MCMC链。
# 跟踪结果的数量为 nwalkers * nsteps

ndim = 3 # number of parameters in the model
nwalkers = 100  # number of MCMC walkers
nburn = 200  # "burn-in" period to let chains stabilize
nsteps = 1000  # number of MCMC steps to take

# set theta near the maximum likelihood, with
np.random.seed(0)
starting_guesses = np.random.random((nwalkers, ndim))
starting_guesses[:,0]+=-1.7
# Here's the function call where all the work happens:
# we'll time it using IPython's %time magic


sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args = [xdata, ydata,yerr])
sampler.run_mcmc(starting_guesses, nsteps)
print("done")

# sampler.chain is of shape (nwalkers, nsteps, ndim)
# we'll throw-out the burn-in points and reshape:

# sampler.chain返回数组维度为(nwalkers, nsteps, ndim)

sampler.chain
emcee_trace = sampler.chain[:, nburn:, :].reshape(-1, ndim).T
# print(emcee_trace.shape)
print(emcee_trace)
plot_MCMC_results(xdata, ydata, emcee_trace)
plt.show()