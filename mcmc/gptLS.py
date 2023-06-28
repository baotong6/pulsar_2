import numpy as np
import matplotlib.pyplot as plt

import numpy as np
from numba import jit

@jit(nopython=True)
def weighted_mean(y, w):
    weighted_sum = np.sum(w * y)
    return weighted_sum / np.sum(w)

@jit(nopython=True)
def off_diagonal_sum(t, omega):
    sin_2t = np.sin(2 * omega * t)
    cos_2t = np.cos(2 * omega * t)
    sin_t = np.sin(omega * t)
    cos_t = np.cos(omega * t)
    return np.sum(sin_2t), np.sum(cos_2t), np.sum(sin_t), np.sum(cos_t)

@jit(nopython=True)
def compute_amplitude(t, y, omega, weighted_y_sum, weighted_cos_sum,
                      weighted_sin_sum, num_data_points, num_max_harmonics):
    sin_2t, cos_2t, sin_t, cos_t = off_diagonal_sum(t, omega)
    cos_n2t = 2 * cos_2t**2 - 1
    sin_n2t = 2 * sin_2t * cos_2t
    cos_nt = np.empty(num_max_harmonics, dtype=np.float64)
    sin_nt = np.empty(num_max_harmonics, dtype=np.float64)
    cos_nt[0] = cos_t
    sin_nt[0] = sin_t
    for n in range(num_max_harmonics-1):
        cos_nt[n+1] = cos_nt[n] * cos_t - sin_nt[n] * sin_t
        sin_nt[n+1] = sin_nt[n] * cos_t + cos_nt[n] * sin_t
    c = ((weighted_y_sum * cos_t - weighted_sin_sum * sin_t)**2 /
         (weighted_cos_sum**2 -
          ((num_data_points / 2) * (cos_n2t - sin_n2t / num_data_points))))
    s = (weighted_y_sum**2 /
         (weighted_cos_sum**2 +
          ((num_data_points / 2) * (sin_n2t + cos_n2t / num_data_points))))
    a = np.sqrt((2 / num_data_points) * (c + s))
    return a

@jit(nopython=True)
def compute_binned_mags_and_counts(t, y, min_period, max_period,
                                   bin_width_factor=5):
    period_grid_size = int(np.ceil((1 / min_period - 1 / max_period) /
                                   (0.5 * bin_width_factor / len(t))))
    period_step_size = 1 / (max_period * 2 * bin_width_factor)
    period_grid = np.arange(1 / max_period + period_step_size / 2,
                            1 / min_period, period_step_size)
    min_per_arr = 1 / (period_grid + period_step_size / 2)
    max_per_arr = 1 / (period_grid - period_step_size / 2)
    y_sum_arr = np.zeros_like(period_grid)
    y_var_arr = np.zeros_like(period_grid)
    counts_arr = np.zeros_like(period_grid)
    for i in range(len(period_grid)):
        within_bin = (t >= t[0] + max_per_arr[i]) & (t <= t[-1] + min_per_arr[i])
        if np.sum(within_bin) == 0:
            continue
        t_bin = t[within_bin]
        y_bin = y[within_bin]
        bin_width = max_per_arr[i] - min_per_arr[i]
        num_bins = int(np.ceil(bin_width / np.median(np.diff(t_bin))))
        binned_t = np.linspace(t_bin[0], t_bin[-1] + bin_width, num_bins+1)
        binned_mean = np.zeros(num_bins)
        for j in range(num_bins):
            bin_mask = (t_bin >= binned_t[j]) & (t_bin < binned_t[j+1])
            binned_mean[j] = np.mean(y_bin[bin_mask])
        y_sum_arr[i] = np.sum(binned_mean)
        y_var_arr[i] = np.var(binned_mean) * num_bins
        counts_arr[i] = num_bins
    return y_sum_arr, y_var_arr, counts_arr, period_grid

def lomb_scargle(t, y, freq, normalization='standard', subtract_mean=True,
                 nyquist_factor=5, num_iterations=3, fit_mean=True,
                 center_data=True):
    """
    Compute the generalized Lomb-Scargle periodogram.

    Parameters
    ----------
    t : array_like
        The input times, assumed to be sorted in ascending order.
    y : array_like
        The input values.
    freq : array_like
        The frequencies to evaluate the periodogram.
    normalization : string, optional
        The normalization to apply to the periodogram. Options include
        "standard", "model" and "loglik". Default is "standard".
    subtract_mean : bool, optional
        If set to True, the mean of the data will be subtracted before computing
        the periodogram. Default is True.
    nyquist_factor : float, optional
        The maximum frequency to evaluate the periodogram will be the Nyquist
        frequency times this value. Default is 5.
    num_iterations : int, optional
        The number of iterations to use in the fitting procedure. Default is 3.
    fit_mean : bool, optional
        If set to True, a constant offset will be fit to the data. Default is
        True.
    center_data : bool, optional
        If set to True, the data will be centered before computing the
        periodogram. Default is True.

    Returns
    -------
    freq : ndarray
        The frequencies used in the periodogram, as input to the function.
    pgram : ndarray
        The periodogram.
    """
    if center_data:
        y = y - np.mean(y)

    if subtract_mean and fit_mean:
        y = y - weighted_mean(y, np.ones_like(y))

    if nyquist_factor is None:
        max_frequency = 0.5 / np.median(np.diff(t))
    else:
        max_frequency = nyquist_factor / np.median(np.diff(t))
    freq = np.asarray(freq)
    min_frequency = freq.min()

    # set up the frequency grid
    df = np.median(np.diff(freq))
    freq_grid = np.arange(min_frequency, max_frequency, df)

    # set up input variables
    num_data_points = len(t)
    weighted_y_sum = np.sum(y)
    weighted_cos_sum = np.sum(np.cos(2 * np.pi * freq[np.newaxis, :] *
                                      t[:, np.newaxis]) *
                              np.ones_like(freq[np.newaxis, :]), 0)
    weighted_sin_sum = np.sum(np.sin(2 * np.pi * freq[np.newaxis, :] *
                                      t[:, np.newaxis]) *
                              np.ones_like(freq[np.newaxis, :]), 0)
    sigma2 = np.sum((y - weighted_mean(y, np.ones_like(y)))**2)

    num_max_harmonics = 1 + (2 * num_data_points / 3)
    num_max_harmonics = np.minimum(num_max_harmonics, len(t))

    all_amps = np.empty((num_iterations, len(freq_grid)))

    for i in range(num_iterations):
        # compute the amplitude at each frequency on the grid
        for j in range(len(freq_grid)):
            omega = 2 * np.pi * freq_grid[j]
            a = compute_amplitude(t, y, omega, weighted_y_sum,
                                   weighted_cos_sum, weighted_sin_sum,
                                   num_data_points, num_max_harmonics)
            all_amps[i, j] = a

        # subtract the highest amplitude model from the data and repeat
        if i < num_iterations - 1:
            next_y = y - all_amps[i,:] @ np.cos(2 * np.pi * freq_grid[np.newaxis, :] *
                                                t[:, np.newaxis]).T
            weighted_y_sum = np.sum(next_y)
            y = next_y

    # compute the periodogram
    pgram = 0.5 * all_amps**2 / sigma2

    # perform the normalization
    if normalization == 'standard':
        norm = np.sum(pgram)
        pgram /= norm
    elif normalization == 'model':
        norm = np.max(pgram)
        pgram /= norm
    elif normalization == 'loglik':
        log_norm = np.sum(np.log(pgram))
        pgram *= np.exp(-log_norm)

    return freq_grid, pgram

if __name__=='__main__':
    # Generate test data
    np.random.seed(1234)
    t = np.linspace(0, 10, 1000)
    y = 0.5 * np.sin(2*np.pi*2*t) + np.random.randn(len(t))+10

    # Compute the Lomb-Scargle periodogram
    freq = np.linspace(0.01, 10, 1000)
    p = lomb_scargle(t, y, freq)
    plt.plot(t,y)
    # Plot the results
    plt.plot(freq, p)
    plt.xlabel('Frequency')
    plt.ylabel('Power')
    plt.show()