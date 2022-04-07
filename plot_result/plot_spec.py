#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
from matplotlib.ticker import FuncFormatter
import pandas as pd

f1 = plt.figure(figsize=(15, 24))

plt.rc('text', usetex=True)
plt.rc('font', family='Helvetica', size=25)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

ax1 = plt.subplot(3, 1, 1, sharex=ax3, position=[0.2, 0.5, 0.75, 0.45])
# ax1 = AxesGrid(f1,311,share_all=True,)

# plt.title('Delta Chi-2 test of obs over model')

plt.ylabel(r'\textbf{normalization}($\rm ph\ s^{-1}\ cm^{-2}$) ')
# plt.semilogx()
# plt.semilogy()
plt.contour(X_E, Y_N, delta_chi, [-9.21, -4.61, -2.3, 0.50], colors=('blue', 'green', 'red', 'black'))
# CS=plt.contour(X_E,Y_N,delta_chi,[-1,0.10],colors = ('red','black'))
# plt.contour(X_E,Y_N,delta_chi2)
# plt.clabel(CS,fontsize=9)
a = []
plt.plot(a, c='blue', label=r'$99\%$')
# plt.plot(a,c = 'blue',label = r'$\Delta C$ = -9.21(99\%)')
plt.plot(a, c='green', label=r'$90\%$')
# plt.plot(a,c = 'green',label = r'$\Delta C$ = -4.61(90\%)')
plt.plot(a, c='red', label=r'$68\%$')
# plt.plot(a,c = 'red',label = r'$\Delta C$ = -2.3(68\%)')
plt.plot(a, c='black', label=r'$\Delta C$ = +0.5')


# plt.yticks(-0.2e-5,1.5e-5)
def formatnum(x, pos):
    return '$%.1f$x$10^{-5}$' % (x * 100000)


formatter = FuncFormatter(formatnum)
ax1.yaxis.set_major_formatter(formatter)
ax1.tick_params(width=1.4, length=5, labelbottom=False)
# ax1.set_xticks([5.8,6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4],[])

## unrecognized lines
# plt.axvline(1.617,c='y')
# plt.axvline(1.662,c='y')
# plt.axvline(1.938,c='y')
# plt.axvline(2.045,c='y')
# plt.axvline(2.483,c='y')
# plt.axvline(3.689,c='y')
# plt.axvline(5.06,c='y')
# plt.axvline(5.65,c='y')
plt.plot([6.39505864598, 6.68829096202, 6.88323121211, 7.04956380678], [5e-06, 3.5e-06, 5e-06, 3e-06], marker='+',
         linewidth=0, markersize=16, markeredgecolor='black', markeredgewidth=1.5)

# boerder line
# plt.axvline(6.6)
# plt.axvline(6.68)

plt.legend(loc=2, fontsize=22)
plt.xlim(5.7, 7.5)

# ax2=f1[1]
ax2 = plt.subplot(3, 1, 2, sharex=ax3, position=[0.2, 0.2, 0.75, 0.3])
# ax1 = AxesGrid(f1,311,share_all=True)
# ax2.plot(energy,ydata,drawstyle='steps-mid',lw = 1.0, c = 'k')
plt.errorbar(energy, ydata, xerr=width, yerr=y_err, lw=0.0, elinewidth=1.2, c='k')
plt.plot(energy, mod_tot, drawstyle='steps-mid', c='blue', lw=1.5, label='Total model')
plt.plot(energy, mod_pl, c='blue', linestyle='--', label='Powerlaw component')
plt.plot(energy, mod_1, c='r', linestyle=':', lw=1.5, label='Gaussian lines')
plt.plot(energy, mod_2, c='r', linestyle=':', lw=1.5)
plt.plot(energy, mod_3, c='r', linestyle=':', lw=1.5)
plt.plot(energy, mod_4, c='r', linestyle=':', lw=1.5)
plt.yticks((0, 0.00005, 0.0001, 0.00015, 0.0002, 0.00025))


def formatnum2(x, pos):
    return '$%.1f$x$10^{-4}$' % (x * 10000)


formatter2 = FuncFormatter(formatnum2)
ax2.yaxis.set_major_formatter(formatter2)
ax2.tick_params(width=1.4, length=5, labelbottom=False)
# ax2.set_xticks([5.8,6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4],[])
plt.xlim(5.7, 7.5)
plt.ylim(0, 3e-4)

plt.ylabel(r'\textbf{normalized counts} $\rm s^{-1}\ keV^{-1}\ cm^{-2}$')
plt.legend(loc=2)

# ax3=f1[2]
ax3 = plt.subplot(3, 1, 3, position=[0.2, 0.05, 0.75, 0.15])
plt.errorbar(energy, y_del, xerr=width, yerr=y_del_err, lw=0.0, elinewidth=1.2, c='k')
# plt.plot(energy,y_del,drawstyle='steps-mid',c='k')
plt.axhline(0, c='g')
plt.xlim(5.7, 7.5)
plt.ylim(-2, 2)
plt.yticks((-1, 1))
# ax3.set_xticks([5.8,6.2,6.6,7.0,7.4],[5.8,6.2,6.6,7.0,7.4])
plt.xlabel(r'\textbf{Energy} (keV)')
plt.ylabel(r'\textbf{(data-model)/error}')
ax3.tick_params(width=1.4, length=5)
ax3.set_xticks([5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4], [5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4])

f1.savefig('Delta_c_200324_sigma10_nobin2.eps', format='eps', papertype='ledger')
plt.show(f1)
# plt.ylim(-1e-5,2e-5)
# plt.savefig('Delta_c_200324_sigma10_nobin2.eps',format='eps',papertype='ledger')

# plt.savefig('blind_line_c_stat.png', dpi = 300)