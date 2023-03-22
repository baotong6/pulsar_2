import numpy as np
import matplotlib.pyplot as plt
import hawkeye as hawk
import rocket as rocket
import rednoise as rednoise
from scipy import integrate
from scipy.optimize import curve_fit
import scipy.stats
from lmfit import Model
from astropy.modeling import models

def check_fD():
    N=100
    path='/Users/baotong/Desktop/period_Tuc/simulation/result_rd_sin/317_amp_0.4/'
    list_id=[];list_prob=[];list_period=[]
    for i in range(N):
        filename=f'result_10h_{i+1}.txt'
        res=np.loadtxt(path+filename)
        list_id.append(res[0])
        list_prob.append(res[2])
        list_period.append(res[3])
    list_prob=np.array(list_prob);list_period=np.array(list_period)
    # fD = len(np.where(list_prob > 0.99)[0]) / N
    real_period=23583.11
    detindex=np.where((list_prob>0.9)&(np.abs(list_period-real_period)<0.05*real_period))[0]
    period=list_period[detindex]
    detrate=len(detindex)/N
    print(detrate)
    plt.scatter(period,period)
    plt.show()
    # plt.hist(period,bins=20,histtype='step')
    # plt.show()
    std=np.std(period)
    print(std)

check_fD()