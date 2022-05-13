#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import poisson_conf_interval
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import constants as c
font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }
def get_R2(M2):
    print(M2)
    if M2<=0.069:
        return 0.118*(M2/0.069)**0.3
    elif 0.069<M2<0.2:
        return 0.225*(M2/0.20)**0.61
    elif 0.2<M2:
        return 0.293*(M2/0.20)**0.69

def get_radius(M2):
    if M2<1.66:a=0.026;b=0.945
    else:
        a=0.124;b=0.555
    R=np.exp(a+b*np.log(M2))
    return R

def read_fromfit():
    (c1, c2, c3, c4, c5) = (0.5126, 0.7388, 0.6710, 0.7349, 0.3983)
    path = '/Users/baotong/Desktop/period_gc/'
    sim_info = fits.open(path + 're_model.fit')
    M1 = sim_info[1].data.field(0)
    # 白矮星质量
    M2 = sim_info[1].data.field(1)
    # 伴星质量
    R2 = sim_info[1].data.field(2)
    # 伴星半径
    Period = sim_info[1].data.field(3)
    Sep = sim_info[1].data.field(4)
    q=M2/M1
    RL_on_a = c1 * q ** c2 / (c3 * q ** c4 + np.log(1 + q ** c5))
    RL=RL_on_a*Sep
    plt.plot(M2,R2)
    plt.plot(M2,RL)
    plt.show()

def plot_M2_R():
    (c1,c2,c3,c4,c5)=(0.5126,0.7388,0.6710,0.7349,0.3983)
    q=np.linspace(0.01,1.0,1000)
    period_list=np.linspace(6,16,5)
    color_list=['r','orange','y','green','blue','g','purple']
    M1=0.6*c.M_sun
    M2=q*M1
    RL_on_a=c1*q**c2/(c3*q**c4+np.log(1+q**c5))
    for i in range(len(period_list)):
        period=period_list[i]*3600*u.s
        # print((period**2*c.G*(M1+M2)/(4*np.pi**2))**(1/3.))
        # period=1*3600*u.s
        a=(period**2*c.G*(M1+M2)/(4*np.pi**2))**(1/3)*2
        R2=np.array([get_radius(M2[i].to(u.M_sun).value) for i in range(len(M2))])
        plt.figure(1,(9,6))
        plt.plot(M2.to(u.M_sun).value,RL_on_a*a.to(u.R_sun).value/R2,label=f'P={period.to(u.hour)}',color=color_list[i])
        plt.xlabel(r'$\rm M_2 (M_\odot)$', font1)
        plt.ylabel(r'$\rm R_{RL}/R_2$', font1)
        plt.plot([M2.to(u.M_sun).value.min(),1],[1,1],'--')
    plt.text(M2.to(u.M_sun).value.min(),1.1,r'$\rm R_{RL}=R_2$',font1)
    plt.tick_params(labelsize=18)
    plt.legend()
    plt.semilogy()
    plt.show()

if __name__ == '__main__':
    plot_M2_R()