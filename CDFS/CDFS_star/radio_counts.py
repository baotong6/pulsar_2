#!/bin/bash
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

def funC_N_dayu_S(S,deg2=0):
    ## S unit must be jy
    ## return num deg^2
    if deg2==0:
        return 60 * S ** (-1.5) * (np.pi / 180) ** 2
    else:
        return 60 * S ** (-1.5) * (np.pi / 180) ** 2*deg2

def test_Huynh_2012():
    dS_edge=np.array([50,79,126,199,315,500,792,1256,1991,3155,5000])
    S_bar=np.array([58,102,164,268,429,658,1069,1671,2196,3755])
    N=np.array([16,30,20,9,13,6,4,3,4])
    Nc=np.array([77.3,51.0,27.4,11.4,11.1,15.8,7.15,4.67,3.46,4.54])
    Nc_on_Nexp=np.array([0.47,0.81,0.90,0.81,1.62,4.24,4.06,5.11,4.73,14.97])*1e-2

    Nexp=Nc/Nc_on_Nexp
    N_test=funC_N_dayu_S(dS_edge[:-1]*1e-6,deg2=0.25)-funC_N_dayu_S(dS_edge[1:]*1e-6,deg2=0.25)
    print(Nexp)
    print(N_test)
    print(funC_N_dayu_S(S_bar*1e-6))

    print(np.sum(N)*3.14/18/0.25)
    print(np.sum(Nc)*3.14/18/0.25)


def test_Huynh_2015():
    def log_NcNexp_logS(N,S):
        ## S unit uJy
        if S<400:
            alpha=0.32
        elif S>=400:
            alpha=0.51
        return alpha*np.log10(S)**(alpha-1)*N

    S=np.arange(10,1000,900)
    Nc_on_Nexp=np.zeros(len(S))
    for i in range(len(S)):
        Nc_on_Nexp[i]=10**(log_NcNexp_logS(-1.65,S[i]))
    Nexp=funC_N_dayu_S(S)
    Nc=Nc_on_Nexp/Nexp
    print(np.sum(Nc))
test_Huynh_2015()
a1=79.3;a2=148.248
b1=89.76;b2=140.42823529411757