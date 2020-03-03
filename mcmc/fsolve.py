import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root,fsolve
def f1(t):
    return np.e**(-t)*(t**2+t+1)-1




sol1_fsolve = fsolve(f1,[4])
print(sol1_fsolve)

vc=1./sol1_fsolve-(1./sol1_fsolve+1)*np.e**(-sol1_fsolve)
print(vc)