import numpy as np
from scipy.optimize import root_scalar
import os

mean_dx = 1
xmax = 200
def mean_fun(xmin):
    return ((1-a)/(2-a)) * (xmax**(2-a)-xmin**(2-a))/(xmax**(1-a)-xmin**(1-a))-mean_dx

colformat = '{:3.1f}'+'{:23.15e}'*3 + '\n'
if not os.path.exists('levy_parameters'):
    os.mkdir('levy_parameters')
for a in np.arange(1.2,4.01,0.1):
    file=open('levy_parameters/A'+str(int(round(10*a)))+'.dat','w')
    file.write('# mean=0.03, x_max='+str(xmax)+'\n')
    file.write('# alpha, A, B, EXP\n')
    xmin=root_scalar(mean_fun, bracket=[1.e-12,1]).root
    A   = xmax**(1-a)-xmin**(1-a)
    B   = xmin**(1-a)
    EXP = 1/(1-a)
    file.write(colformat.format(a,A,B,EXP))
