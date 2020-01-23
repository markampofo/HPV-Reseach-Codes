# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 09:43:02 2017

@author: Mark Ampofo
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 09:29:49 2017

@author: Mark Ampofo
"""

import scipy.integrate as spi
import numpy as np
import pylab as pl

r     = 0.02#rate at which the virus enters the uninfected basal layer
delta = 0.0082 # rate of free variance in the basal layer
a     = 0.01 # rate cell differentiation from the uninfected  basal layer to the uninfected epidermal layer 
n     = 10000 # rate cell differentiation from the uninfected  epidermal layer to the uninfected keratinous layer 
alpha = 0.0001 # rate cell de-differentiation from the uninfected  epidermal layer to the uninfected basal layer 
c     = 50  # rate cell differentiation from the infected  basal layer to the infected epidermal layer 
p     = 13.44 # rate cell differentiation from the infected  epidermal layer to the infected basal layer 
b     = 1.00 # rate cell de-differentiation from the infected  epidermal layer to the infected basal layer 
theta = 2.03 # rate cell de-differentiation from the infected  keratinous layer  to the infected epidermal layer 
k     = 1.01  # rate cell differentiation from the free variance basal layer to the free variance epidermal layer  
u_1   = 0.2
TS    = 1.0
ND    = 10000
S0    = 1
I0    = 0
V0    = 1
P0    = 0
C0    = 0
INPUT = (S0, I0, V0, P0, C0)


def diff_eqs(INP,t):  
        '''The main set of equations'''
        Y1 = np.zeros((5))
        F = INP    
        Y1[0] = r * F[0] * (1 - (F[0] + F[1])) - alpha * F[0] * F[2] * (1 - u_1)
        Y1[1] = alpha * F[0] * F[2] * (1 - u_1) - a * F[1] - delta * F[1]
        Y1[2] = n * F[1] - c * F[2]
        Y1[3] = delta * p * F[1] + b * F[3] - theta * ((F[3]**2)/(1 + F[3]**2))
        Y1[4] = theta * ((F[3]**2)/(1 + F[3]**2)) - k * F[4]
        return Y1   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc)
RES = spi.odeint(diff_eqs,INPUT,t_range)

print (RES)









pl.subplot(111)
pl.plot(RES[:,3], '-g', label='Pre-cancer Cells with vaccination of 0.2')
#pl.plot(RES[:,1], '-y', label='Infected Cells with vaccination of 0.2')


u_1   = 0.4

def diff_eqs(INP,t):  
        '''The main set of equations'''
        Y = np.zeros((5))
        F = INP    
        Y[0] = r * F[0] * (1 - (F[0] + F[1])) - alpha * F[0] * F[2] * (1 - u_1)
        Y[1] = alpha * F[0] * F[2] * (1 - u_1) - a * F[1] - delta * F[1]
        Y[2] = n * F[1] - c * F[2]
        Y[3] = delta * p * F[1] + b * F[3] - theta * ((F[3]**2)/(1 + F[3]**2))
        Y[4] = theta * ((F[3]**2)/(1 + F[3]**2)) - k * F[4]
        return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc)
RES = spi.odeint(diff_eqs,INPUT,t_range)

pl.plot(RES[:,3], '-r', label='Pre-cancer Cells with vaccination of 0.4')


u_1   = 0.6

def diff_eqs(INP,t):  
        '''The main set of equations'''
        Y = np.zeros((5))
        F = INP    
        Y[0] = r * F[0] * (1 - (F[0] + F[1])) - alpha * F[0] * F[2] * (1 - u_1)
        Y[1] = alpha * F[0] * F[2] * (1 - u_1) - a * F[1] - delta * F[1]
        Y[2] = n * F[1] - c * F[2]
        Y[3] = delta * p * F[1] + b * F[3] - theta * ((F[3]**2)/(1 + F[3]**2))
        Y[4] = theta * ((F[3]**2)/(1 + F[3]**2)) - k * F[4]
        return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc)
RES = spi.odeint(diff_eqs,INPUT,t_range)

pl.plot(RES[:,3], '-b', label='Pre-cancer Cells with vaccination of 0.6')


u_1   = 0.8

def diff_eqs(INP,t):  
        '''The main set of equations'''
        Y = np.zeros((5))
        F = INP    
        Y[0] = r * F[0] * (1 - (F[0] + F[1])) - alpha * F[0] * F[2] * (1 - u_1)
        Y[1] = alpha * F[0] * F[2] * (1 - u_1) - a * F[1] - delta * F[1]
        Y[2] = n * F[1] - c * F[2]
        Y[3] = delta * p * F[1] + b * F[3] - theta * ((F[3]**2)/(1 + F[3]**2))
        Y[4] = theta * ((F[3]**2)/(1 + F[3]**2)) - k * F[4]
        return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc)
RES = spi.odeint(diff_eqs,INPUT,t_range)

pl.plot(RES[:,3], '-k', label='Pre-cancer Cells with vaccination of 0.8')






u_1   = 1

def diff_eqs(INP,t):  
        '''The main set of equations'''
        Y = np.zeros((5))
        F = INP    
        Y[0] = r * F[0] * (1 - (F[0] + F[1])) - alpha * F[0] * F[2] * (1 - u_1)
        Y[1] = alpha * F[0] * F[2] * (1 - u_1) - a * F[1] - delta * F[1]
        Y[2] = n * F[1] - c * F[2]
        Y[3] = delta * p * F[1] + b * F[3] - theta * ((F[3]**2)/(1 + F[3]**2))
        Y[4] = theta * ((F[3]**2)/(1 + F[3]**2)) - k * F[4]
        return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc)
RES = spi.odeint(diff_eqs,INPUT,t_range)

pl.plot(RES[:,3], '-y', label='Pre-cancer Cells with vaccination of 1')
pl.legend(loc=0)
#pl.title('Program_2_1.py')
pl.xlabel('Time (days)')
pl.ylabel('Pre-cancer Cells (P)')
pl.show()