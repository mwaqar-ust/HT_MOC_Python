# Packages
import math
import numpy as np
import matplotlib.pyplot as plt


# %% constant variables
nu = 0.4
fDW = 0.02
g = 9.81
K = 2.1*10**9
rho = 1000

# %% Pipe Properties
L = 2000
Di = 0.0792
Do = 0.09
e = (Do - Di)/2
A = math.pi*0.25*Di**2
E = 120*10**9
a = 1000; #(rho*(1/K+(1-nu**2)*2/E*(Di/2)/e))**(-0.5)

# %% Computation Constant Parameters
Ca = g*A/a
Z = 1/Ca
R = fDW/(2*Di*A)
mio=fDW/(2*g*Di*A**2)


# %% Steady-state pressure head and flow-rate at the downstream end
Qv = 2/1000
Hv = 40
V0 = Qv/A

# %% Data Sampling & Discretization
Fs = 10
dt = 1/Fs
period = 4*L/a
T = 4*period
t = np.arange(0,T+dt,dt)
nt = len(t)

dx = a*dt;
x = np.arange(0,L+dx,dx)
nx = len(x)

t_lim = np.arange(1,nt,1)
x_lim = np.arange(0,nx,1)

# %% Define Steady State Values
Hi = np.zeros((nx,1))
Qi = np.zeros((nx,1))

Hi[-1] = Hv
Qi[-1] = Qv


for i in np.flip(x_lim[0:-1],0):
    Hi[i,0] = Hi[i+1,0]+mio*dx*Qi[i+1,0]**2
    Qi[i,0] = Qi[i+1,0]
    

# %% Define Matrices
h = np.zeros((nx,nt))
Q = np.zeros((nx,nt))

h[:,0] = Hi[:,0]
Q[:,0] = Qi[:,0]



for i in t_lim:
    for j in x_lim:
        if j == x_lim[0]:
            QB=Q[j+1,i-1]
            HB=h[j+1,i-1]
            CN = QB - Ca*HB - R * dt *abs(QB)*QB
            h[j,i] = h[j,0]
            Q[j,i] = CN + Ca * h[j,i] 

        elif j == x_lim[-1]:
            QA=Q[j-1,i-1]
            HA=h[j-1,i-1]
            CP = QA + Ca*HA - R * dt *abs(QA)*QA
            Q[j,i] = 0
            h[j,i] = (CP - Q[j,i])/Ca
            
        else:
            QA=Q[j-1,i-1]
            HA=h[j-1,i-1]
            QB=Q[j+1,i-1]
            HB=h[j+1,i-1]
            CP = QA + Ca*HA - R * dt *abs(QA)*QA
            CN = QB - Ca*HB - R * dt *abs(QB)*QB
            Q[j,i] = 0.5*(CP + CN)
            h[j,i] = (CP - Q[j,i])/Ca

# %% Extract Data from matrix and plot

Hvalve = h[-1,:]

plt.plot(t, Hvalve)
plt.xlabel("Time (s)")
plt.ylabel("Pressure head (m)")
plt.title('MOC - Simple')
plt.show()





