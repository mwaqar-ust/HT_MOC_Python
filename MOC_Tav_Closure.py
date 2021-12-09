# Packages
import math
import numpy as np
import matplotlib.pyplot as plt


# %% constant variables
nu = 0
fDW = 0
g = 9.81
K = 2*10**9
rho = 1000

# %% Pipe Properties
L = 100
Di = 0.15
e = 0.0025
A = math.pi*0.25*Di**2
E = 120*10**9
a = (rho*(1/K+(1-nu**2)*2/E*(Di/2)/e))**(-0.5)

# %% Computation Constant Parameters
Ca = g*A/a
Z = 1/Ca
R = fDW/(2*Di*A)
mio=fDW/(2*g*Di*A**2)


# %% Steady-state pressure head and flow-rate at the downstream end
Qv = 2/1000
Hv = 40
V0 = Qv/A
HJ = a*V0/g;
# %% Data Sampling & Discretization
Fs = 250
dt = 1/Fs
period = 4*L/a
T = 5*period
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
    
# %% Tau Curve for Valve Closure
Tclose = 0.04

closure_type = "pulse" # linear, nlinear, gaussian, pulse

if closure_type == "linear":
    t1 = np.arange(0,Tclose+dt,dt)
    tav1 = (1 - t1/Tclose)
    Tav = np.zeros((nt))
    Tav[0:len(tav1)] = tav1
elif closure_type == "nlinear":
    t1 = np.arange(0,Tclose+dt,dt)
    tav1 = (1 - t1/Tclose)**2
    Tav = np.zeros((nt))
    Tav[0:len(tav1)] = tav1
elif closure_type == "gaussian":
    sig     = np.sqrt(-Tclose**2/2/np.log(0.001))
    Tav     = np.exp(-t**2/2/sig**2)
elif closure_type == "pulse":
    Tav = -1* abs(np.exp(-0.5 * ((t - Tclose)/((Tclose/5)))**2))+1
    # However, this is giving some dissipation error...!!!
    
    

# %% Define Matrices
h = np.zeros((nx,nt))
Q = np.zeros((nx,nt))

h[:,0] = Hi[:,0]
Q[:,0] = Qi[:,0]



for i in t_lim:
    print(np.round(i/nt*100,2))
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
            
            Qo_tav = Q[j,0]*Tav[i]
            hout = 0
            Aev = Qo_tav/np.sqrt(h[j,0])
            h[j,i] = CP/Ca + (Aev*(Aev - (Aev**2 - 4*hout*Ca**2 + 4*CP*Ca)**(1/2)))/(2*Ca**2)
            Q[j,i] = CP - Ca * h[j,i] 
            
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

Hvalve = (h[-1,:]-h[-1,0])/HJ

plt.plot(t, Hvalve)
plt.xlabel("Time (s)")
plt.ylabel("Pressure head (m)")
plt.title('MOC - Simple')
plt.show()





