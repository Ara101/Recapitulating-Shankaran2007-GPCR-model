"""Recapitulating the model: Receptor downregulation 
and desensitization enhance the information processing ability of signalling receptors"""
# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tellurium as te
import numpy as np 

# The intial EGFR Model
shankaran_model_EGFR = te.loada("""
Ligand + Surface_receptor -> Surface_complex; Ligand*Surface_receptor*kon
Surface_complex -> Ligand + Surface_receptor; Surface_complex*koff

Interalized_receptor -> Surface_receptor; Interalized_receptor*Vr
Surface_receptor -> Interalized_receptor; Surface_receptor*kt

Surface_complex -> Internalized_complex; Surface_complex*ke

Ligand  = 0.01*Kd + 5961
Surface_receptor = Rt
Surface_complex = 0

kon = 0.097
koff = 0.24
Kd = koff/kon

kt = 0.02
ke = 0.15
Vr = kt*Rt 

Rt = 200000
V = 0.0000000004
""")

shankaran_model_EGFR.reset()
m = shankaran_model_EGFR.simulate(0, 1000, 100)
plt.figure(1) 
shankaran_model_EGFR.plot()

# change a species value to recreate the graph
shankaran_model_EGFR.resetAll()
count = 0
plt.figure(2) 
for i in [1,2,5,7.6,10,20]: 
    shankaran_model_EGFR.resetToOrigin()
    shankaran_model_EGFR.ke = 0.15*i
    result = shankaran_model_EGFR.simulate(0, 50, 100,['time','Surface_complex'])
    plt.plot(result[:,0],result[:,1], label = '%s ke = '%i) 
    plt.title('EGFR w varied ke')
    count += 1

plt.legend()
plt.show()



#The GPCR Model
# seperate ligand and receptor GPCR
shankaran_model_GPCR = te.loada(""" 
$L + $R -> C; L*R*kon
C -> R; C*koff
C -> L; koff*C

C -> Ca; C*kfr
Ca -> C; Ca*krr

Ca -> Cd; Ca*kds

Ga -> G; Ga*ki
G -> Ga; G*ka
Ca + G -> Ga; ka*G*Ca

R = RT
G = GT
C = 0
Ca = 0
Cd =0
Ga = 0

kon = 8.4*10^7
koff = 0.37
Kd = koff/kon
kfr = 10
krr = 10
kds = 0.065
ka = 10^-7
ki = 2*10^-1

Nav = 6.022*10^23
V = 4*10^-10
RT = (5.5*10^4)
GT = (1*10^5)

at (time > 1): $L = 0.01*Kd
""")
shankaran_model_GPCR.resetToOrigin()
GPCR_model = shankaran_model_GPCR.simulate(0,250,100)
plt.figure(3) 
shankaran_model_GPCR.plot()

shankaran_model_GPCR.resetAll()
plt.figure(4) 
for i in [10**-2, 10**-1, 10**0, 10**1, 10**2]:
    shankaran_model_GPCR.resetToOrigin()
    shankaran_model_GPCR.kds = 0.065*i
    result = shankaran_model_GPCR.simulate(0, 250, 200,['time','Ga'])
    print(max(result[:,1]))
    print("Kds = ", shankaran_model_GPCR.kds)
    plt.plot(result[:,0],result[:,1]*1000, label = '%s kds = ' %i)
    plt.title('GPCR with varied kds')

plt.legend()
plt.show()

# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tellurium as te
import numpy as np 
from scipy.integrate import odeint
#The GPCR Model with differential eq.
# seperate ligand and receptor GPCR

# funtction that returns dX/dt
#The GPCR Model with differential eq.
# seperate ligand and receptor GPCR

# funtction that returns dX/dt
def model(X, t):
    
    # parameters
    kon = 8.4e7
    koff = 0.37
    Kd = koff/kon
    kfr = 10
    krr = 10
    kds = 0.065
    ka = 1e-7
    ki = 2e-1
    Nav = 6.022e23
    V = 4e-10

    # mathematical form of the model
    R = X[0]
    L = X[1]
    C = X[2]
    Ca = X[3]
    Cd = X[4]
    G = X[5]
    Ga = X[6]
    ft = 0 # or  0.01*Kd

    dRdt = -koff*R*L + koff*C
    dLdt = (-koff*R*L + koff*C)/(Nav*V) + ft
    dCdt = koff*R*L - koff*C - kfr*C + krr*Ca
    dCadt = kfr*C - krr*Ca - kds*Ca
    dCddt = kds*Ca
    dGdt = -ka*G*Ca + ki*Ga
    dGadt = ka*G*Ca - ki*Ga

    return [dRdt, dLdt, dCdt, dCadt, dCddt, dGdt, dGadt]

# Inital conidtions

Kd = 4.404761904761905e-09
ft = 0.01*Kd

RT = (5.5e4)
GT = (1e5)
R = RT
G = GT
L = ft
C = 0
Ca = 0
Cd = 0
Ga = 0
IC = [R, L, C, Ca, Cd, G, Ga]

# simulation time points
t = np.linspace(0,250,100)

# solve ODE
results = odeint(model, IC, t)

plt.figure(5)
plt.plot(t,results)
plt.xlabel('time')
plt.ylabel('Concentation')
plt.legend(['R', 'L', 'C', 'Ca', 'Cd', 'G', 'Ga'])
plt.title('GPCR with ODE solver')
plt.show()

plt.figure(6)
species = IC.index(Ga)
plt.plot(t,results[:,species])
plt.xlabel('time')
plt.ylabel('Concentation')
plt.title('GPCR w ODE solver and Ga focus')
plt.show()