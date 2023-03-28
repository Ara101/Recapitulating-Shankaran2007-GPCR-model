
Here I Recapitulate the mechanistic GPCR model in the paper “Receptor downregulation and desensitization enhance the information processing ability of signalling receptors,”. doi: 10.1186/1752-0509-1-48.


The main aspects of the model are the system of ODEs:

    dR/dt = –konRL + koffC 
    dL/dt =  (–kon RL + koffC)/(NavV) + f(t)
    dC/dt = kon RL – koffC − kfrC + krrCa
    dCa/dt = kfrC − krrCa − kdsCa
    dCd/dt = kds Ca
    dG/dt = −kaGCa + kiGa 
    dGa/dt = kaGCa − kiGa

The parameters and inital conditions:

    kon = 8.4e7 /M/s
    koff = 0.37 /s
    kfr = 10 /s
    krr = 10 /s
    kds = 0.065 /s
    ka = 10e−7 /s
    ki = 2e−1 /s
    RT = 5.5e4
    GT = 1e5
    V = 4e−10 liters/cell


I attempted recreate this model with ODE solvers and multiple modeling software the result was a partial recreation with most methods
