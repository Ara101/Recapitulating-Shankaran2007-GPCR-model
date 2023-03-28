
Here I Recapitulate the mechanistic GPCR model in the paper “Receptor downregulation and desensitization enhance the information processing ability of signalling receptors,” .

The main aspects of the model are the system of ODEs:

    dR/dt = –konRL + koffC 
    dL/dt =  (–kon RL + koffC)/(NavV) + f(t)
    dC/dt = kon RL – koffC − kfrC + krrCa
    dCa/dt = kfrC − krrCa − kdsCa
    dCd/dt = kds Ca
    dG/dt = −kaGCa + kiGa 
    dGa/dt = kaGCa − kiGa

The parameters and inital conditions:

    kon = 8.4×107 /M/s
    koff = 0.37 /s
    kfr = 10 /s
    krr = 10 /s
    kds = 0.065 /s
    ka = 10−7 /s
    ki = 2×10−1 /s
    RT = 5.5×104
    GT = 1×105
    V = 4×10−10 liters/cell


I attempted recreate this model with ODE solvers and multiple modeling software the result was a partial recreation with most methods
