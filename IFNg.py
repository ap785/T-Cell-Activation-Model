import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import seaborn as sns

# Model without costimulation
def model_no_costim(x, t, params):
    kon, koff, ksyn, kbasal, kdown, kact, kdeg = params

    # assign each ODE to a vector element
    C = x[0]
    R = x[1]
    M = x[2] # Signal from TCR activation
    O = x[3]
    L = LT - C

    # Define ODEs (without costimulation)
    dRdt = ksyn - (kbasal * R) - (kon * L * R) + (koff * C)
    dCdt = (kon * L * R) - ((koff + kdown) * C)
    dMdt = kact * (1 - M) * (np.heaviside(C - Cstar, 1)) - kdeg * M # (1-M) is self-regulating (inhibitory) mechanism - as M increases contribution of kact to M production decreases
    dOdt = kprod * M

    return [dCdt, dRdt, dMdt, dOdt]


#####################################################################################
# Model for CD2 Costimulation
#####################################################################################

# Model with CD2 costimulation (affects TCR-pMHC binding )
def model_with_cd2(x, t, params):
    kon, koff, ksyn, kbasal, kdown, kact, kdeg, kon_2, koff_2, kcd2a, kcd2b, kdown_2, kbasal_2, ksyn_2 = params

    # assign each ODE to a vector element
    C = x[0]
    R = x[1]
    M = x[2]
    O = x[3]
    R2 = x[4]
    C2 = x[5]
    L = LT - C

    dRdt = ksyn - (kbasal * R) - (kon * L * R) + (koff * C) - (kcd2a * C2) * (kon * L * R) # have to offset production of R 
    dCdt = (kon * L * R) - ((koff + kdown) * C) + (kcd2a * C2) * (kon * L * R)  # TCR-pMHC complex formation modulated by costimulation
    dMdt = kact * (1 - M) * (np.heaviside(C - Cstar, 1)) - kdeg * M
    dOdt = kprod * M *(1 + kcd2b*C2)
    dR2dt = ksyn_2 - (kbasal_2 * R2) - (kon_2 * L2 * R2) + (koff_2 * C2)
    dC2dt = (kon_2 * L2 * R2) - ((koff_2 + kdown_2) * C2)

    return [dCdt, dRdt, dMdt, dOdt, dR2dt, dC2dt]


#####################################################################################
# Model for ICAM Costimulation
#####################################################################################

def model_with_icam(x, t, params):
    kon, koff, ksyn, kbasal, kdown, kact, kdeg, kon_2, koff_2, kicama, kicamb, kdown_2, kbasal_2, ksyn_2 = params

    # assign each ODE to a vector element
    C = x[0]
    R = x[1]
    M = x[2]
    O = x[3]
    R2 = x[4]
    C2 = x[5]
    L = LT - C

    dRdt = ksyn - (kbasal * R) - (kon * L * R) + (koff * C) - (kicama * C2) * (kon * L * R) # have to offset production of R 
    dCdt = (kon * L * R) - ((koff + kdown) * C) + (kicama * C2) * (kon * L * R)  # TCR-pMHC complex formation modulated by costimulation
    dMdt = kact * (1 - M) * (np.heaviside(C - Cstar, 1)) - kdeg * M
    dOdt = kprod * M * (1 + kicamb*C2)
    dR2dt = ksyn_2 - (kbasal_2 * R2) - (kon_2 * L2 * R2) + (koff_2 * C2)
    dC2dt = (kon_2 * L2 * R2) - ((koff_2 + kdown_2) * C2)

    return [dCdt, dRdt, dMdt, dOdt, dR2dt, dC2dt]


#####################################################################################
# Model for CD28 Costimulation
#####################################################################################

# Model with CD28 costimulation
def model_with_cd28(x, t, params):
    kon, koff, ksyn, kbasal, kdown, kact, kdeg, kon_2, koff_2, k28, kdown_2, kbasal_2, ksyn_2 = params

    # assign each ODE to a vector element
    C = x[0]
    R = x[1]
    M = x[2]
    O = x[3]
    R2 = x[4]
    C2 = x[5]
    L = LT - C
    
    dRdt = ksyn - (kbasal * R) - (kon * L * R) + (koff * C)
    dCdt = (kon * L * R) - ((koff + kdown) * C) 
    dMdt = kact * (1 - M) * (np.heaviside(C - Cstar, 1)) - kdeg * M
    dOdt = kprod * (1 + (k28 * C2)) * M # Costimulation integrating after threshold switch at kprod. Kprod scaled by
    dR2dt = ksyn_2 - (kbasal_2 * R2) - (kon_2 * L2 * R2) + (koff_2 * C2)
    dC2dt = (kon_2 * L2 * R2) - ((koff_2 + kdown_2) * C2)

    return [dCdt, dRdt, dMdt, dOdt, dR2dt, dC2dt]

#####################################################################################
# Model for CD27 Costimulation
#####################################################################################

def model_with_cd27(x, t, params):
    kon, koff, ksyn, kbasal, kdown, kact, kdeg, kon_2, koff_2, kdown_2, kbasal_2, ksyn_2, Cstar_2 = params

    # assign each ODE to a vector element
    C = x[0]
    R = x[1]
    M = x[2]
    O = x[3]
    R2 = x[4]
    C2 = x[5]
    L = LT - C
    Y = Cstar / (1 + (C2 / Cstar_2))
    
    dRdt = ksyn - (kbasal * R) - (kon * L * R) + (koff * C)
    dCdt = (kon * L * R) - ((koff + kdown) * C) 
    dMdt = kact * (1 - M) * (np.heaviside(C - Y, 1)) - kdeg * M
    dOdt = (kprod * M)
    dR2dt = ksyn_2 - (kbasal_2 * R2) - (kon_2 * L2 * R2) + (koff_2 * C2)
    dC2dt = (kon_2 * L2 * R2) - ((koff_2 + kdown_2) * C2)

    return [dCdt, dRdt, dMdt, dOdt, dR2dt, dC2dt]

#####################################################################################
# Model for 4-1BB Costimulation
#####################################################################################

# Model with 4-1BB costimulation - affects downstream signalling (switch)
def model_with_41bb(x, t, params):
    kon, koff, ksyn, kbasal, kdown, kact, kdeg, kon_2, koff_2, kdown_2, kbasal_2, ksyn_2, time_delay, Cstar_2 = params

    # assign each ODE to a vector element
    C = x[0]
    R = x[1]
    M = x[2]
    O = x[3]
    R2 = x[4]
    C2 = x[5]
    mRNA = x[6]
    L = LT - C
    Y = Cstar / (1 + (C2 / Cstar_2))
    # Cstar_2 represents the number of 4-1BB molecules needed to lower threshold of cstar

    dRdt = ksyn - (kbasal * R) - (kon * L * R) + (koff * C)
    dCdt = (kon * L * R) - ((koff + kdown) * C) 
    dMdt = (kact * (1 - M) * (np.heaviside(C - Y, 1)) - kdeg * M)
    dOdt = kprod * M
    # 4-1BB Upregulation dependent on formation of C from T cell activation. 
    # Introduce a time delay using a Heaviside function for 41BB upregulation.
    dR2dt = (ksyn_2 * mRNA *  np.heaviside(t - time_delay, 1) * O) + (koff_2 * C2) - (kbasal_2 * R2) - (kon_2 * L2 * R2)
    dmRNAdt = np.heaviside(C - Y, 1)
    dC2dt = (kon_2 * L2 * R2) - ((koff_2 + kdown_2) * C2)

    return [dCdt, dRdt, dMdt, dOdt, dR2dt, dC2dt, dmRNAdt]

#####################################################################################
# Model Initialisation
#####################################################################################

# Define initial Conditions

x0_no_costim = [0, 100, 0, 0]  # No R2, L2, C2 for no costimulation
x0_with_cd2 = [0, 100, 0, 0, 100, 0]
x0_with_icam = [0, 100, 0, 0, 100, 0]
x0_with_cd28 = [0, 100, 0, 0, 100, 0] 
x0_with_cd27 = [0, 100, 0, 0, 100, 0]  
x0_with_41bb = [0, 100, 0, 0, 0, 0, 0]  # N.B R2 starts at 0 for 41BB
L2 = 1000  # Constant ligand concentration for L2
LT = 1000
Cstar = 0.05
kprod = 0.00002

t = np.linspace(0, 20 * 3600, 100)

t_hours = t / 3600

# Define parameter sets for both models
                # kon, koff, ksyn, kbasal, kdown, kact, kdeg

params_no_costim_s = (0.0001, 1, 1e-2, 1, 0.1, 0.5, 0.2)

#kon, koff, ksyn, kbasal, kdown, kact, kdeg, kon_2, koff_2, kcd2a, kcd2b, kdown_2, kbasal_2, ksyn_2
params_with_cd2_s = (0.0001, 1, 1e-2, 1, 0.1, 0.5, 0.2, 0.005, 1, 18 * 3600, 0.004 * 3600, 0.1, 1, 0.02)

params_with_icam_s = (0.0001, 1, 1e-2, 1, 0.1, 0.5, 0.2, 0.001, 1, 5 * 3600, 0.01 * 3600, 0.1, 1, 0.02)

#k28 value lower than kcd2 / kbb
params_with_cd28_s = (0.0001, 1, 1e-2, 1, 0.1, 0.5, 0.2, 0.001, 1, 0.01 * 3600, 0.1, 1, 0.02)


params_with_cd27_s = (0.0001, 1, 1e-2, 1, 0.1, 0.5, 0.2, 0.001, 1, 0.1, 0.0001, 0.02, 0.001 * 3600 * 3600)


# kon, koff, ksyn, kbasal, kdown, kact, kdeg,  
params_with_41bb_s = (0.0001, 1, 1e-2, 1, 0.1, 0.5, 0.2, 
                    # kon_2, koff_2, kdown_2, kbasal_2, ksyn_2, time_delay, Cstar_2
                    0.001, 1, 0.0001, 0.0001, 100, 9 * 3600 * 3600, 0.001 * 3600 * 3600)

params_no_costim = [param / 3600 for param in params_no_costim_s]
params_with_cd2 = [param / 3600 for param in params_with_cd2_s]
params_with_icam = [param / 3600 for param in params_with_icam_s]
params_with_cd28 = [param / 3600 for param in params_with_cd28_s]
params_with_cd27 = [param / 3600 for param in params_with_cd27_s]
params_with_41bb = [param / 3600 for param in params_with_41bb_s]

# Solve ODEs for both models
x_no_costim = odeint(model_no_costim, x0_no_costim, t, args=(params_no_costim,))
x_with_cd2 = odeint(model_with_cd2, x0_with_cd2, t, args=(params_with_cd2,))
x_with_icam = odeint(model_with_icam, x0_with_icam, t, args=(params_with_icam,))
x_with_cd28 = odeint(model_with_cd28, x0_with_cd28, t, args=(params_with_cd28,))
x_with_cd27 = odeint(model_with_cd27, x0_with_cd27, t, args=(params_with_cd27,))
x_with_41bb = odeint(model_with_41bb, x0_with_41bb, t, args=(params_with_41bb,))

cumulative_O_no_costim = np.cumsum(x_no_costim[:, 3])
cumulative_O_cd2 = np.cumsum(x_with_cd2[:, 3])
cumulative_O_icam = np.cumsum(x_with_icam[:, 3])
cumulative_O_cd28 = np.cumsum(x_with_cd28[:, 3])
cumulative_O_cd27 = np.cumsum(x_with_cd27[:, 3])
cumulative_O_41bb = np.cumsum(x_with_41bb[:, 3])
