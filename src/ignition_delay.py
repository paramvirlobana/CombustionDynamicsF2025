import numpy as np
import matplotlib.pyplot as plt

import time

import cantera as ct
print('Running Cantera version: ' + ct.__version__)

from inputs import *

gas = ct.Solution(os.path.join(MODELS, "ar22.yaml"))
gas.set_equivalence_ratio(phi=1.0, fuel={'CH4':1.0}, oxidizer={"O2":0.21, "N2":0.79})
gas.TP = T_30, P_30


# init the reactor
r = ct.IdealGasConstPressureMoleReactor(gas, name='R1')
sim = ct.ReactorNet([r])
sim.derivative_settings = {"skip-third-bodies":True, "skip-falloff":True}
sim.preconditioner = ct.AdaptivePreconditioner() # Manage the numerical convergence / iteration by "conditioning" the matrix.


states = ct.SolutionArray(gas, extra='time')

T_old = 0
T_new = 0
t = 0.0
maxtime = 1e+6
counter = 0

while t < maxtime:
    T_old = r.T
    t = sim.step()
    T_new = r.T

    if(counter%1 == 0):
        states.append(r.thermo.state, time=t)

        if (np.absolute(T_new - T_old) < 0.0001 and T_old > (T_30 + 800)):
            print(f"T30 = {T_30:4.2f} K, Tf = {r.T:.2f} K and Ignition Delay = {t:.4f} seconds")
            break

        counter += 1