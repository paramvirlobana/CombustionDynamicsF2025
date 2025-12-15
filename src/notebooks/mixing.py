import os
import numpy as np
import pandas as pd
import cantera as ct
from time import time
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

mech = 'ar22.yaml'
fuel:dict       = {"CH4":1}
oxidizer:dict   = {"O2":0.21 ,"N2":0.79}

gas_a = ct.Solution(mech)
gas_a.TPX = 750, 24 * ct.one_atm, oxidizer
rho_a = gas_a.density

gas_b = ct.Solution(mech)
gas_b.TPX = 300.0, ct.one_atm, fuel
rho_b = gas_b.density

res_a = ct.Reservoir(gas_a, name='Air Reservoir')
res_b = ct.Reservoir(gas_b, name='Fuel Reservoir')
downstream = ct.Reservoir(gas_a, name='Outlet Reservoir')

gas_b.TPX = 300.0, ct.one_atm, 'O2:0.21, N2:0.79'
mixer = ct.IdealGasReactor(gas_b, name='Mixer')

mfc1 = ct.MassFlowController(res_a, mixer, mdot=61, name="Air Inlet")
mfc2 = ct.MassFlowController(res_b, mixer, mdot=1.66, name="Fuel Inlet")

outlet = ct.Valve(mixer, downstream, K=10.0, name="Valve")

sim = ct.ReactorNet([mixer])
sim.advance_to_steady_state()

print(mixer.thermo.report())

try:
    diagram = sim.draw(print_state=True, species="X")
except ImportError as err:
    print(f"Unable to show network structure:\n{err}")