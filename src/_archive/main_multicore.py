#===========================================================
# Program written for  Combustion Dynamics - Fall 2025
#                        PROJECT
# Authors:
#   --   PL   --
#   --      --
#===========================================================

import sys
import numpy as np
import pandas as pd
import cantera as ct
from time import time
import argparse
import matplotlib.pyplot as plt

from multiprocessing import Pool
from functools import partial

# User imports
from inputs import *
import modules.thermodynamics as thermo
from modules.utilities import *
from modules.borghi import BorghiPlot


# set mechanism
mech:str    = 'ar22.yaml'

# Set fuel and oxidizer
fuel:dict       = {"CH4":1}
oxidizer:dict   = {"O2":0.21 ,"N2":0.79}


def evaluate_design_point(U0, PHI, xM):
    """Evaluate a single design point"""
    
    # Initialize gas object
    gas = ct.Solution(mech)
    gas.TP = T_30, P_30
    gas.set_equivalence_ratio(phi=PHI, fuel=fuel, oxidizer=oxidizer)
    
    # Store initial properties
    rho_initial = gas.density_mass
    visc_initial = gas.viscosity
    
    # Calculate adiabatic temperature
    gas.equilibrate('HP', solver='gibbs', max_steps=10000)
    t_ad = gas.T
    
    # Mass flow and diameter calculations
    m_dot_air = (AFR_STOIC * mdot_fuel_total) / (PHI * n_premixers)
    D = np.sqrt((m_dot_air * 4) / (rho_initial * np.pi * U0))
    
    # Reynolds number and turbulence intensity
    Re = (rho_initial * U0 * D) / visc_initial
    I0 = 0.16 * Re**(-1/8)
    
    # Reset gas for flame calculation
    gas.TP = T_30, P_30
    gas.set_equivalence_ratio(phi=PHI, fuel=fuel, oxidizer=oxidizer)
    
    # Compute laminar flame speed
    flame = ct.FreeFlame(gas, width=D)
    flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)
    flame.solve(loglevel=0, auto=True)
    S_L = flame.velocity[0]
    
    # Compute flame thickness
    x = flame.grid
    T = flame.T
    delta_f = (T[-1] - T[0]) / (max(np.gradient(T, x)))
    
    # Borghi diagram parameters
    lt_delta_f = (0.1 * D) / delta_f
    u_prime_S_L = (0.035 * U0) / S_L
    
    return {
        'U0': U0,
        'PHI': PHI,
        'xM': xM,
        'T_ad': t_ad,
        'm_dot_air': m_dot_air,
        'D': D,
        'S_L': S_L,
        'delta_f': delta_f,
        'Re': Re,
        'I0': I0,
        'borghi_x': lt_delta_f,
        'borghi_y': u_prime_S_L
    }


def main(args):
    """Parallel version - much faster for many design points"""
    from multiprocessing import Pool
    
    startTime = time()
    
    U0_range    = np.arange(40, 241, 10)            # Bulk flow velocity
    PHI_range   = np.arange(0.45, 0.48, 0.01)      # Equivalence ratio
    xM_range = np.array([24])

    
    params_list = [(u, p, x) for u in U0_range for p in PHI_range for x in xM_range]

    print(f"Processing {len(params_list)} design points in parallel...")
    
    with Pool() as pool:
        results = pool.starmap(evaluate_design_point, params_list)
    
    df = pd.DataFrame(results)
    print(df)
    
    if args.data:
        df.to_csv('test_data.csv', index=False)
    
    print_stats(startTime, time())
    return df


if __name__ == "__main__":

    # Main arguments
    parser = argparse.ArgumentParser(prog='combustor', description='Cantera model for premixed flames.')

    # Settings
    parser.add_argument('-d', '--data',    action='store_true', help='Save all output data to csv.')
    parser.add_argument('-p', '--plot',    action='store_true', help='Save all plots.')

    arguments = parser.parse_args()
    main(arguments)