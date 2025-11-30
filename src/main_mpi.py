#===========================================================
# Program written for  Combustion Dynamics - Fall 2025
#                        PROJECT
# Authors:
#   --   Paramvir Lobana   --
#   --   Fouad Al Laham    --
#===========================================================

import sys
import numpy as np
import pandas as pd
import cantera as ct
from time import time
import argparse
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
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


def main(args):

    startTime = time()

    # Need to define a loop that can iterate over several input parameters.
    # What can we define as the input variable?

    # NOTE: Design Variables
    U0_range    = np.arange(60, 180, 10)            # Bulk flow velocity
    PHI_range   = np.arange(0.45, 0.48, 0.01)      # Equivalence ratio
    M_range     = np.arange(0.005, 0.02, 0.004)         # Turbulence intensity (as a fraction of U0)
    L_range     = np.arange(0.06, 0.1201, 0.02)        # Premixer length

    design_variables:tuple = (U0_range, PHI_range, M_range, L_range) # order matters.

    design_points = len_product(*design_variables)
    
    print(f"Total design points: {design_points}")
    print(f"Using {cpu_count()} CPU cores for parallel processing") 
    print("------------------")
    
    points_per_u0 = len(PHI_range) * len(M_range) * len(L_range)
    design_params = [(U0, i * points_per_u0) for i, U0 in enumerate(U0_range)]
    
    process_func = partial(computeDesignPoints, PHI_range=PHI_range, M_range=M_range, L_range=L_range)
    

    # For multiprocessing
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(process_func, design_params)
    
    # Combine results from all processes
    data = {
        key: np.concatenate([r[key] for r in results])
        for key in results[0].keys()
    }

    df = pd.DataFrame({
        'U0': data['U0'],
        'PHI': data['PHI'],
        'M' : data['M'],
        'L' : data['L'],
        'T_ad': data['T_ad'],
        'm_dot_air': data['m_dot_air'],
        'D': data['D'],
        'tau_mix': data['tau_mix'],
        'eta': data['eta'],
        'L_premixer': data['L_premixer'],
        'x_injector': data['x_injector'],
        'S_L': data['S_L'],
        'delta_f': data['delta_f'],
        'lt_Lf': data['lt_Lf'],
        'ui_SL': data['ui_SL'],
        'ST_SL_Peters': data['ST_SL_Peters'],
        'ST_SL_Zimont': data['ST_SL_Zimont'],
        'dP_grid': data['dP_grid'],
        'dP_fric': data['dP_fric'],
        'dP_total': data['dP_total'],
        'dP': data['dP'],
    })
    
    if args.verbose:
        print(df)

    if args.data:
        df.to_csv('test_data.csv', index=False)

    print_stats(startTime, time())



def computeDesignPoints(design_params, PHI_range, M_range, L_range):
    """
    Process a single design point for a given U0 value.
    This function handles all PHI and xM combinations for one U0.
    """
    U0, iter_start = design_params
    
    # Calculate number of points for this U0
    n_points = len(PHI_range) * len(M_range)
    
    results = {
        'U0':           np.zeros(n_points),
        'PHI':          np.zeros(n_points),
        'M':            np.zeros(n_points),
        'L':            np.zeros(n_points),
        'T_ad':         np.zeros(n_points),
        'm_dot_air':    np.zeros(n_points),
        'D':            np.zeros(n_points),
        'L_premixer':   np.zeros(n_points),
        'tau_mix':      np.zeros(n_points),
        'eta':          np.zeros(n_points),
        'S_L':          np.zeros(n_points),
        'delta_f':      np.zeros(n_points),
        'lt_Lf':        np.zeros(n_points),
        'ui_SL':        np.zeros(n_points),
        'x_injector':   np.zeros(n_points),
        'dP_grid':      np.zeros(n_points),
        'dP_fric':      np.zeros(n_points),
        'dP_total':     np.zeros(n_points),
        'dP':           np.zeros(n_points),
        'ST_SL_Peters': np.zeros(n_points),
        'ST_SL_Zimont': np.zeros(n_points)
    }
    
    # design iterator
    for (PHI, M, L), (i_phi, i_M, i_L) in product(PHI_range, M_range, L_range):
        iter = i_phi * (len(M_range) * len(L_range)) + i_M * len(L_range) + i_L

        results['U0'][iter]     = U0
        results['PHI'][iter]    = PHI
        results['M'][iter]      = M
        results['L'][iter]      = L 

        # Initialize gas object for equilibrium calculations
        gas_A = ct.Solution(mech)
        gas_A.TP = T_30, P_30
        gas_A.set_equivalence_ratio(phi=PHI, fuel=fuel, oxidizer=oxidizer)

        # Calculate t_ad for each design.
        gas_A.equilibrate('HP', solver='gibbs', max_steps=10000)
        results['T_ad'][iter] = gas_A.T

        # Using equivalence ratio to get mass flow of air. 
        # NOTE: This is the mass flow of air in each premix chamber.
        results['m_dot_air'][iter] = (AFR_STOIC * (mdot_fuel_total/n_premixers)) / (PHI)

        # Using conservation of mass, we can calculate the required diameter.
        results['D'][iter] = np.sqrt((results['m_dot_air'][iter] * 4) / (gas_A.density_mass * np.pi * U0))



        # NOTE Find premixer length based on turbulence intensity decay
        Irange = np.sqrt(y1) * 100
        I_to_x = interp1d(Irange, x1, kind='linear', fill_value='extrapolate')
        x_at_val = I_to_x(2.95)
        idx = np.argmin(np.abs(x1 - x_at_val))

        cutoff_xM = x1[idx]
        results['L_premixer'][iter] = np.round(cutoff_xM * m_grid_ig, 6)
        results['tau_mix'][iter] = results['L_premixer'][iter] / U0

        # Calculate the mixing timescale (This is objective 3).
        #eta_range = np.linspace(0.02, 0.29, 1001)
        #sigma_range = f2(eta_range)
        #max_idx = np.argmax(sigma_range)
        
        #results['x_injector'][iter] = (results['D'][iter]) / (2 * eta_range[max_idx])

        # Calculate the actual eta for the current design.
        # NOTE Later in the post processing step, if eta > 0.28, we drop the design point.
        results['eta'][iter] = (results['D'][iter] / 2) / results['L_premixer'][iter]
        
        # New gas object for pressure loss calculations.
        gas_B = ct.Solution(mech)
        gas_B.TP = T_30, P_30
        gas_B.set_equivalence_ratio(phi=PHI, fuel=fuel, oxidizer=oxidizer)

        # Calculate pressure losses using the pressure_loss function
        pressure_results = pressureLoss(gas_B, U0, results['D'][iter], results['L_premixer'][iter])
        
        results['dP_grid'][iter] = pressure_results['dP_grid']
        results['dP_fric'][iter] = pressure_results['dP_fric']
        results['dP_total'][iter] = pressure_results['dP_total']
        results['dP'][iter] = pressure_results['dP']

        # Make sure gas_B is reset for flame speed calculations.
        gas_B.TP = T_30, P_30
        gas_B.set_equivalence_ratio(phi=PHI, fuel=fuel, oxidizer=oxidizer)

        flame = ct.FreeFlame(gas_B, width=results['D'][iter])
        flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)
        flame.solve(loglevel=0, auto=True)
        results['S_L'][iter] = flame.velocity[0]

        # Compute the flame thickness
        x = flame.grid
        T = flame.T
        results['delta_f'][iter] = (T[-1] - T[0]) / (max(np.gradient(T, x)))

        # Plot results on the Borghi diagram.
        results['lt_Lf'][iter] = (0.1 * results['D'][iter]) / results['delta_f'][iter]
        results['ui_SL'][iter] = (0.030 * U0) / results['S_L'][iter]

        # NOTE Objective 5a
        # Using corelations to calculate the turbulent flame speed.
        # Burner Rim Stabilized Flames

        # Calculating the turbulent Re
        # TODO to be verified because using a dummy value for intensity at the moment.
        Re_T = (0.030 * U0 * 0.1 * results['D'][iter]) / (results['S_L'][iter] * results['delta_f'][iter])

        # NOTE Peters Model

        A = (Re_T / (results['ui_SL'][iter])) + 1
        E1 = 1 * (a_peters/2) * A
        results['ST_SL_Peters'][iter] = -1 * E1 + np.sqrt((E1)**2 + (a_peters * A * results['ui_SL'][iter]) + a_peters + 1)

        # NOTE Zimont Model
        results['ST_SL_Zimont'][iter] = A_zimont * (Pr_zimont)**(1/4) * (Re_T)**(1/4) * np.sqrt(results['ui_SL'][iter])

        global_iter = iter_start + iter
        if (global_iter+1) % 4 == 0: # we print every 4th iteration.
            print(f"Iteration {global_iter+1:>4} | U0 = {U0:>7.2f} m/s | PHI = {PHI:>6.3f} | M = {M:>7.4f}")
            print(f"Iteration {global_iter+1:>4}| ST/SL P = {results['ST_SL_Peters'][iter]:>7.4f} | ST/SL Z = {results['ST_SL_Zimont'][iter]:>7.4f}")
            print("------------------")
    
    return results



def postProcessign():
    """
    Author: PL
    Returns:
    -------
    Processed dataframe based on design conditions.
    """
    pass

def pressureLoss(gas, U0, D_val, L_val):
    """
    Author: FL
    
    Returns:
    --------
    dict : Dictionary containing pressure loss results
        - Re: Reynolds number
        - f_D: Darcy friction factor
        - dP_grid: Grid pressure drop [Pa]
        - dP_fric: Frictional pressure drop [Pa]
        - dP_total: Total pressure drop [Pa]
    """
    # NOTE Pressure Loss Calculations
    rho = gas.density_mass
    mu = gas.viscosity

    Re = (rho * U0 * D_val) / mu

    # Darcy friction factor (Blasius)
    # valid for 3000 < Re < 100,000 but used widely up to 10^5–10^6 for estimates
    f_D = 0.3164 * Re ** (-0.25)

    # Grid solidity (NOTE ASSUMPTION: sigma = 0.40)
    sigma = 0.40

    # Grid drag coefficient (NOTE ASSUMPTION: C_D = 1.5)
    C_D = 1.5

    # Pressure drop across the grid
    DeltaP_grid = 0.5 * rho * U0**2 * C_D * sigma

    # Frictional pressure drop along premixer length L
    DeltaP_fric = f_D * (L_val / D_val) * 0.5 * rho * U0**2

    # Total pressure drop through injector/premixer
    DeltaP_total = DeltaP_grid + DeltaP_fric
    dP_percent = (DeltaP_total / (P_30)) * 100

    return {
        'dP_grid': DeltaP_grid,
        'dP_fric': DeltaP_fric,
        'dP_total': DeltaP_total,
        'dP': dP_percent,
    }


if __name__ == "__main__":

    # Main arguments
    parser = argparse.ArgumentParser(prog='combustor', description='Cantera model for premixed flames.')

    # Settings
    parser.add_argument('-d', '--data',    action='store_true', help='Save all output data to csv.')
    parser.add_argument('-p', '--plot',    action='store_true', help='Save all plots.')
    parser.add_argument('-v', '--verbose', action='store_true', help="Prints results to console if flag invoked.")

    arguments = parser.parse_args()
    main(arguments)