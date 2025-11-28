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


def process_design_point(design_params, PHI_range, xM_range):

    U0, iter_start = design_params
    
    # Calculate number of points for this U0
    n_points = len(PHI_range) * len(xM_range)
    
    # Initialize arrays for this U0
    U0_arr      = np.zeros(n_points)
    phi_arr     = np.zeros(n_points)
    xM_arr      = np.zeros(n_points)
    t_ad        = np.zeros(n_points)
    m_dot_air   = np.zeros(n_points)
    D           = np.zeros(n_points)
    S_L         = np.zeros(n_points)
    delta_f     = np.zeros(n_points)
    lt_delta_f  = np.zeros(n_points)
    u_prime_S_L = np.zeros(n_points)
    L_min       = np.zeros(n_points)
    L           = np.zeros(n_points)
    ST_SL_P     = np.zeros(n_points)
    ST_SL_Z     = np.zeros(n_points)
    
    iter = 0
    
    for PHI in PHI_range:
        for xM in xM_range:
            # To store the iteration value.
            U0_arr[iter]      = U0
            phi_arr[iter]     = PHI
            xM_arr[iter]      = xM

            gas_A = ct.Solution(mech)
            gas_A.TP = T_30, P_30
            gas_A.set_equivalence_ratio(phi=PHI, fuel=fuel, oxidizer=oxidizer)

            # Calculate t_ad for each design.
            gas_A.equilibrate('HP', solver='gibbs', max_steps=10000)
            t_ad[iter] = gas_A.T

            # Using equivalence ratio to get mass flow of air. 
            # NOTE: This is the mass flow of air in each premix chamber.
            m_dot_air[iter] = (AFR_STOIC * mdot_fuel_total) / (PHI * n_premixers)

            # Using conservation of mass, we can calculate the required diameter.
            D[iter] = np.sqrt((m_dot_air[iter] * 4) / (gas_A.density_mass * np.pi * U0))

            # TODO
            # Need to add base code for objective 2 and objective 3.
            # .
            # .
            # .
            # .
            # .
            # TODO
                
            # Turbulence Grid and L_min
            # We calculate the minimum length required where turbulent intensity is above 3%.
            grid_range = np.linspace(20, 30, 1000)
            y_range = f1(grid_range)

            # Compute where turbulence drops below threshold.
            idx = np.where(y_range < 0.03**2)[0][0]
            cutoff_xM = grid_range[idx]

            L_min[iter] = np.round(cutoff_xM * m_grid_ig, 6)
            L[iter]     = xM * m_grid_ig


            # Calculate the mixing timescale (This is objective 3).
            eta = (D[iter] / 2) / L[iter]

            # Compute the laminar flame speed.
            # We initialize a new gas object because the previous object is
            # already solved and in equilibrium state.

            gas_B = ct.Solution(mech)
            gas_B.TP = T_30, P_30
            gas_B.set_equivalence_ratio(phi=PHI, fuel=fuel, oxidizer=oxidizer)

            flame = ct.FreeFlame(gas_B, width=D[iter])
            flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)
            flame.solve(loglevel=0, auto=True)
            S_L[iter] = flame.velocity[0]


            # Compute the flame thickness
            x = flame.grid
            T = flame.T
            delta_f[iter] = (T[-1] - T[0]) / (max(np.gradient(T, x)))


            # Plot results on the Borghi diagram.
            lt_delta_f[iter]    = (0.1 * D[iter]) / delta_f[iter]
            u_prime_S_L[iter]   = (0.030 * U0) / S_L[iter]


            # NOTE Objective 5a
            # Using corelations to calculate the turbulent flame speed.
            # Burner Rim Stabilized Flames

            # Calculating the turbulent Re
            # TODO to be verified because using a dummy value for intensity at the moment.
            Re_T = (0.030 * U0 * 0.1 * D[iter]) / (S_L[iter] * delta_f[iter])

            # NOTE Peters Model

            A = (Re_T / (u_prime_S_L[iter])) + 1
            E1 = 1 * (a_peters/2) * A
            ST_SL_P[iter] = -1 * E1 + np.sqrt((E1)**2 + (a_peters * A * u_prime_S_L[iter]) + a_peters + 1)

            # NOTE Zimont Model
            ST_SL_Z[iter] = A_zimont * (Pr_zimont)**(1/4) * (Re_T)**(1/4) * np.sqrt(u_prime_S_L[iter])

            global_iter = iter_start + iter
            if (global_iter+1)%4 == 0: # we print every 10th iteration.
                print(f"Iteration {global_iter+1:>4} | U0 = {U0:>7.2f} m/s | PHI = {PHI:>6.3f} | xM = {xM:>7.4f}")
                print(f"Iteration {global_iter+1:>4}| ST/SL P = {ST_SL_P[iter]:>7.4f} | ST/SL Z = {ST_SL_Z[iter]:>7.4f}")
                print("------------------")
            iter += 1
    
    return {
        'U0': U0_arr,
        'PHI': phi_arr,
        'xM': xM_arr,
        'T_ad': t_ad,
        'm_dot_air': m_dot_air,
        'D': D,
        'L_min': L_min,
        'L_actual': L,
        'S_L': S_L,
        'delta_f': delta_f,
        'lt_Lf': lt_delta_f,
        'ui_SL': u_prime_S_L,
        'ST_SL_Peters': ST_SL_P,
        'ST_SL_Zimont': ST_SL_Z
    }


def main(args):

    startTime = time()

    # Get command line flags
    save_data:bool  = args.data
    save_plot:bool  = args.plot

    # Need to define a loop that can iterate over several input parameters.
    # What can we define as the input variable?

    # NOTE: Design Variables
    U0_range    = np.arange(60, 241, 10)            # Bulk flow velocity
    PHI_range   = np.arange(0.45, 0.501, 0.01)      # Equivalence ratio
    xM_range    = np.arange(20, 26.01, 0.1)           # Non-dimensional location for turbulence intensity

    U0_range    = np.arange(60, 241, 10)            # Bulk flow velocity
    PHI_range   = np.arange(0.45, 0.48, 0.01)      # Equivalence ratio
    xM_range    = np.arange(20, 26.01, 2.0)           # Non-dimensional location for turbulence intensity


    design_variables:tuple = (U0_range, PHI_range, xM_range) # order matters.

    design_points = len_product(*design_variables)
    
    print(f"Total design points: {design_points}")
    print(f"Using {cpu_count() - 4} CPU cores for parallel processing")
    print("------------------")
    
    points_per_u0 = len(PHI_range) * len(xM_range)
    design_params = [(U0, i * points_per_u0) for i, U0 in enumerate(U0_range)]
    
    # Create partial function with fixed PHI_range and xM_range
    process_func = partial(process_design_point, PHI_range=PHI_range, xM_range=xM_range)
    
    with Pool(processes=cpu_count() - 4) as pool:
        results = pool.map(process_func, design_params)
    
    # Combine results from all processes
    combined_data = {
        'U0': np.concatenate([r['U0'] for r in results]),
        'PHI': np.concatenate([r['PHI'] for r in results]),
        'xM': np.concatenate([r['xM'] for r in results]),
        'T_ad': np.concatenate([r['T_ad'] for r in results]),
        'm_dot_air': np.concatenate([r['m_dot_air'] for r in results]),
        'D': np.concatenate([r['D'] for r in results]),
        'L_min': np.concatenate([r['L_min'] for r in results]),
        'L_actual': np.concatenate([r['L_actual'] for r in results]),
        'S_L': np.concatenate([r['S_L'] for r in results]),
        'delta_f': np.concatenate([r['delta_f'] for r in results]),
        'lt_Lf': np.concatenate([r['lt_Lf'] for r in results]),
        'ui_SL': np.concatenate([r['ui_SL'] for r in results]),
        'ST_SL_Peters': np.concatenate([r['ST_SL_Peters'] for r in results]),
        'ST_SL_Zimont': np.concatenate([r['ST_SL_Zimont'] for r in results])
    }
    
    df = pd.DataFrame(combined_data)

    print(df)
    df.to_csv('test_data.csv', index=False)
    print_stats(startTime, time())


if __name__ == "__main__":

    # Main arguments
    parser = argparse.ArgumentParser(prog='combustor', description='Cantera model for premixed flames.')

    # Settings
    parser.add_argument('-d', '--data',    action='store_true', help='Save all output data to csv.')
    parser.add_argument('-p', '--plot',    action='store_true', help='Save all plots.')

    arguments = parser.parse_args()
    main(arguments)