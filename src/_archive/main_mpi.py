#===========================================================
# Program written for  Combustion Dynamics - Fall 2025
#                        PROJECT
# Authors:
#   --   PL   --
#   --   FL   --
#   --   LM   --
#===========================================================

import sys
import numpy as np
import pandas as pd
import cantera as ct
from time import time
import argparse
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count, Manager
from functools import partial
from tqdm import tqdm
import threading

# User imports
from inputs import *
from modules.utilities import *

# set mechanism
mech:str    = 'ar22.yaml'

# Set fuel and oxidizer
fuel:dict       = {"CH4":1}
oxidizer:dict   = {"O2":0.21 ,"N2":0.79}

gas_air = ct.Solution("air.yaml")
gas_air.TP = T_30, P_30
density_air = gas_air.density_mass
viscosity_air = gas_air.viscosity

def main(args):
    startTime = time()

    compute_design_space(args)

    print_stats(startTime, time())


def process_design_point(design_params, PHI_range, M_range, progress_queue=None):
    """
    Process a single design point for a given U0 value.
    This function handles all PHI and xM combinations for one U0.
    """
    U0, iter_start, process_id = design_params
    
    # Calculate number of points for this U0
    n_points = len(PHI_range) * len(M_range)
    
    # Pre-allocate arrays for this U0 using a dictionary for cleaner code
    results = {
        'U0': np.zeros(n_points),
        'PHI': np.zeros(n_points),
        'M': np.zeros(n_points),
        'T_ad': np.zeros(n_points),
        'm_dot_air': np.zeros(n_points),
        'D': np.zeros(n_points),
        'S_L': np.zeros(n_points),
        'delta_f': np.zeros(n_points),
        'lt_Lf': np.zeros(n_points),
        'ui_SL': np.zeros(n_points),
        'L_premixer': np.zeros(n_points),
        'x_injector': np.zeros(n_points),
        'dP_grid': np.zeros(n_points),
        'dP_fric': np.zeros(n_points),
        'dP_total': np.zeros(n_points),
        'dP_percent': np.zeros(n_points),
        'ST_SL_Peters': np.zeros(n_points),
        'ST_SL_Zimont': np.zeros(n_points)
    }
    
    # design iterator
    for (PHI, M), (i_phi, i_m) in product(PHI_range, M_range):
        local_iter = i_phi * len(M_range) + i_m
        # To store the iteration value.
        results['U0'][local_iter] = U0
        results['PHI'][local_iter] = PHI
        results['M'][local_iter] = M

        # Initialize gas object for equilibrium calculations
        gas_A = ct.Solution(mech)
        gas_A.TP = T_30, P_30
        gas_A.set_equivalence_ratio(phi=PHI, fuel=fuel, oxidizer=oxidizer)

        # Calculate t_ad for each design.
        gas_A.equilibrate('HP', solver='gibbs', max_steps=10000)
        results['T_ad'][local_iter] = gas_A.T

        # Using equivalence ratio to get mass flow of air. 
        # NOTE: This is the mass flow of air in each premix chamber.
        results['m_dot_air'][local_iter] = (AFR_STOIC * mdot_fuel_total) / (PHI * n_premixers)

        # Using conservation of mass, we can calculate the required diameter.
        D_val = np.sqrt((results['m_dot_air'][local_iter] * 4) / (density_air * np.pi * U0))
        results['D'][local_iter] = D_val

        # Calculating the length of the premixer where the intensity falls < 3%.
        Irange = np.sqrt(y1) * 100
        I_to_x = interp1d(Irange, x1, kind='linear', fill_value='extrapolate')
        x_at_val = I_to_x(2.95)     # We target for 2.95% turbulent intensity.
        cutoff_xM = x1[np.argmin(np.abs(x1 - x_at_val))]
        results['L_premixer'][local_iter] = np.round(cutoff_xM * M, 6)

        # TODO
        # Injector location is where variance is highest.
        results['x_injector'][local_iter] = D_val / (2 * 0.12) # assuming 85% mixing efficiency before injector.

        # Calculate the mixing timescale (This is objective 3).
        eta = (D_val / 2) / results['L_premixer'][local_iter]


        # New gas object for pressure loss calculations.
        gas_B = ct.Solution(mech)
        gas_B.TP = T_30, P_30
        gas_B.set_equivalence_ratio(phi=PHI, fuel=fuel, oxidizer=oxidizer)

        # Calculate pressure losses using the pressure_loss function
        pressure_results = pressureLoss(density_air, viscosity_air, U0, D_val, results['L_premixer'][local_iter])
        
        results['dP_grid'][local_iter] = pressure_results['dP_grid']
        results['dP_fric'][local_iter] = pressure_results['dP_fric']
        results['dP_total'][local_iter] = pressure_results['dP_total']
        results['dP_percent'][local_iter] = pressure_results['dP_percent']

        # Compute the laminar flame speed.
        # We initialize a new gas object because the previous object is
        # already solved and in equilibrium state.

        gas_C = ct.Solution(mech)
        gas_C.TP = T_30, P_30
        gas_C.set_equivalence_ratio(phi=PHI, fuel=fuel, oxidizer=oxidizer)

        flame = ct.FreeFlame(gas_C, width=D_val)
        flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)
        flame.solve(loglevel=0, auto=True)
        results['S_L'][local_iter] = flame.velocity[0]

        # Compute the flame thickness
        x = flame.grid
        T = flame.T
        results['delta_f'][local_iter] = (T[-1] - T[0]) / (max(np.gradient(T, x)))

        # Plot results on the Borghi diagram.
        results['lt_Lf'][local_iter] = (0.1 * D_val) / results['delta_f'][local_iter]
        results['ui_SL'][local_iter] = (0.0295 * U0) / results['S_L'][local_iter]

        # NOTE Objective 5a
        # Using corelations to calculate the turbulent flame speed.
        # Burner Rim Stabilized Flames

        # Calculating the turbulent Re
        # TODO to be verified because using a dummy value for intensity at the moment.
        Re_T = (0.0295 * U0 * 0.1 * D_val) / (results['S_L'][local_iter] * results['delta_f'][local_iter])

        # NOTE Peters Model

        A = (Re_T / (results['ui_SL'][local_iter])) + 1
        E1 = 1 * (a_peters/2) * A
        results['ST_SL_Peters'][local_iter] = -1 * E1 + np.sqrt((E1)**2 + (a_peters * A * results['ui_SL'][local_iter]) + a_peters + 1)

        # NOTE Zimont Model
        results['ST_SL_Zimont'][local_iter] = A_zimont * (Pr_zimont)**(1/4) * (Re_T)**(1/4) * np.sqrt(results['ui_SL'][local_iter])

        global_iter = iter_start + local_iter
        
        # Update progress for this process
        if progress_queue is not None:
            progress_queue.put((process_id, 1))
            
    return results

def compute_design_space(args):

    # Get command line flags
    save_data:bool  = args.data
    save_plot:bool  = args.plot

    # Need to define a loop that can iterate over several input parameters.
    # What can we define as the input variable?

    # NOTE: Design Variables
    U0_range    = np.arange(60, 101, 10)            # Bulk flow velocity
    PHI_range   = np.arange(0.46, 0.4801, 0.01)      # Equivalence ratio
    M_range     = np.arange(0.003, 0.01, 0.001)     # 

    design_variables:tuple = (U0_range, PHI_range, M_range) # order matters.

    design_points = len_product(*design_variables)
    
    print(f"Total design points: {design_points}")
    num_processes = min(cpu_count(), len(U0_range))
    print(f"Using {num_processes} CPU cores.")
    print("------------------")
    
    # Prepare parameters for each U0 value
    points_per_u0 = len(PHI_range) * len(M_range)
    design_params = [(U0, i * points_per_u0, i) for i, U0 in enumerate(U0_range)]
    
    manager = Manager()
    progress_queue = manager.Queue()
    
    progress_bars = {}
    for i, U0 in enumerate(U0_range):
        progress_bars[i] = tqdm(
            total=points_per_u0,
            desc=f"U0={U0:>6.1f} m/s",
            position=i,
            leave=True,
            ncols=100,
            bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]'
        )
    
    stop_event = threading.Event()
    monitor_thread = threading.Thread(
        target=progress_monitor,
        args=(progress_queue, progress_bars, stop_event)
    )
    monitor_thread.start()
    
    process_func = partial(process_design_point, PHI_range=PHI_range, M_range=M_range, progress_queue=progress_queue)
    
    with Pool(processes=num_processes) as pool:
        results = pool.map(process_func, design_params)
    
    stop_event.set()
    monitor_thread.join()
    
    for pbar in progress_bars.values():
        pbar.close()
    
    print("\nAll processes completed!")
    print("------------------")
    
    combined_data = {
        key: np.concatenate([r[key] for r in results])
        for key in results[0].keys()
    }
    
    df = pd.DataFrame({
        'U0': combined_data['U0'],
        'PHI': combined_data['PHI'],
        'M': combined_data['M'],
        'T_ad': combined_data['T_ad'],
        'm_dot_air': combined_data['m_dot_air'],
        'D': combined_data['D'],
        'L_premixer': combined_data['L_premixer'],
        'x_injector': combined_data['x_injector'],
        'S_L': combined_data['S_L'],
        'delta_f': combined_data['delta_f'],
        'lt_Lf': combined_data['lt_Lf'],
        'ui_SL': combined_data['ui_SL'],
        'ST_SL_Peters': combined_data['ST_SL_Peters'],
        'ST_SL_Zimont': combined_data['ST_SL_Zimont'],
        'dP_grid': combined_data['dP_grid'],
        'dP_fric': combined_data['dP_fric'],
        'dP_total': combined_data['dP_total'],
        'dP_percent': combined_data['dP_percent']
    })

    print(df)
    if save_data:
        df.to_csv('test_data.csv', index=False)


def progress_monitor(progress_queue, progress_bars, stop_event):
    """
    Monitor progress from multiple processes and update progress bars.
    
    Parameters:
    -----------
    progress_queue : multiprocessing.Queue
        Queue receiving progress updates from worker processes
    progress_bars : dict
        Dictionary of tqdm progress bars for each process
    stop_event : threading.Event
        Event to signal when monitoring should stop
    """
    while not stop_event.is_set():
        try:
            # Get progress update with timeout
            process_id, increment = progress_queue.get(timeout=0.1)
            if process_id in progress_bars:
                progress_bars[process_id].update(increment)
        except:
            # Timeout or empty queue
            continue
    
    # Process any remaining items in queue
    while not progress_queue.empty():
        try:
            process_id, increment = progress_queue.get_nowait()
            if process_id in progress_bars:
                progress_bars[process_id].update(increment)
        except:
            break

def pressureLoss(rho, mu, U0, D_val, L_val):
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
        'dP_percent': dP_percent,
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