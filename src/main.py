#===========================================================
# Program written for  Combustion Dynamics - Fall 2025
#                        PROJECT
# Authors:
#   --  Paramvir Lobana      --
#   --  Fouad    Al Laham    --
#   --  Lilou    Mandray     --
#===========================================================

import numpy as np
import pandas as pd
import cantera as ct
from time import time
import argparse
import matplotlib.pyplot as plt

# Input variables
# Variable                      Value           Units

# Operating conditions
P_30:float              =       24*ct.one_atm   # Pa
T_30:float              =       750             # K
power_out:int           =       34              # MW
eta_th:float            =       0.41            # %

# Pressure losses
P_loss_max:float        =       0.06            # %
P_loss_inlet:float      =       0.03            # %
P_loss_combustor        =       0.01            # %

# Mass flows
combustion_air:float    =       0.75            # %
cooling_air:float       =       0.11            # % 

# Init gases
gas                    =       ct.Solution('gri30.yaml')

def main(
    args
    ):
    startTime = time()

    # We define the fuel and oxidizer temperature and pressure.
    gas.TP = T_30, P_30

    time_array, temp_array, tau_ign = subroutine_velocity_analysis(gas, phi=1.0)

    print(f"\n--- Results ---")
    print(f"Initial Temperature: {T_30} K")
    print(f"Initial Pressure: {P_30 / ct.one_atm:.2f} bar")
    print(f"Mixture: Stoichiometric Methane-Air")
    print(f"Autoignition Delay Time (tau_ign): {tau_ign * 1000:.4f} ms")

    plt.figure(figsize=(10, 6))
    plt.plot(time_array * 1000, temp_array, lw=2)
    plt.xlabel('Time (ms)')
    plt.ylabel('Temperature (K)')
    plt.title(f'Autoignition at {T_30} K, {P_30 / ct.one_atm:.0f} bar')
    plt.grid(True, which='both', linestyle='--', alpha=0.7)

    plt.legend()
    plt.show()


    endTime = time()
    print("") 
    print("STATS:")
    print("-"*6)
    print(f"Program took {(endTime - startTime):10.03f}s to execute.")


def subroutine_velocity_analysis(
    gas:ct.Solution, phi:float
    ) -> tuple:

    """
    Analysis on the premixer to ensure autoignition does not occur.
    Need to prove this using 0D analysis, using phi = 1.
    Code similar to project 1.
    """
    # init storage arrays
    time_store:list = []
    temp_store:list = []

    gas.set_equivalence_ratio(phi=phi, fuel="CH4:1", oxidizer="O2:0.21,N2:0.79")

    # Create reactor
    r = ct.IdealGasReactor(gas)
    sim = ct.ReactorNet([r])

    t_max:float = 0.1   # max simulation time
    t_now:float = 0.0   # current time

    while t_now < t_max:
        time_store.append(t_now)
        temp_store.append(r.T)
        
        # Integrate forward in time
        t_now = sim.step()

    time_array = np.array(time_store)
    temp_array = np.array(temp_store)

    dTdt = np.gradient(temp_array, time_array) # temp gradient
    ign_index = np.argmax(dTdt)                # index max gradient
    tau_ign = time_array[ign_index]

    return time_array, temp_array, tau_ign



if __name__ == "__main__":

    # Main arguments
    parser = argparse.ArgumentParser(prog='combustor', description='Cantera model for cylindrical Bunsen-like jet flames.')
    parser.add_argument('-p', '--primary',      action='store_true', help='Runs initial calculations for the primary region.')


    # Settings
    parser.add_argument('-d', '--data',    action='store_true', help='Save all output data to csv.')

    arguments = parser.parse_args()
    main(arguments)
