#===========================================================
# Program written for  Combustion Dynamics - Fall 2025
#                        PROJECT
# Authors:
#   --  Paramvir Lobana      --
#   --  Lilou    Mandray     --
#===========================================================

import numpy as np
import pandas as pd
import cantera as ct
from time import time
import argparse
import matplotlib.pyplot as plt

# User imports
from inputs import *
import modules.thermodynamics as thermo

# Init gases
gas                    =       ct.Solution('gri30.yaml')

# Design variables:
EQR = 1.0

def main(
    args
    ):

    startTime = time()

    # 0. Init gas state
    gas_premixer = ct.Solution(os.path.join(MODELS, 'aramco3.yaml'))
    gas_premixer.TP = T_30, P_30
    gas_premixer.set_equivalence_ratio(phi=EQR, fuel="CH4:1", oxidizer="O2:0.21,N2:0.79")
    gas_premixer_density = gas_premixer.density_mass

    # 1. Obtain initial parameters: mass flow fuel, mass flow air, etc
    gas_combustion = ct.Solution(os.path.join(MODELS, 'aramco3.yaml'))
    gas_combustion.TP = T_30, P_30
    gas_combustion.set_equivalence_ratio(phi=EQR, fuel="CH4:1", oxidizer="O2:0.21,N2:0.79")
    t_ad = thermo.calc_AdiabaticTemperature(gas_combustion)
    mdot_air = func_mdot_air(EQR)

    # Add a check for t_ad. If between range, mark are good.
    if t_ad < 1750 or t_ad > 1850:
        print(f"❌ Adiabatic flame temperature {t_ad:.2f}K is out of the desired range (1750-1850K). Please adjust EQR.\n")
    else:
        print(f"✅ Adiabatic flame temperature is {t_ad:.2f}K, within the desired range.\n")


    # Define the mass flow rates again based on the number of premixers. Default set to 12, update in inputs.py file
    mdotAir = mdot_air/n_premixers
    mdotFuel = mdot_fuel/n_premixers
    mdotTotal = mdotAir + mdotFuel


    print("Mass flow stats (per premixer):")
    print(f"- Mass flow fuel:    {mdotFuel:10.3f} kg/s")
    print(f"- Mass flow air:     {mdotAir:10.3f} kg/s")
    print(f"- Total mass flow:   {mdotTotal:10.3f} kg/s")

    D = np.sqrt((mdotTotal * 4) / (gas_premixer_density * np.pi * initial_vel))
    print(D)


    T, t = thermo.computeIgnitionDelay(gas_premixer, T30=750, save=True)

    endTime = time()
    print("") 
    print("STATS:")
    print("-"*6)
    print(f"Program took {(endTime - startTime):10.03f}s to execute.")


def func_mdot_air(EQR:float):

    AFR = AFR_STOIC / EQR
    mdot_air = AFR * mdot_fuel

    return mdot_air



if __name__ == "__main__":

    # Main arguments
    parser = argparse.ArgumentParser(prog='combustor', description='Cantera model for cylindrical Bunsen-like jet flames.')
    parser.add_argument('-p', '--primary',      action='store_true', help='Runs initial calculations for the primary region.')


    # Settings
    parser.add_argument('-d', '--data',    action='store_true', help='Save all output data to csv.')

    arguments = parser.parse_args()
    main(arguments)
