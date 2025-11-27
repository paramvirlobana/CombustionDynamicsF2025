import cantera as ct
from molmass import Formula
import numpy as np

def calc_AdiabaticTemperature(gas:ct.Solution, phi:float=None, **state) -> float:
    gas_temp = gas

    if phi is not None:
        fuel, oxidizer, T_init, P_init = state['T_init'], state['P_init']

        # Set gas state for given phi
        gas_temp.set_equivalence_ratio(phi, fuel=fuel, oxidizer=oxidizer)
        gas_temp.TP = T_init, P_init

    # Equilibrate the mixture adiabatically at constant pressure
    gas_temp.equilibrate('HP', solver='gibbs', max_steps=10000)
    tad = gas_temp.T

    return tad


def computeIgnitionDelay(gas:ct.Solution, T30:float, save:bool=False) -> tuple:

    """
    Function to compute and return the ignition delay.   
    """
    #oxidizer = {'O2':1.0, 'N2':3.76} # Moles
    #temp_gas = gas # Defined a temporary variable so the calculations do not update the original 'gas' object.
    #temp_gas.set_equivalence_ratio(phi=EQR, fuel=fuel, oxidizer=oxidizer, basis='mole')
    #temp_gas.TP = T30, P30

    # init the reactor
    r = ct.IdealGasConstPressureMoleReactor(gas, name='R1')
    sim = ct.ReactorNet([r])
    sim.derivative_settings = {"skip-third-bodies":True, "skip-falloff":True}
    sim.preconditioner = ct.AdaptivePreconditioner() # Manage the numerical convergence / iteration by "conditioning" the matrix.


    # Define the solution array # Append time:
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

            if (np.absolute(T_new - T_old) < 0.001 and T_old > (T30 + 800)):
                print(f"T30 = {T30:4.2f} K, Tf = {r.T:.2f} K and Ignition Delay = {t:.4f} seconds")
                break

            counter += 1      

    if save:
        states.to_pandas().to_csv("ignition_delay_output.csv")

    return r.T, t


def calc_AirFuelRatio(MW_AIR:float, MW_FUEL:float, stoic:bool=False) -> float:
    """
    Calculates the air to fuel ratio
    """
    if stoic:
        AFR_STOIC = (MW_AIR / MW_FUEL) * 4.76 * (3.0 / 4.0)
        print(f"Stoichiometric AFR: {AFR_STOIC:>.6f}")
        return AFR_STOIC
    else:
        AFR = MW_AIR / MW_FUEL
        return AFR

def calc_EquivalenceRatio(AFR_STOIC:float, AFR:float) -> float:
    """
    Calculates the equivalence ratio.
    """
    EQR = AFR_STOIC / AFR
    return EQR

if __name__ == '__main__':
    pass