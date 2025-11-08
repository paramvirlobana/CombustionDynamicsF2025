import cantera as ct
from molmass import Formula

def calc_AdiabeticTemperature(phi, T_init:float=298.15, P_init:float=101325.0) -> float:

    # Define the cantera reaction mechanim
    gas = ct.Solution('gri30.yaml')

    # Fuel species changed to hydrogen
    fuel_species = 'CH4'

    # Set gas state for given phi
    gas.set_equivalence_ratio(phi, fuel_species, 'O2:0.21, N2:0.79')
    gas.TP = T_init, P_init

    # Equilibrate the mixture adiabatically at constant pressure
    gas.equilibrate('HP', solver='gibbs', max_steps=10000)

    tad = gas.T

    return tad


def calc_AirProperties():
    # Given
    T_FUEL: float = 500  # [K]
    P_FUEL: float = 1.5e+6 # [Pa] = 15 bar = 1.5 MPa

    T_AIR: float = 800
    P_AIR: float = 1.2e+6    # [Pa] = 12 bar = 1.2 MPa


    # Constants
    R_UNIV = 8.314462618  # J/(mol·K)

    # Molecular weights (kg/mol)
    MW_O2       = Formula('O2').mass * 1e-3
    MW_N2       = Formula('N2').mass * 1e-3
    MW_NH3      = Formula('NH3').mass * 1e-3

    MW_AIR  = 0.21 * MW_O2 * 1e-3 + 0.79 * MW_N2 * 1e-3
    MW_FUEL = MW_NH3 * 1e-3

    # Specific gas constants (J/kg.K)
    R_AIR  = R_UNIV / MW_AIR
    R_FUEL = R_UNIV / MW_FUEL


    # Densities [kg/m^3]
    DENSITY_AIR  = P_AIR  / (R_AIR  * T_AIR)
    DENSITY_FUEL = P_FUEL / (R_FUEL * T_FUEL)

    print(f"{'Quantity':<20} {'Value':>15} {'Unit':<10}")
    print("-" * 50)
    print(f"{'Air Density':<20} {DENSITY_AIR:>15.6f} {'kg/m3':<10}")
    print(f"{'Air Temperature':<20} {T_AIR:>15.1f} {'K':<10}")
    print(f"{'Air Pressure':<20} {P_AIR:>15.1f} {'Pa':<10}")
    print("-" * 50)
    print(f"{'Fuel Density':<20} {DENSITY_FUEL:>15.6f} {'kg/m3':<10}")
    print(f"{'Fuel Temperature':<20} {T_FUEL:>15.1f} {'K':<10}")
    print(f"{'Fuel Pressure':<20} {P_FUEL:>15.1f} {'Pa':<10}")
    print("-" * 50)

    return DENSITY_AIR, DENSITY_FUEL, MW_AIR, MW_FUEL


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
    calc_AirProperties()