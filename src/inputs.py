import os
import cantera as ct

# Directories
DIR     = os.path.dirname(os.path.realpath(__file__))
DATA    = os.path.join(DIR, "data", "csv")
MODELS  = os.path.join(DIR, "data", "models")
FIGS    = os.path.join(DIR, "figs")


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

# Fuel
AFR_STOIC:float         =       17.11           # For methane - air combustion
mdot_fuel:float         =       1.66            # kg/s (calculated)

# Premixers
n_premixers:int         =       6.0             # Total number of premixers
initial_vel:float       =       180             # m/s