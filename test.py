import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

print(f"Running Cantera version: {ct.__version__}")

# --- 1. Define Initial Conditions ---
# From project: T = 750 K, P = 24 bar 
T_in = 750.0  # K
P_in = 24.0 * ct.one_atm  # Pa

# We need a stoichiometric CH4/Air mixture 
# Stoichiometric reaction: CH4 + 2(O2 + 3.76N2) -> CO2 + 2H2O + 7.52N2
# Mole fractions: CH4: 1, O2: 2, N2: 2 * 3.76 = 7.52
mixture = 'CH4:1.0, O2:2.0, N2:7.52'

# Use a detailed mechanism for methane combustion
gas = ct.Solution('gri30.yaml')
gas.TPX = T_in, P_in, mixture

# --- 2. Set up the 0-D Reactor ---
# This is a 0-D analysis 
# We model this as a constant-pressure, adiabatic "batch" of gas
r = ct.IdealGasReactor(gas)
sim = ct.ReactorNet([r])

# --- 3. Run the Simulation ---
# We will store the time and temperature to find the delay
time_history = []
temp_history = []

# Set a total simulation time (e.g., 0.1 seconds).
# Adjust if ignition happens faster or slower.
t_end = 0.1  # s

current_time = 0.0
while current_time < t_end:
    time_history.append(current_time)
    temp_history.append(r.T)
    
    # Integrate forward in time
    current_time = sim.step()

# --- 4. Find the Autoignition Delay Time ---
# We define ignition delay (tau_ign) as the time at which the
# temperature gradient (dT/dt) is at its maximum.

# Convert lists to numpy arrays for easier processing
time_array = np.array(time_history)
temp_array = np.array(temp_history)

# Calculate the temperature gradient
dTdt = np.gradient(temp_array, time_array)

# Find the index of the maximum gradient
ign_index = np.argmax(dTdt)

# Find the corresponding time
tau_ign = time_array[ign_index]

print(f"\n--- Results ---")
print(f"Initial Temperature: {T_in} K")
print(f"Initial Pressure: {P_in / ct.one_atm:.2f} bar")
print(f"Mixture: Stoichiometric Methane-Air")
print(f"Autoignition Delay Time (tau_ign): {tau_ign * 1000:.4f} ms")

# --- 5. Plot the Results ---
plt.figure(figsize=(10, 6))
plt.plot(time_array * 1000, temp_array, lw=2)
plt.xlabel('Time (ms)')
plt.ylabel('Temperature (K)')
plt.title(f'Autoignition at {T_in} K, {P_in / ct.one_atm:.0f} bar')
plt.grid(True, which='both', linestyle='--', alpha=0.7)

# Mark the ignition delay time
plt.axvline(x=tau_ign * 1000, color='r', linestyle='--', label=f'$\\tau_{{ign}}$ = {tau_ign * 1000:.4f} ms')
plt.legend()
plt.show()