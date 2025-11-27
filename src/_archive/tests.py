import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# User imports
import modules.utilities as ut
from inputs import *

D = 0.05 # [m]

def main():
    
    data = os.path.join(DATA, "fig5-35.csv")
    x, y = process_variance_data(data)
    
def process_variance_data(path):

    x, y = ut.reader(path)

    y_max = max(y)
    max_diff = y_max - 0.25

    y_new = y - max_diff

    return x, y_new



    gas_air = ct.Solution(os.path.join(MODELS, 'ar22.yaml'))
    gas_air.TPX = T_30, P_30, "O2:0.21,N2:0.79"
    gas_air = gas_air.density_mass

    gas_premixer = ct.Solution(os.path.join(MODELS, 'ar22.yaml'))
    gas_premixer.TP = T_30, P_30
    gas_premixer.set_equivalence_ratio(phi=EQR, fuel="CH4:1", oxidizer="O2:0.21,N2:0.79")
    gas_premixer_density = gas_premixer.density_mass

    # 1. Obtain initial parameters: mass flow fuel, mass flow air, etc
    gas_combustion = ct.Solution(os.path.join(MODELS, 'ar22.yaml'))
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
    mdotAir = mdot_air/n_premixers      # turbulent length scale 0.380 m or mm?
    mdotFuel = mdot_fuel/n_premixers    
    mdotTotal = mdotAir + mdotFuel

    print("Mass flow stats (per premixer):")
    print(f"- Mass flow fuel:    {mdotFuel:10.4f} kg/s")
    print(f"- Mass flow air:     {mdotAir:10.4f} kg/s")
    print(f"- Total mass flow:   {mdotTotal:10.4f} kg/s")

    D = np.sqrt((mdotAir * 4) / (gas_premixer_density * np.pi * initial_vel))
    print(f"- Diameter:          {D:10.4f} m\n\n")


    T, t = thermo.computeIgnitionDelay(gas_premixer, T30=750, save=True)

if __name__ == "__main__":
    main()



def texty():
        # Create log-spaced axis values
    x = np.logspace(-1, 2, 200)     # l_t / δ
    y = np.logspace(-2, 2, 200)     # u'/S_L^0

    # Regime boundaries (schematic power laws)
    Re_t = 1 / x                      # Re_t = 1 line
    Ka_1 = x**0.7                     # Ka = 1 (arbitrary slope for schematic)
    Ka_100 = 100 * x**0.7             # Ka = 100
    Da_1 = x**0.3                     # Da = 1

    fig = go.Figure()

    # --- Shaded regions ---------------------------------------------------------

    # Laminar region (bottom-left)
    fig.add_trace(go.Scatter(
        x=[1e-1, 1, 1, 1e-1],
        y=[1e-2, 1e-2, 1, 1],
        line=dict(color="rgba(0,0,0,0)"),
        name="Laminar"
    ))

    # Wrinkled flamelets (mid-right)
    fig.add_trace(go.Scatter(
        x=[1, 1e2, 1e2, 1],
        y=[1e-2, 1e-2, 1, 1],
        line=dict(color="rgba(0,0,0,0)"),
        name="Wrinkled flamelets"
    ))

    # Thickened / thickened-wrinkled (upper region)
    fig.add_trace(go.Scatter(
        x=[1e-1, 1e2, 1e2, 1e-1],
        y=[1, 1e2, 1e2, 1],
        line=dict(color="rgba(0,0,0,0)"),
        name="Thickened / Thickened-wrinkled"
    ))


    # --- Boundary lines ---------------------------------------------------------

    fig.add_trace(go.Scatter(x=x, y=Re_t, mode="lines",
        line=dict(color="black", width=2),
        name="Re_t = 1"))

    fig.add_trace(go.Scatter(x=x, y=Ka_1, mode="lines",
        line=dict(color="black", width=2, dash="dash"),
        name="Ka = 1"))

    fig.add_trace(go.Scatter(x=x, y=Ka_100, mode="lines",
        line=dict(color="black", width=2, dash="dot"),
        name="Ka = 100"))

    fig.add_trace(go.Scatter(x=x, y=Da_1, mode="lines",
        line=dict(color="black", width=2, dash="longdash"),
        name="Da = 1"))


    # --- Labels ---------------------------------------------------------

    fig.add_annotation(x=0.3, y=0.3, text="Laminar combustion",
                    showarrow=False, font=dict(size=14))

    fig.add_annotation(x=10, y=0.3, text="Wrinkled flamelets",
                    showarrow=False, font=dict(size=14))

    fig.add_annotation(x=5, y=20, text="Thickened-wrinkled flame",
                    showarrow=False, font=dict(size=14))

    fig.add_annotation(x=0.3, y=10, text="Thickened flame",
                    showarrow=False, font=dict(size=14))


    # --- Axes and layout ---------------------------------------------------------

    fig.update_layout(
        xaxis=dict(
            type="log",
            title="Integral length scale / flame thickness (l_t / δ)",
            range=[-1, 2]
        ),
        yaxis=dict(
            type="log",
            title="RMS velocity / flame speed (u' / S_L^0)",
            range=[-2, 2]
        ),
        width=1000,
        height=1000*1.294,
        showlegend=False
    )

    fig.show()
