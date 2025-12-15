import cantera as ct
import numpy as np
import plotly.graph_objects as go

# User imports
from inputs import *

# Input variables
# Variable                      Value           Units

# Operating conditions
P_30:float              =       24*ct.one_atm   # Pa
T_30:float              =       750             # K

# Fuel

def calc_AdiabaticTemperature(gas:ct.Solution, phi:float=None, **state) -> float:
    gas_temp = gas
    T_init, P_init = state['T_init'], state['P_init']

    # Set gas state for given phi
    gas_temp.set_equivalence_ratio(phi, fuel={"CH4":1}, oxidizer={"O2":0.21 ,"N2":0.79})
    gas_temp.TP = T_init, P_init

    # Equilibrate the mixture adiabatically at constant pressure
    gas_temp.equilibrate('HP', solver='gibbs', max_steps=10000)
    tad = gas_temp.T

    return tad



t_ad_store = []
mdot_air_store = []

phi_range = np.linspace(0.2, 1.1, 101)

for phival in phi_range:
    gas = ct.Solution('ar22.yaml')
    
    t_ad = calc_AdiabaticTemperature(gas, phi=phival, T_init=T_30, P_init=P_30)

    AFR         = AFR_STOIC / phival
    mdot_air    = AFR * mdot_fuel_total

    t_ad_store.append(t_ad)
    mdot_air_store.append(mdot_air)


    print(f"phi = {phival:.4f}, T_ad = {t_ad:.4f}, mdot air = {mdot_air:4f}")

t_ad_arr = np.array(t_ad_store)
mdot_air_arr = np.array(mdot_air_store)

import plotly.graph_objects as go

fig = go.Figure(data=[
    go.Scatter(x=phi_range, y=t_ad_arr, mode='lines', name='My Line')
])

# Add shaded region between x = 0.46 and 0.48
fig.add_shape(
    type="rect",
    x0=0.45, x1=0.5,
    y0=min(t_ad_arr),  # full vertical span
    y1=max(t_ad_arr),
    fillcolor="blue",
    opacity=0.5,
    line_width=0
)

fig.update_layout(template='simple_white')
fig.update_yaxes(showgrid=True, gridwidth=1)
fig.update_xaxes(showgrid=True, gridwidth=1)

fig.update_layout(
    height=600,
    width=500,
    xaxis_title="Equivalence ratio [-]",
    yaxis_title=r"$T_{ad}$",
    xaxis=dict(mirror=True, ticks='outside', showline=True),
    yaxis=dict(mirror=True, ticks='outside', showline=True),
    legend=dict(x=0.5, y=-0.15, orientation="h", xanchor="center"),
    margin=dict(l=50, r=20, t=50, b=50),
    font=dict(size=14)
)

fig.write_image("t_ad_vs_eqr.pdf")
fig.show()
