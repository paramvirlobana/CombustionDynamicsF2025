import cantera as ct
import numpy as np
import plotly.graph_objects as go

import modules.thermodynamics as thermo


# Input variables
# Variable                      Value           Units

# Operating conditions
P_30:float              =       24*ct.one_atm   # Pa
T_30:float              =       750             # K


t_ad_store = []

phi_range = np.linspace(0.3, 0.5, 101)
for phival in phi_range:
    t_ad = thermo.calc_AdiabeticTemperature(phi=phival, T_init=T_30, P_init=P_30)
    t_ad_store.append(t_ad)
    print(f"phi = {phival:.4f}, T_ad = {t_ad:.4f}")

t_ad_arr = np.array(t_ad_store)

fig = go.Figure()


fig.add_trace(
    go.Scatter(x=phi_range, y=t_ad_arr, line=dict(color='black', width=2)),
)

fig.update_layout(
    title_text="Mixture Fraction and Scalar Dissipation",
    title_x=0.5,
    legend=dict(x=0.02, y=0.98),
    xaxis=dict(title_text=r"$\phi$"),
    yaxis=dict(title_text=r"$t_{ad}$", color='blue'),
    template="plotly_white"
)

fig.show()
