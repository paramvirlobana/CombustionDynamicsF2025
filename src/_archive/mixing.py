import os
import sys
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# User imports
import modules.utilities as ut
from inputs import *

def main():

    x, y = ut.reader(os.path.join(DATA, "fig5-33.csv"))
    f    = ut.func_interpolate(x, y)
    
    x_M_range   = np.linspace(min(x), max(x), 11)
    u_U_squared = np.zeros_like(x_M_range)
    I_range     = np.zeros_like(x_M_range)

    for i, x_M in enumerate(x_M_range):
        u_U_squared[i]  = f(x_M)
        I_range[i]      = np.sqrt(f(x_M))

    # NOTE 
    # Our target is to reach a value of I above 3%
    # This is because no mixing happens below the 3% threshold.

    mask = I_range >= 0.03
    x_M_filtered = x_M_range[mask]
    I_filtered = I_range[mask]

    if '-s' in sys.argv:

        sym=['circle-open','square-open','diamond-open','star-open']
        colours = ['red','blue','green','orange','black','pink','grey','purple','lightblue','magenta','brown','lightgreen']

        fig = go.Figure()

        fig.add_trace(
            go.Scatter(x=x, y=y, 
                    line=dict(color='blue', width=1),
                    name='AR22', 
                    showlegend=False),
        )

        fig.add_trace(
            go.Scatter(x=x_M_filtered, y=I_filtered, 
                    line=dict(color='blue', width=1),
                    name='AR22', 
                    showlegend=False),
        )

        fig.update_xaxes(
            range=[10, 500], autorange=True,
            type='log',
            showgrid=True,
            showline=True,
            mirror=True,
            automargin=True,
            zeroline=False,
            title="x/M"
        )

        fig.update_yaxes(
            type='log',  # Set log scale
            showgrid=True,
            showline=True,
            mirror=True,
            automargin=True,
            zeroline=False,
            title=r"$\frac{<u^2>}{U_0}$"
        )

        fig.update_layout(
            height=600, 
            width=600, 
            template='none',
            showlegend=True, 
            title='',
            legend=dict(x=0.5, y=-0.15, orientation="h", xanchor="center")
        )

        fig.show()

if __name__ == "__main__":
    main()