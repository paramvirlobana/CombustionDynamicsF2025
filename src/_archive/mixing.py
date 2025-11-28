import os
import sys
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.interpolate import interp1d

# User imports
from inputs import *

def main():

    x, y = reader("/Users/raja/Library/CloudStorage/OneDrive-SharedLibraries-McGillUniversity/CombustionDynamicsF2025 - General/04. Code - git repo/src/data/csv/fig5-33.csv")
    f    = func_interpolate(x, y)
    
    x_M_range   = np.linspace(min(x), max(x), 101)
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



def reader(csv:str):
    data = np.loadtxt(csv, delimiter=',', unpack=True)
    return data

def func_interpolate(xData, yData, kind='linear'):
    return interp1d(xData, yData, kind=kind)

def len_product(*args:np.ndarray) -> int:
    """
    This function takes in arrays as inputs.

    It returns the total number of combinations possible
    by taking the product of the lengths of each array.
    """

    total_length = 1
    for array in args:
        total_length *= len(array)
    return total_length

def sum(*args:int) -> int:
    """
    This function takes in integers as inputs.
    It returns the sum of all the integers.
    """

    total = 0
    for number in args:
        total += number
    return total

def product(*iterables, repeat=1):
    """
    Modified the code from itertools library from python.
    This function returns the index values along with the actual value as a tuple.
    Returns two tuple:
        tuple(values) and tuple(indices)
    """
    if repeat < 0:
        raise ValueError('repeat argument cannot be negative')
    
    pools = [tuple(pool) for pool in iterables] * repeat
    result = [([], [])]  # (values, indices)

    for pool in pools:
        result = [
            (xv + [y], xi + [i])
            for xv, xi in result
            for i, y in enumerate(pool)
        ]

    for values, indices in result:
        yield tuple(values), tuple(indices)

def print_stats(start_time, end_time):
    print("") 
    print("STATS:")
    print("-"*6)
    print(f"Program took {(end_time - start_time):10.03f}s to execute.")



if __name__ == "__main__":
    main()