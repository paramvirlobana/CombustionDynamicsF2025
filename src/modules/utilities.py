import numpy as np
from scipy.interpolate import interp1d

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
    import os
    cwd = os.getcwd()
    file = os.path.join(cwd, "04. Code - git repo", "src", "data", "fig5-33.csv")
    df1, df2 = reader(file)
    f = func_interpolate(df1, df2)
    print(f(100))