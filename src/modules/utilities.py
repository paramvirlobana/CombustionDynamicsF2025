import numpy as np
from scipy.interpolate import interp1d

def reader(csv:str):
    data = np.loadtxt(csv, delimiter=',', unpack=True)

    return data


def func_interpolate(xData, yData, kind='linear'):
    return interp1d(xData, yData, kind=kind)


# Run for testing
if __name__ == "__main__":
    import os
    cwd = os.getcwd()
    file = os.path.join(cwd, "04. Code - git repo", "src", "data", "fig5-33.csv")
    df1, df2 = reader(file)
    f = func_interpolate(df1, df2)
    print(f(100))