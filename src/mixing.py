import os
import numpy as np

# User imports
import modules.utilities as ut
from inputs import *

x, y = ut.reader(os.path.join(DATA, "fig5-33.csv"))
f = ut.func_interpolate(x, y)
print(f(100))