from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

c = []
start = 0
end = 20
for i in range(0,2,1):
    x = np.arange(start,end,0.01)
    c = np.append(c,x)
    start = end
    end+=end
print c
