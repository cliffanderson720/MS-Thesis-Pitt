#!/usr/bin/env python
import velocity_functions as vf
import matplotlib.pyplot as plt
import numpy as np
from driver import driver
H = driver(H0 = 0.05,H_loc=0.5)
H.plot()
H.plot('gif')
