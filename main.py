import scipy.constants
import numpy as np
import matplotlib.pyplot as plt

# Parameters

# radius of the port [m]
r = 0.0006
# volume of the pumping drop [m^3]
V = 0.5 * 1e-9 * 11
V_reservoir = 0.5 * 1e-9 * 2
# a = 3r^2
a = 3 * (r ** 2)
# b = 6V/pi
b = (6 * V) / scipy.constants.pi
b_reservoir = (6 * V_reservoir) / scipy.constants.pi
# length [m], height [m], width [m] of microchannel
L, H, W = 0.02, 0.00014, 0.001
# dynamic viscosity [Pa * s]
dynamic_viscosity = 8.9 * 1e-4

# reservoir drop density [kg * m^-3]
p = 998
# gravity [m * s^-2]
g = scipy.constants.g
# surface free energy of the liquid [J * m^-2]
gamma = 72.8 * 1e-3

# Height, Drop Radius, Resistance, Flow Rate
# height of a spherical cap of volume V [m]
h = (((1 / 6) * ((108 * b + 12 * ((12 * (a ** 3) + 81 * (b ** 2)) ** (1 / 2))) ** (1 / 3))) -
     ((2 * a) / ((108 * b + 12 * ((12 * (a ** 3) + 81 * (b ** 2)) ** (1 / 2))) ** (1 / 3))))
h_reservoir = (
            ((1 / 6) * ((108 * b_reservoir + 12 * ((12 * (a ** 3) + 81 * (b_reservoir ** 2)) ** (1 / 2))) ** (1 / 3))) -
            ((2 * a) / ((108 * b_reservoir + 12 * ((12 * (a ** 3) + 81 * (b_reservoir ** 2)) ** (1 / 2))) ** (1 / 3))))
# the drop radius [m]
R = (((3 * V) / scipy.constants.pi) + (h ** 3)) * (1 / (3 * (h ** 2)))
R_reservoir = (((3 * V_reservoir) / scipy.constants.pi) + (h_reservoir ** 3)) * (1 / (3 * (h_reservoir ** 2)))
# resistance of the microchannel [Pa * s * uL^-1]
Z = ((1 / (1 - 0.63 * (H / W))) * ((12 * dynamic_viscosity * L) / ((H ** 3) * W))) * 1e-9
# the volumetric flow rate [uL * s^-1]
dV_dt = (1 / Z) * ((p * g * h_reservoir) - ((2 * gamma) / R))

print(h, R, dV_dt)

# pressure difference
gamma = 72.8 * 1e-3
delta_p_pump = (2 * gamma) / R
gamma = 40 * 1e-3
delta_p_res = (2 * gamma) / R_reservoir
print(delta_p_pump, delta_p_res)
print(delta_p_pump - delta_p_res)
# 55 uL water drop as pumping vs 10 uL collagen drop as reservoir?
