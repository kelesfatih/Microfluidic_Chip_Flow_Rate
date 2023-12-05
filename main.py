# A passive pumping method for microfluidic devices, https://doi.org/10.1039/B204381E
import scipy.constants
import numpy as np
import matplotlib.pyplot as plt

# Parameters

# radius of the port [m]
r = 0.0006
# volume of the pumping drop [m^3]
V = 0.5 * 1e-9
# a = 3r^2
a = 3 * (r ** 2)
# b = 6V/pi
b = (6 * V) / scipy.constants.pi

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
# the drop radius [m]
R = (((3 * V) / scipy.constants.pi) + (h ** 3)) * (1 / (3 * (h ** 2)))
# resistance of the microchannel [Pa * s * uL^-1]
Z = ((1 / (1 - 0.63 * (H / W))) * ((12 * dynamic_viscosity * L) / ((H ** 3) * W))) * 1e-9
# the volumetric flow rate [uL * s^-1]
dV_dt = (1 / Z) * ((p * g * L) - ((2 * gamma) / R))

print(h, R, Z, dV_dt)


def dV_dt(t, V):
    # Define the equation for dV_dt here
    # For example:
    # dV_dt = (1 / Z) * ((p * g * L) - (2 * gamma) / R)
    return (1 / Z) * ((p * g * L) - (2 * gamma) / V)


def rungeKutta(t0, V0, t, h):
    n = int((t - t0) / h)
    V = V0
    for i in range(1, n + 1):
        k1 = h * dV_dt(t0, V)
        k2 = h * dV_dt(t0 + 0.5 * h, V + 0.5 * k1)
        k3 = h * dV_dt(t0 + 0.5 * h, V + 0.5 * k2)
        k4 = h * dV_dt(t0 + h, V + k3)
        V = V + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        t0 = t0 + h
    return V


# Example usage
t0 = 0
V0 = 0.5  # Replace with the initial value of V
t = 0.05  # Replace with the final time
h = 0.01  # Replace with the step size
result = rungeKutta(t0, V0, t, h)
print('The value of V at t is:', result)
