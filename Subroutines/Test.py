from math import sqrt
import numpy as np
import sys

c = np.array([1,2,3,4])
bad = [0,1]
print(c[...,~0])

sys.exit()

T = 328
kB = 1.3806503e-23
y = 7 / 5
mass = 142 * 1.6605390666e-27 #kg
v = (((2 * kB) / mass) * (y / (y - 1)) * T) ** 0.5
print(v)

J, K, M =  2.0, 1.0, 2.0 
epi, omega = 1.0, 0.5
A, B = 0.25021, 5.17340
#convert cm-1 to m-1
Am, Bm = A * 100, B * 100

h = 6.626068e-34 #m^2 kg / s
c = 299792458  #m / s

dipole = 0.1585 * 3.33654e-30 #Cm, C = Amp (A) * second (s)
V = 2000 # volts (V) = m2 kg s-3 A-1
n = 3
r0 = 7 / 1000 #meters

# omega has units of meters
w = (dipole * V) / ((h * c * Bm) + (h * c * (Am - Bm) * K * K))
lt = ((J * J) - (K * K)) * ((J * J) - (M * M))
lb = (J * J * J) * (2 * J - 1) * (2 * J + 1)
rt = ((J + 1) * (J + 1) - (K * K)) * ((J + 1) * (J + 1) - (M * M))
rb = ((J + 1) * (J + 1) * (J + 1)) * (2 * J + 1) * (2 * J + 3)
dW2dE = dipole * w * ((lt / lb) - (rt / rb))

if K == 0: K = epi * omega

dEdx = (n * (n - 1) * V) / (r0 * r0 * r0)
dW1dE = (dipole * M * K) / (J * (J + 1))
c = (dEdx * dW1dE) / mass
c2 = (dEdx * dW2dE) / (mass)
print(sqrt(c),c2)