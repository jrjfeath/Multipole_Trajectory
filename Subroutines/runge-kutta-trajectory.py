import numpy as np

J, K, M = 1, 1, 1
n = 3
mass = 30
dipole = 0.1585
E = 10300
r0 = 7.0 / 1000
r03 = r0 * r0 * r0

atomic_mass = 1.6605390666e-27 #(kg)
debye = 3.33654e-30 #A s m, electric dipole moment
dipole *= debye
mass *= atomic_mass

l = (-1 / mass)
m = (-dipole * M * K) / (J * (J + 1))
r = (n * (n-1) * E) / (mass * r03)
d2xdt2 = l * m * r

def f(t,p,v):
    #returns final velocity
    return d2xdt2 * np.array(p) * t + np.array(v)

C1 = f(0.001,(0.001,0.1),(200,300))
print(C1)