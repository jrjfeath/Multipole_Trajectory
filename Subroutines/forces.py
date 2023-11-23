from math import sqrt

#Speed of light
c = 299792458  #m * s_1
#Plank's constant
h = 6.626068e-34 #m^2 * kg * s_1
hbar = 1.05457148e-34 #m^2 * kg * s_1 = (h/2pi)
hc = h * c #m^3 * kg * s_2

#References:
#1. Anderson, R. W. (1997). Tracks of Symmetric Top Molecules in Hexapole Electric Fields. 
#   The Journal of Physical Chemistry A, 101(41), 7664–7673. https://doi.org/10.1021/jp971313s
#2. Zare, R. N. (1988). Angular Momentum: Understanding Spatial Aspects in Chemistry and 
#   Physics - Baker Lecture Series. John Wiley & Sons Inc.

def first_order(d,n,r0):
    #Pull variables from dict
    j, k, m, dip, mass = d['j'], d['k'], d['m'], d['dipole'], d['mass']
    
    #Check if molecule is linear
    if k == 0: 
        epi, omega = d['epi'], d['omega']
        k = epi * omega
    
    w_e = {}
    #Equation from Ref. 2 pg 124
    #interaction energy (w) in units of cm^-1 (cm_1)
    #E = (m^2 * kg * s_3 * A_1) / cm_1, dip = C * m => A * s * m
    for E in d['V']:
        #Does the molecule have inversion or doublet splitting?
        try: 
            #w_inv in units of cm_1, multiply by hc to get m^3 * kg * s_2 * cm_1
            w_inv = d['w_inv'] * hc * 100 # In Joules
            l = (w_inv / 2) ** 2
            #r in units of m^3 * kg * s_2 * cm_1
            r = ((dip * E * 100) * ((m * k) / (j * (j + 1)))) ** 2
            w = (l + r)**0.5
            w = (-w_inv / 2) + w
            w_e[E] = w / hc / 100 # convert J to cm-1

        #If no splitting use generic first order stark equation
        except KeyError:
            #w in units of m^3 * kg * s_2 * cm_1
            w = abs(dip * E * m * k) / (j * (j+1))
            #divide by hc to obtain units in cm_1
            w_e[E] = w / hc

    return w_e

def second_order(d,n,r0):
    #Pull variables from dict
    j, k, m, dip, b, mass = d['j'], d['k'], d['m'], d['dipole'], d['b'], d['mass']

    #Determine second order perturbation theory interaction energy (cm^-1)
    #Equations are split into chunks to make it easier to follow
    #Check if molecule is linear
    if j < 1:
        w2 = 0

    elif k == 0: 
        top = (j * (j + 1)) - (3 * m * m)
        bottom = (j * (j + 1)) * (2 * j - 1) * (2 * j + 3)
        w2 = top / bottom

    else:
        lt = ((j * j) - (k * k)) * ((j * j) - (m * m))
        lb = (j * j * j) * (2 * j - 1) * (2 * j + 1)
        rt = ((j + 1) * (j + 1) - (k * k)) * ((j + 1) * (j + 1) - (m * m))
        rb = ((j + 1) * (j + 1) * (j + 1)) * (2 * j + 1) * (2 * j + 3)
        w2 = (lt / lb) - (rt / rb)

    dip2 = dip * dip
    thcb = 2 * h * c * b
    w_e = {}
    #Equation from Ref. 2 pg 124
    #interaction energy (w) in units of cm^-1 (cm_1)
    #E = (m^2 * kg * s_3 * A_1) / cm_1, dip = C * m => A * s * m
    for E in d['V']:
        #force is in units of m^3 * kg * s_2 * cm_1 (J * cm-1)
        w = abs(((dip2 * E * E) / thcb) * w2)
        #divide by hc to obtain units in cm_1
        w_e[E] = w / hc
    return w_e

def calculate_a(d,n,r0):
    '''
    Calculate the force applied on an ion done by the multipole.
    '''
    q, mass, w = d['charge'], d['mass'], d['Vfrequency']

    #make empty lists for saving results
    sqc, invc = [], []
    for V in d['V']:
        a = (2 * 2**(n-1) * V * q) / (mass * w * w * r0**n)
        sqc.append(sqrt(a))
        invc.append(1/sqc[-1])
    return sqc, invc