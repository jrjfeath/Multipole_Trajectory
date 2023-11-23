def calculate_velocity(m,y,T=328):
    '''
    Calculate the velocity of your molecule/ion using it's heat capacity ratio
    and mass. The ratio for monoatomic species is 5/3 and for polyatomic species
    is 7/5.\n The temperature is assumed to be 328 K based off comparison to 
    experimental data.

    367 - 100% CH3I
    977 - 100% ND3
    584 - 100% Ar
    403 - 100% Kr
    '''
    T = 328
    kB = 1.3806503e-23
    mass = m * 1.6605390666e-27 #kg
    v = (((2 * kB) / mass) * (y / (y - 1)) * T) ** 0.5
    print(v)

if __name__ == "__main__":
    calculate_velocity(30,7/5,328)