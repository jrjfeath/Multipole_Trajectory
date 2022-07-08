import os 

from Subroutines.read_input import setup
from Subroutines.scan_voltage_opt import scan_trajectory
from Subroutines.scan_voltage_opt import scan_voltage

cwd = os.getcwd()

filename = '/Samples/NO-Parameters.yaml'

if __name__ == "__main__":
    d = setup(f'{cwd}/{filename}')
    if d['calcv'] == True: scan_voltage(d)
    if d['calct'] == True: scan_trajectory(d)
