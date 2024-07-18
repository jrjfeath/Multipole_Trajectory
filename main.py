import os 

from Subroutines.read_input import setup
from Subroutines.scan_voltage_trajectory import scan_voltage

cwd = os.path.dirname(__file__)

filename = '/Samples/ND3-Parameters-Double.yaml'

if __name__ == "__main__":
    d = setup(f'{cwd}/{filename}')
    scan_voltage(d)
