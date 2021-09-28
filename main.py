from Subroutines.read_input import setup
from Subroutines.scan_trajectory import scan_trajectory
from Subroutines.scan_voltage import scan_voltage

filename = '/home/josh/Dropbox/Python/Brouard/Multipole_Design_Estimator/Code/Samples/CH3I-Double-Parameters.yaml'

if __name__ == "__main__":
    d = setup(filename)
    if d['calcv'] == True: scan_voltage(d)
    if d['calct'] == True: scan_trajectory(d)
