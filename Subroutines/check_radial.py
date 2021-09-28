def check_radial_position(r0,ry,rz):
    '''
    Check if molecule is within bounds.\n
    r0 = radius of the multipole\n
    ry,rz = off-axis position of molecule
    '''
    rint = (rz * rz) + (ry * ry)
    #if molecule is outside bounds return true
    if(rint >= (r0 * r0)): return 1
    return 0 