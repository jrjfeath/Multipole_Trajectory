def check_collision(vel,vy,vz,ry,rz,lcol,D=0.001):
    '''
    Calculate if molecule hits collision area.
    '''
    #Calculate position at distance d from end of hexapole & increment image intensity
    tcol = lcol / vel
    ryf = (vy * tcol) + ry
    rzf = (vz * tcol) + rz   

    # is the trajectory in the interaction region ( 1mm x 2mm x 1mm ) --> 2mm in beam velocity direction
    if (-D <= rzf <= D) and (-D <= ryf <= D):
        return 1
    return 0