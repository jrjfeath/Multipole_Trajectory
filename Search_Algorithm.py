import numpy as np

def randomize_distances(d):
    print (np.random.dirichlet(np.ones(d['counts']),size=1))
    '''md = d['lmax'] #max distance
    for key in d['multipole']:
        d['multipole'][key]['lpole']'''
    return d

def setup_search(d):
    
    return d