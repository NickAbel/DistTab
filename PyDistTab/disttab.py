#!/usr/bin/env python

from mpi4py import MPI
import numpy as np
import itertools
from math import floor

# Returns a baloney lookup table.
# Hardcoded for 2 CVs, 1 SV.
#TODO
# First few lines contain idea on how to split array btwn ranks
# which ought to be verified and generalized.
# Need to extend the tables so there are no 'gaps' between them
# (for example, need 0.49 but end of one table has 0.48 and
# beginning of other has 0.50; don't let this happen, perhaps
# just give the lower table 0.50)
def table_array(CVAR,CMEAN,state_vars,rank,nprocs):
    lower_rank_bound_cvar  = floor((rank)*CVAR.shape[0]/nprocs)
    upper_rank_bound_cvar  = floor((rank+1)*CVAR.shape[0]/nprocs)
    dim_rank_cvar          = upper_rank_bound_cvar - lower_rank_bound_cvar
    #If not the last rank, extend tables by 1 entry to fill gaps in
    if (rank < nprocs-1):
        dim_rank_cvar += 1
        upper_rank_bound_cvar +=1
    array = np.zeros([dim_rank_cvar, CMEAN.shape[0], 2+state_vars])
    k=0
    l=0
    for i in CVAR[lower_rank_bound_cvar:upper_rank_bound_cvar]:
        l=0
        for j in CMEAN:
            array[k,l,0] = i
            array[k,l,1] = j
            l=l+1
        k=k+1
    array[:,:,2] = np.random.rand(dim_rank_cvar, CMEAN.shape[0])
    print(rank,array.shape,array)
    return array

def create_window(array,location,rank,comm):
    if comm.Get_rank() == rank:
        win = MPI.Win.Create(array[location], disp_unit=1, comm=comm)
    else:
        win = MPI.Win.Create(None, comm=comm)
    win.Fence()
    return win

# Gets and returns a cloud of lookup table points for interpolation.
# Assumes that all CVs are in range on table, TODO return error if not so
def get_point_cloud(array,coords):
    print("Rank ",MPI.COMM_WORLD.Get_rank(),coords," fetching local point cloud")

    t = []

    a = 0
    b = array.shape[0]-1
    c = floor((a+b)/2)

    while (b-a > 1):
        midpoint_indices = np.zeros(array.ndim,dtype=np.int8)
        midpoint_indices[0] = c
        print(midpoint_indices,array.ndim,"####")
        if (array[tuple(midpoint_indices)] < coords[0]):
            a = c
        elif (array[tuple(midpoint_indices)] > coords[0]):
            b = c
        c = floor((a+b)/2)

    assert(b-a==1)

    t.append([a,b])

    #TODO "exact match" logic 
    #With 2 CVs need four points: [a,d], [b,d], [a,e], [b,e]
    for i in range(1, len(coords)):
        a = 0
        b = array.shape[i]-1
        c = floor((a+b)/2)
        while (b-a > 1):
            midpoint_indices = np.zeros(array.ndim,dtype=np.int8)
            midpoint_indices[0]  = t[i-1][0]
            midpoint_indices[i]  = c
            midpoint_indices[-1] = i
            if (array[tuple(midpoint_indices)] < coords[i]):
                a = c
            elif (array[tuple(midpoint_indices)] > coords[i]):
                b = c
            c = floor((a+b)/2)

        assert(b-a==1)
    
        t.append([a,b])

    print(array)

    for element in itertools.product(*t):
        print(element)

    return 0

# Checks if the lookup values needed are going to be in the table.
# Hardcoded for 2 CVs.
def coords_are_local(array,coords):
    cvar_is_local = array[0,0,0] <= coords[0] and array[-1,0,0] >= coords[0]
    cmean_is_local = array[0,0,1] <= coords[1] and array[0,-1,1] >= coords[1]
    return cvar_is_local and cmean_is_local

def main():
    # MPI Objects
    comm      = MPI.COMM_WORLD
    infonull  = MPI.INFO_NULL
    nprocs    = comm.Get_size()
    rank      = comm.Get_rank()

    # Lookup table concoction objects
    CVAR       = np.linspace(0,1,num=9)
    CMEAN      = np.linspace(0,1,num=5)
    ctrl_vars  = 2
    state_vars = 1

    # Window storage objects
    window_location = (1,1)
    windows = dict()

    # Create a baloney array
    array = table_array(CVAR,CMEAN,state_vars,rank,nprocs)

    # Some location to use for lookup in [0,1]x[0,1]
    req_coords = np.random.rand(ctrl_vars)

    if coords_are_local(array, req_coords):
        get_point_cloud(array, req_coords)


    windows[window_location] = create_window(array,window_location,0,comm)

    #print("Before get at rank ", rank ,": ", array[window_location])
    if rank != 0:
        windows[window_location].Get([array[window_location], MPI.DOUBLE], 0)
    windows[window_location].Fence()

    #print("After get at rank ", rank ,":  ", array[window_location])
    windows[window_location].Free()

if __name__=='__main__':
    main()
