import numpy as np
import scipy.sparse as sp
from collections import defaultdict
import cPickle as pkl
#from networkx.generators.random_graphs import fast_gnp_random_graph
#from networkx.convert import to_scipy_sparse_matrix

#def makedata(time,size,prob):
#    graphs = [to_scipy_sparse_matrix(fast_gnp_random_graph(size,prob),dtype=int) for ii in range(time)]
#    return graphs

def group_edges_by_times(edges, maxtime, mintime=0):
    """ returns list of tupels: [(d,[(u,v),...]),...] """
    dct = defaultdict(list)
    for u, v, d in edges:
        dct[d].append((u, v))
    dct_s = dict.fromkeys(range(0, maxtime-mintime), [])
    for d in dct:
        dct_s[d - mintime] = dct[d]
    return dct_s.items()

def readdata(filename, symmetrisch=False):
    """ Reading data from a textfile formatted as follows:
        first_node | second_note | timestep.
        Returns a Dictionary (key=timestep value=adjacency-matrix)
    """

    # Auslesen der Daten als Array
    netarray_tmp = np.loadtxt(filename, dtype = int)
    netarray = netarray_tmp[:,0:3].copy()
    # Entpacktes Auslesen, um schneller Knoten und Zeiten zu bekommen
    if netarray_tmp.shape[1]>3:
        u, v, times, weight = np.loadtxt(filename, dtype = int, unpack=True)
    else:
        u, v, times = np.loadtxt(filename, dtype = int, unpack=True)
    # Anzahl der Knoten: groesster aller Indizes
    lower_bound = min(min(u), min(v))
    u -= lower_bound
    v -= lower_bound
    netarray[:,0:2] -= lower_bound
    number_of_nodes = max(max(u)+1, max(v)+1)
    #print number_of_nodes

    # das dict selbst
    # edges in bessere Form bringen
    edges = group_edges_by_times(netarray, max(times), min(times))

    # das dict selbst
    point = {}
    for d, es in edges:
        us = [u for u,v in es]
        vs = [v for u,v in es]
        bs = [True for i in range(len(es))]
        
        m = sp.csr_matrix((bs,(vs,us)),
                shape=(number_of_nodes, number_of_nodes), dtype=float)
        # !!!!!!!!! Annahme, dass Kante: u -> v und p(t+1) = Ap(t) !!!!!!!!
        if symmetrisch:
            point[d] = m+m.transpose() # !!!!!!!! symmetrisieren !!!!!!
        else:
            point[d] = m
    return point

def lonely_nodes(mlist):
    dim = mlist[0].shape[0]
    runtime = len(mlist)
    contact = np.zeros((dim,))
    for ii in range(dim):
        jj = 0
        while contact[ii]==0 and jj<runtime:
            contact[ii] += mlist[jj][ii,:].nnz + mlist[jj][:,ii].nnz
            jj += 1
    return contact > 0

def remove_nodes(mlist, mask):
    runtime = len(mlist)
    for  ii in range(runtime):
        mlist[ii] = mlist[ii][mask]
        mlist[ii] = mlist[ii].transpose()[mask].transpose()
    return mlist

def clean_data(mlist):
    ind = lonely_nodes(mlist)
    mlist = remove_nodes(mlist,ind)
    return mlist

def save_sparse_csr(filename,array):
    time    = len(array)
    data    = {}
    indices = {}
    indptr  = {}
    shape   = {}

    for ii in range(time):
        data[ii]    = array[ii].data
        indices[ii] = array[ii].indices
        indptr[ii]  = array[ii].indptr
        shape[ii]   = array[ii].shape

    newarray = [data, indices, indptr, shape]
    pkl.dump( newarray, open( filename, "wb" ) )

def load_sparse_csr(filename):

    array    = pkl.load( open (filename, "rb") )
    time     = len(array[0])
    newarray = {}

    for ii in range(time):
        data         = array[0][ii]
        indices      = array[1][ii]
        indptr       = array[2][ii]
        shape        = array[3][ii]
        matrix       = sp.csr_matrix( (data,indices,indptr), shape )
        newarray[ii] = matrix

    return newarray

if __name__ == "__main__": # (Als Testumgebung)
    A = readdata('sexual_contacts.dat')
    print A[2232]

