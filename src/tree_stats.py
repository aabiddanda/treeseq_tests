"""
    Computing measures of similarity on sets of tree sequences
"""
import tskit
import numpy as np 
from tqdm import tqdm


def ntrees_consist(ts1,ts2):
    """
        Given a set of two tree sequences 
        check whether they have the same number of trees
    """
    return(ts1.num_trees == ts2.num_trees)


def breakpoint_consist(ts1,ts2):
    """
        Given two tree sequences make sure that 
        the breakpoints are the same
    """
    break_pts1 = ts1.breakpoints(as_array=True)
    break_pts2 = ts2.breakpoints(as_array=True)
    eq = np.all(break_pts1 == break_pts2)
    return(eq)


def pairwise_tmrcas(ts, id1, id2):

    # Find the node ids for the individual IDs
    node_ids = np.array([i.nodes for i in ts.individuals()])
    indiv_ids = np.array([eval(i.metadata)['name'] for i in ts.individuals()])
    idx1 = np.where(indiv_ids == id1)[0]
    assert(idx1.size == 1)
    idx2 = np.where(indiv_ids == id2)[0]
    assert(idx2.size == 1)
    node_ids_1 = node_ids[idx1].flatten()
    node_ids_2 = node_ids[idx2].flatten()
    nnodes1 = node_ids_1.size
    nnodes2 = node_ids_2.size
    # Now we iterate through the trees
    tmrcas = np.zeros(shape=(ts.num_trees, nnodes1, nnodes2))
    progress = tqdm(total=ts.num_trees)
    k = 0
    for t in ts.trees():
        for i in range(nnodes1):
            for j in range(nnodes2):
                tmrcas[k,i,j] = t.tmrca(node_ids_1[i], node_ids_2[j])
        k += 1
        progress.update()
    progress.close()
    return(tmrcas)

def tmrca_consistency(ts1, ts2, id1, id2):
    """
        Given a set of ids calculate the tmrca for these nodes
        NOTE: this requires you have metadata for the individuals  
    """
    # Check that the pairwise tmrcas look very similar
    tmrcas1 = pairwise_tmrcas(ts1, id1, id2)
    tmrcas2 = pairwise_tmrcas(ts2, id1, id2)
    return(np.all(tmrcas1 == tmrcas2))


