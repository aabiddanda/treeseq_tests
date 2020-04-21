"""
    Computing measures of similarity on sets of tree sequences
"""
import tskit
import numpy as np 
from tqdm import tqdm

def gen_ld_scores(ts, nsnps=100):
    """
        Generate LD-Scores for variants in a tree-sequence at a specified width
        NOTE : we do not return position or allelic information which we might
        want to at some point
    """
    ld_calc = tskit.LdCalculator(ts) 
    ld_scores = np.zeros(ts.num_mutations)
    for i in tqdm(range(ts.num_mutations)):
        r2_forward = ld_calc.r2_array(a=i, max_mutations=nsnps,
                direction=tskit.FORWARD)
        r2_reverse = ld_calc.r2_array(a=i, max_mutations=nsnps,
                direction=tskit.REVERSE)
        ld_scores[i] = np.sum(r2_forward) + np.sum(r2_reverse)
    return(ld_scores)

