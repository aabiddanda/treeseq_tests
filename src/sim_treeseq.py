#!/usr/local/bin/python3

'''
    Simulating sampled datasets
'''

import click
import tsinfer 
import msprime as msp
from tqdm import tqdm 

@click.command()
@click.option('--samples', default=1000, help='Number of haplotypes')
@click.option('--length', default=10, help='Length of sequence (Mb)')
@click.option('--Ne', default=10000, help='Effective Population Size')
@click.option('--mu',  default=1e-8, help='Mutation Rate')
@click.option('--rho', default=1.2e-8, help='Recombination Rate')
@click.option('--seed', default=1.2e-8, help='Recombination Rate')
@click.option('--treeout', default='out.trees', help='Output Trees File')
@click.option('--sampleout', default='out.samples', help='Output Samples File')
def main(samples, length, ne, mu, rho, seed, treeout, sampleout):
    # 1. Simulating msprime 
    ts = msp.simulate(sample_size=samples, Ne=ne, recombination_rate=rho, mutation_rate=mu, length=length*(10**6), random_seed=seed)
    ts.dump(treeout)
    # 2. Converting to a sample file
    progress = tqdm(total=ts.num_sites)
    with tsinfer.SampleData(
        path=sampleout, sequence_length=ts.sequence_length,
        num_flush_threads=2) as sample_data:
        for var in ts.variants():
            sample_data.add_site(var.site.position, var.genotypes, var.alleles)
            progress.update()
    progress.close()

if __name__ =='__main__':
    main()
