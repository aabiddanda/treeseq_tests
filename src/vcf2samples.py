#!/usr/local/bin/python3

import click
import tsinfer
import allel
import numpy as np 
from tqdm import tqdm


@click.command()
@click.option('--vcf', required=True, help='Phased VCF as input')
@click.option('--nalt', default=1, help='Ploidy of the dataset')
@click.option('--shuf', default=False, help='shuffle the order of the individuals')
@click.option('--seed', default=42, help='shuffle the order of the individuals')
@click.option('--threads', default=2, help='num. threads')
@click.option('--ploidy', default=2, help='Ploidy of your samples')
@click.option('--sampleout', default='out.samples', help='sample output')
def main(vcf, nalt, shuf, seed, ploidy, threads, sampleout):
    # 1. Read in the VCF
    vcf_data = allel.read_vcf(vcf, alt_number=nalt)
    gt = vcf_data['calldata/GT'] 
    pos = vcf_data['variants/POS']
    ref = vcf_data['variants/REF']
    alt = vcf_data['variants/ALT']
    ids = vcf_data['samples']
    alleles = np.vstack([ref,alt]).T
    # 2. Converting to a sample file
    np.random.seed(seed)
    nsnps = pos.size
    nids = ids.size
    ind_idx = np.arange(nids)
    if shuf:
        # set the shuffling here
        ind_idx = np.random.permutation(ind_idx)
    progress = tqdm(total=nsnps)
    with tsinfer.SampleData(
        path=sampleout, sequence_length=np.max(pos)+1,
        num_flush_threads=threads) as sample_data:
        for x in ind_idx:
            sample_data.add_individual(ploidy=ploidy, metadata={"name": ids[x]})
        for i in range(nsnps):
            sample_data.add_site(pos[i], gt[i,ind_idx,:].flatten().tolist(), alleles[i,:].tolist())
            progress.update()
    progress.close()
    
if __name__ =='__main__':
    main()
