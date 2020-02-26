#!/usr/local/bin/python3

import click
import tsinfer
import allel
import numpy as np 
from tqdm import tqdm


@click.command()
@click.option('--vcf', required=True, help='Phased VCF as input')
@click.option('--nalt', default=1, help='Ploidy of the dataset')
@click.option('--threads', default=2, help='Ploidy of the dataset')
@click.option('--sampleout', default='out.samples', help='sample output')
def main(vcf, nalt, threads, sampleout):
    # 1. Read in the VCF
    vcf_data = allel.read_vcf(vcf, alt_number=nalt)
    gt = vcf_data['calldata/GT'] 
    pos = vcf_data['variants/POS']
    ref = vcf_data['variants/REF']
    alt = vcf_data['variants/ALT']
    alleles = np.vstack([ref,alt]).T
    # 2. Converting to a sample file
    nsnps = pos.size
    progress = tqdm(total=nsnps)
    with tsinfer.SampleData(
        path=sampleout, sequence_length=np.max(pos)+1,
        num_flush_threads=threads) as sample_data:
        for i in range(nsnps):
            sample_data.add_site(pos[i], gt[i,:,:].flatten().tolist(), alleles[i,:].tolist())
            progress.update()
    progress.close()
    
if __name__ =='__main__':
    main()
