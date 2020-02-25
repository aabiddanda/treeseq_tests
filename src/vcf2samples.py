#!/usr/local/bin/python3

import click
import tsinfer
import allel
from tqdm import tqdm


@click.command()
@click.option('--vcf', required=True, help='Phased VCF as input')
@click.option('--samplesout', default='out.samples', help='')
def main(vcf, samplesout):
    # 1. Read in the VCF
    vcf_data = allel.read_vcf(vcf)

    
if __name__ =='__main__':
	main()
