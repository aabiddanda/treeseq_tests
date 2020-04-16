#!python3

import sys
import tskit
sys.path.append('src/')
from tree_stats import *

# Setting the test VCF file
VCF_DIR = '/home/abiddanda/novembre_lab2/data/external_public/1kg_phase3/haps/'


## ------ 1. Generating Samples ------ ##
rule gen_samples:
  """
    Given a VCF file generate an msprime samples file
  """
  input:
    vcf = VCF_DIR +
    'ALL.chr{CHROM}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
  output:
    'data/real_data/kg_phase3_chr{CHROM,\d+}.samples'
  shell:
    """
      python3 src/vcf2samples.py --vcf {input.vcf}  --threads 5 --sampleout {output}
    """

rule run_tsinfer:
  """
    Run the tsinfer algorithm on our datasets
  """
  input:
    samples = rules.gen_samples.output
  output:
    tmp_trees = temp('data/real_data/kg_phase3_chr{CHROM,\d+}.trees'),
    out_trees = 'data/real_data/kg_phase3_chr{CHROM,\d+}.trees.tsz'
  shell:
    """
      tsinfer infer {input.samples} -t 5 -p  
      tszip -k -f --suffix tsz {output.tmp_trees}
      rm {input.samples}
    """

rule gen_1kg_tree_seq:
  input:
    expand('data/real_data/kg_phase3_chr{CHROM}.trees.tsz', CHROM=22)

           
