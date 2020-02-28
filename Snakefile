#!python3

import sys
import tskit
sys.path.append('src/')
from tree_stats import *

# Setting the test VCF file
vcf_file = 'data/real_data/ceu_panel_chr22.nonmono.vcf.gz'

## ------ 1. Generating Samples ------ ##
rule gen_samples:
  """
    Given a VCF file generate an msprime samples file
  """
  input:
    vcf = vcf_file
  wildcard_constraints:
    shuf = '(True|False)'
  output:
    'data/real_data/rep{rep,\d+}_noshuf{shuf}_seed{seed,\d+}.samples'
  shell:
    """
      python3 src/vcf2samples.py --vcf {input.vcf} --shuf {wildcards.shuf} --seed {wildcards.seed} --threads 2 --sampleout {output}
    """

rule run_tsinfer:
  """
    Run the tsinfer algorithm on our datasets
  """
  input:
    samples = rules.gen_samples.output
  output:
    'data/real_data/rep{rep,\d+}_noshuf{shuf}_seed{seed,\d+}.trees'
  shell:
    """
      tsinfer infer {input.samples} -t 2 -p 
    """

## ----- 2. Measuring similarity between trees  ------ ##

rule sim_report:
    """
        Generate a report on tree-similarity 
    """
    input:
        tree1 = 'data/real_data/rep{rep1}_noshuf{shuf1}_seed{seed}.trees',
        tree2 = 'data/real_data/rep{rep2}_noshuf{shuf2}_seed{seed}.trees'
    output:
       sim_report = 'data/sim_report/rep_{rep1}_{rep2}_shuf_{shuf1}_{shuf2}_seed{seed}.txt'
    run:
      ts1 = tskit.load(input.tree1)
      ts2 = tskit.load(input.tree2)
      with open(output.sim_report, 'w+') as f:
        # Checking for similarity in breakpoints
        num_trees = ntrees_consist(ts1,ts2)
        print('Trees:', num_trees, file=f)
        bp_sim = breakpoint_consist(ts1,ts2)
        print('Breakpoints:',bp_sim, file=f)
            
rule final_test:
    """
        Rule to perform checks on tree-sequence structure
    """
    input:
        expand('data/sim_report/rep_{rep1}_{rep2}_shuf_{shuf1}_{shuf2}_seed{seed}.txt', rep1=1, rep2=[1,2], shuf1='False', shuf2=['True','False'], seed=[42])
 


