#!python3

import sys
import tskit, tszip
sys.path.append('src/')
from tree_stats import *
from tree_ld_score import gen_ld_scores


# Setting the test VCF file
VCF_DIR = '/home/abiddanda/novembre_lab2/data/external_public/1kg_phase3/haps/'

## ------ 1. Generating Samples ------ ##
rule filter_biallelic:
    input:
        vcf = VCF_DIR + 'ALL.chr{CHROM}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
    output:
        filt_biallelic_vcf = 'data/real_data/vcf/ALL.chr{CHROM}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.biallelic.vcf.gz'
    shell:
        """
        bcftools norm -m+ --threads 5 {input.vcf} | bcftools view -v snps -m2 -M2 | bgzip -@5 > {output.filt_biallelic_vcf}
        """

rule gen_samples:
  """
    Given a VCF file generate an msprime samples file from only biallelic snps
  """
  input:
    filt_vcf = rules.filter_biallelic.output.filt_biallelic_vcf
  output:
    samples = 'data/real_data/kg_phase3_chr{CHROM,\d+}.samples'
  shell:
    """
      python3 src/vcf2samples.py --vcf {input.filt_vcf} --threads 5  --sampleout {output.samples}
    """

rule run_tsinfer_total:
  """
    Run the tsinfer algorithm on our datasets
  """
  input:
    samples = rules.gen_samples.output.samples
  output:
    tmp_trees = temp('data/real_data/kg_phase3_chr{CHROM,\d+}.trees'),
    out_trees = 'data/real_data/kg_phase3_chr{CHROM,\d+}.trees.tsz'
  shell:
    """
      tsinfer infer {input.samples} -t 5 -p 
      tszip -k -f --suffix .tsz {output.tmp_trees}
    """

rule gen_1kg_tree_seq:
  input:
    'data/real_data/kg_phase3_chr22.trees.tsz'


## ----- 2. Setting up filtering procedures ------ ##
pop_panel_file = 'params/integrated_call_samples_v3.20130502.ALL.panel' 
pop_panel = np.loadtxt(pop_panel_file, dtype=str)
indiv_ids = pop_panel[1:,0]
pop_ids = pop_panel[1:,1]
popdict = {}
for i in range(indiv_ids.size):
    popdict[indiv_ids[i]] = pop_ids[i]


rule filter_pop:
    input:
        trees = rules.run_tsinfer_total.output.out_trees
    output:
        trees_filt =
        'data/real_data/filt_pops/kg_phase3_chr{CHROM,\d+}.{POP}.trees.tsz'
    run:
        ts = tszip.decompress(input.trees)
        filt_ts = filter_pop(ts, popdict=popdict, pop=wildcards.POP)
        tszip.compress(filt_ts, output.trees_filt)


rule calc_pop_ld_scores:
    """
        Calculate population-specific LD Scores
    """
    input:
        treeseq = rules.filter_pop.output.trees_filt
    output:
        pop_ld_scores = 
        'data/real_data/filt_pops/ld_scores/kg_phase3_chr{CHROM,\d+}.{POP}.ld.npy'
    run:
        ts = tszip.decompress(input.treeseq)
        ld_scores = gen_ld_scores(ts)
        np.save(output.pop_ld_scores, ld_scores)

rule test:
    input:
        'data/real_data/filt_pops/ld_scores/kg_phase3_chr22.CEU.ld.npy'

        
