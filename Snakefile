#!python3

import sys
import tskit, tszip
import msprime as msp
sys.path.append('src/')
from tree_stats import *
from tree_ld_score import gen_ld_scores


# Setting the test VCF file
VCF_DIR = '/home/abiddanda/novembre_lab2/data/external_public/1kg_phase3/haps/'

## ------ 1. Generating Samples ------ ##
rule filter_biallelic:
    """
        Filter to biallelic variants that are at least doubletons 
    """
    input:
        vcf = VCF_DIR + 'ALL.chr{CHROM}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
    output:
        filt_biallelic_vcf =
        'data/real_data/vcf/ALL.chr{CHROM}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.filtmac10.biallelic.vcf.gz'
    shell:
        """
        bcftools norm -m+ --threads 5 {input.vcf} | bcftools view -v snps -m2 -M2 -e \'ALT[*]~\"CN\"\' -c 10:minor  | bgzip -@5 > {output.filt_biallelic_vcf}
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

rule run_tsinfer:
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

TEST_POPS = ['ACB', 'ASW', 'BEB', 'CDX','CEU','CHB', 'CHS', 'CLM', 'ESN', 'FIN', 'GBR', 'GIH', 'GWD', 'IBS', 'ITU', 'JPT','KHV', 'LWK', 'MSL', 'MXL', 'PEL', 'PJL', 'PUR', 'STU', 'TSI', 'YRI']
# print(TEST_POPS)

rule filter_pop:
    input:
        trees = rules.run_tsinfer.output.out_trees
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
        'data/real_data/filt_pops/ld_scores/kg_phase3_chr{CHROM,\d+}.{POP}.ld_score.npy'
    run:
        ts = tszip.decompress(input.treeseq)
        non_inf_sites = [site for site in ts.sites() if len(site.mutations) > 1]
        ts_inf_only = ts.delete_sites([i.id for i in non_inf_sites]).simplify()
        tot_sites = [s for s in ts_inf_only.sites()]
        site_data = np.array([t.position for t in tot_sites], dtype=np.uint32)
        # NOTE : the default is to look at the closest 100 snps in both directions
        ld_scores = gen_ld_scores(ts_inf_only)
        tot_data = np.vstack([site_data, ld_scores]).T
        np.save(output.pop_ld_scores, tot_data)

rule test_ld_scores:
    input:
        expand('data/real_data/filt_pops/ld_scores/kg_phase3_chr22.{POP}.ld_score.npy', POP=TEST_POPS)

# ------ Testing with simulated data ----------- # 
rule sim_constant_data:
  output:
    treeseq = 'data/sim_data/sim_n{n,\d+}_{L,\d+}Mb.simple.rep_{rep,\d+}.trees.tsz'
  run:
    ts = msp.simulate(sample_size=int(wildcards.n), length=int(int(wildcards.L)*1e6), mutation_rate=1e-8, recombination_rate=1e-8, Ne=1e4)
    tszip.compress(ts, output.treeseq)
    
    
rule calc_sim_constant_ld_scores:
    """
        Calculate population-specific LD Scores
    """
    input:
        treeseq = rules.sim_constant_data.output.treeseq
    output:
        pop_ld_scores = 
        'data/sim_data/ld_scores/sim_{n}_{L}.rep_{rep}.ld.npy'
    run:
        ts = tszip.decompress(input.treeseq)
        ld_scores = gen_ld_scores(ts)
        np.save(output.pop_ld_scores, ld_scores)

        
rule sim_data_ld_test:
  input:
    expand('data/sim_data/ld_scores/sim_{n}_{L}.rep_{rep}.ld.npy', n=5000, L=[10],rep=[0])
    
    
