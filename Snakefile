#!python3

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
      python3 src/vcf2samples.py --vcf {input.vcf} --seed {wildcards.seed} --threads 2 --sampleout {output}
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

rule gen_trees_all:
    input:
        expand('data/real_data/rep{rep}_noshuf{shuf}_seed{seed}.trees', rep=[1,2], shuf=['False'], seed=42),
        expand('data/real_data/rep{rep}_noshuf{shuf}_seed{seed}.trees', rep=[1,2], shuf=['True'], seed=42),
        expand('data/real_data/rep{rep}_noshuf{shuf}_seed{seed}.trees', rep=[1], shuf=['True'], seed=[100,72]),


## ----- Measuring Similarity ------ ##



