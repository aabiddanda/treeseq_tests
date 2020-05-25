# Performance Testing for LD generation

## Installation / packages

In order to run the python and bash scripts you will need to install the following python packages via pip:

```
pip install tskit tszip click --user --upgrade
```

And load the `plink` package using `module load plink`

## Running LD inference via plink

In order to obtain times you can run the following:

```
time bash local_ld_mat_plink.sh ../../data/real_data/vcf/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.filtmac10.biallelic.vcf.gz <test.pos> 200000 > /dev/null
```

Note that `<test.pos>` is a test_position from the `test.pos` file (which contains 50 randomly chosen snp positions to center a window on)

## Running LD inference via tskit

If you are running from this directory you can run the following:

```
time python local_ld_mat_tskit.py --tsz ../../data/real_data/kg_phase3_chr22.trees.tsz --position <test.pos> --width 200000 
```

Where again `<test.pos>` is a position from the randomly `test.pos` file.

## Simulation experiment

For each position in the `test.pos`, we should generate two additional columns: (1) for the mean time it takes to run the LD estimation in plink and (2) for the mean time it takes to run the LD estimation using tskit















