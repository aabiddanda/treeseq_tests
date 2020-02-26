# tree_seq_testing
----------------------------

Testing of tree sequence operations on both simulated and real datasets

### Simulated Data

To simulate both `.trees` and `.samples` files that are necessary for establishing ground truth as well inferring tree topology and tree length. 

The command to simulate trees using the python script is:

```
python3 src/sim_treeseq.py --samples <num samples> --length <length in Mb> --mu <per bp mutation rate> --rho <per bp recombination rate> --treeout <output tree file> --sampleout <sample output file>
```

### Real Data

To convert a VCF to a samples file we have:

```
python3 src/vcf2samples.py --vcf <input VCF> --sampleout <output samples>
```
