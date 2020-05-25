#!/usr/local/bin/python3

'''
	Template python script with arguments
'''

import tszip, tskit
import click

@click.command()
@click.option('--tsz', required=True, help='Tskit tree sequence file')
@click.option('--position', type=float, required=True, help='position for middle of window')
@click.option('--width', type=float, required=True, default=2e5, help='Window position')
def main(tsz, position, width):
  # Read in the compressed tree-sequence
  ts = tszip.decompress(tsz)
  # Filter to a particular region to compute LD
  region = [[(position - width/2), (position + width/2)]]
  ts_region = ts.keep_intervals(region)
  # Compute LD (r2) matrix between all variants in the region
  Ldcalc = tskit.LdCalculator(ts_region)
  r2_mat = Ldcalc.r2_matrix()
  
  
if __name__ =='__main__':
	main()
