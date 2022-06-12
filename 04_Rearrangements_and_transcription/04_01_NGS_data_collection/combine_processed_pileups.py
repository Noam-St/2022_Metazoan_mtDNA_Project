#!/gpfs0/biores/apps/Miniconda3/Miniconda_v4.3.21/bin/python3.6

"""
Welcome to combine_processed_pileups.py script!
This script combines two processed_pileup.py output csvs together by moving to_return amount of bases from the fixed_path csv to the beginning of the normal_path csv file.
The main idea is to create a single csv file in output that is fixed for the circularity of the mtDNA. 
"""
#%%
import pandas as pd
import numpy as np
import os, argparse
import matplotlib.pyplot as plt
PATH = os.getcwd()
#%%

def combine_pileups(n_path, c_path, output, n_return, n = 500):
    """
    Combine two different processed_pileup csv files and output it to a given path

    Parameters
    ----------
    n_path : str
      Path to the normal processed_pileup csv
    c_path : str
      Path to the fixed (for circularity) processed_pileup csv
    output : str
      Path to output the new combined csv
    n_return : int
      Amount of bases to move from one csv to another (fixed to normal)
    n : int
      Amount of bases moved from the beginning of the initial fasta ref to its end
    
    Raises
    ------
    ValueError :
      For the rare case when the initial amount (n) happens to be bigger than n_return, the amount of bases to move from fixed processed pileup csv to the normal processed pileup csv
    """
    if n < n_return:
        raise ValueError(f'The amount of bases returned must be smaller or equal to the amount of bases initially translocated (default = 500)\nValue provided\t{n_return}')
    n_diff = n - n_return
    normal = pd.read_csv(n_path, index_col=0)
    circ = pd.read_csv(c_path, index_col=0)

    combined = circ.loc[(len(circ) - n < circ.Position) & (circ.Position <= len(circ) - n_diff)]
    combined.Position = [i for i in range(1, n_return + 1)]
    combined = combined.append(normal.loc[normal.Position > n_return]).reset_index().drop(columns = 'index')
    to_fix = ['Feature', 'Gene', 'Length', 'Strand']
    normal = normal.set_index('Position')
    combined = combined.set_index('Position')
    for feature in to_fix:
      combined.loc[:, feature] = normal[feature]
    combined.reset_index().to_csv(output)

def main(args):
    """
    """
    combine_pileups(args.normal_path, args.fixed_path, args.output, args.to_return, args.n)

def path_check(path : str) -> str:
  """
  Check whether the path is valid
  """
  if os.path.exists(path):
    return path
  else:
    raise argparse.ArgumentTypeError(f'{path} is not a valid path!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, prog = 'combine_processed_pileups')

    parser.add_argument('--normal_path', help = 'Path to the location of the normal processed pileup csv', type = path_check, metavar = 'normal_csv_path')
    parser.add_argument('--fixed_path', help = 'Path to the location of the fixed processed pileup csv', type = path_check, metavar = 'fixed_csv_path')
    parser.add_argument('--to_return', help = 'Amount of bases to move from the fixed csv to the normal csv beginning', type = int, metavar = 'n_bases', default = 400)
    parser.add_argument('-n', help = 'Amount of bases initially translocated from the end of the fasta file to the beginning', type = int, metavar = 'n_bases', default = 500)
    parser.add_argument('--output', help = 'Path to output', type = str, metavar = 'output')
    args = parser.parse_args()
    main(args)



