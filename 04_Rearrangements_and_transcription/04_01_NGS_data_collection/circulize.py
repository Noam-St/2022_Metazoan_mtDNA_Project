#!/gpfs0/biores/apps/Miniconda3/Miniconda_v4.3.21/bin/python3.6

"""
Welcome to circulize, a simple script that creates a new fasta reference from a given fasta reference by translocating the first N (default = 500) bases from the start of the reference to the end, in order to deal with circular genome mapping.
"""
# Change the working directory to the script's location

#%%
import pandas as pd
import os, argparse
PATH = os.getcwd()

def circulize(fasta, output, n = 500):
    """
    Translocate the final 500 characters of a given DNA to the beginning of the file to later combine it (circularity fix)
    """
    with open(fasta, 'r') as fastafile:
        fastalines = fastafile.readlines()
        
        header = fastalines[0].rstrip('\n')
        body = ''.join([i.rstrip('\n') for i in fastalines[1:]])
        body = body[500:] + body[:500]
        output = open(output, 'w')
        output.write(header + '\n')
        output.write(body)

def path_check(path : str) -> str:
  """
  Check whether the path is valid
  """
  if os.path.exists(path):
    return path
  else:
    raise argparse.ArgumentTypeError(f'{path} is not a path to an existing file!')

def abs_path_checl(path : str) -> str:
    """
    Check whether a given path is an absolute path
    """

def main(args):
    circulize(args.fasta, args.output, args.name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = __doc__, prog = 'CirculizeFasta')
    parser.add_argument('--fasta', help = 'full path to input fasta file', metavar = 'fasta_file', type = path_check)
    parser.add_argument('--output', help = 'full path to output fasta file location', metavar = 'output_path', type = str)
    parser.add_argument('-n', '--name', help = 'amount of bases to translocate from the file end to the beginning', metavar = 'N', type = int, default = 500)
    args = parser.parse_args()
    main(args)

    
    
