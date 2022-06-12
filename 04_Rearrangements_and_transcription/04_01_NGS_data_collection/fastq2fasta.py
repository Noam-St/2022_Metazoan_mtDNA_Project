#!/gpfs0/biores/apps/Miniconda3/Miniconda_v4.3.21/bin/python3.6

"""
A simple tool to convert FASTQ files to FASTA files using Bio module.
"""

#%%
from Bio import SeqIO
import os, argparse, sys, shutil
PATH = os.getcwd()

#%%
def check_fasta(ref, out, n):
  """
  """
  with open(out, 'r') as outfile:
    with open(ref, 'r') as outref:
      fasta = ''.join([i.rstrip('\n') for i in outfile.readlines()[1:]])
      fasta_ref = ''.join([i.rstrip('\n') for i in outref.readlines()[1:]])
      Ncount = fasta.upper().count('N')
      prop = Ncount/len(fasta)
  if prop >= n:
    print(f'ERROR! The consensus fasta sequence contains too many N characters!\nProportion of N\t{prop}\nOutputting the original sequence instead\n', file = sys.stderr)
    shutil.copy(src = ref, dst = out)
  elif len(fasta_ref) != len(fasta):
    print(f'ERROR! The length of the consensus fasta is shorter than the length of the reference fasta!\nREF LENGTH\t{len(fasta_ref)}\nCON LENGTH\t{len(fasta)}', file = sys.stderr)
  else:
    print('SUCCESS!', file = sys.stderr)

def fastq2fasta(fastq, output):
  """
  Convert fastq file to fasta
  """    
  SeqIO.convert(fastq, 'fastq', output, 'fasta')

def main(args):
  """
  """
  fastq2fasta(args.fastq, args.output)
  if args.fasta_ref != '':
    check_fasta(args.fasta_ref, args.output, args.n)

def path_check(path : str) -> str:
  """
  Check whether the path is valid
  """
  if os.path.exists(path):
    return path
  else:
    raise argparse.ArgumentTypeError(f'{path} is not a valid path!')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description=__doc__, prog = 'fastq2fasta')
  
  parser.add_argument('--fastq', '-f', help = 'fastq input', metavar = 'path_to_fastq', type = path_check)
  parser.add_argument('--output', '-o', help = 'output path', metavar = 'output_path')
  parser.add_argument('--fasta_ref', help = 'original organism fasta file to use incase there are too many n\'s', metavar = 'fasta_path', type = str, default = '')
  parser.add_argument('-n', help = 'Minumum proportion of N in fasta file to cancel the operation and output the normal reference fasta', default = 0.05, type = float)

  args = parser.parse_args()
  main(args)

