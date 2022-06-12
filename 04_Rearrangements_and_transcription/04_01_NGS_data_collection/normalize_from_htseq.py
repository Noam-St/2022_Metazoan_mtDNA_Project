"""
Calculate TPM based on a htseq-counts output and a GFF3 file
THIS CODE CURRENTLY ONLY WORKS WITH GFF3 IN NCBI FORMAT
usage: python tpm_from_htseq.py --counts [path to htseq output] --gff [path to gff file] output [path_to_output]
"""

import argparse, os
import pandas as pd
PATH = os.getcwd()

def import_htseq(counts_path : str) -> pd.DataFrame:
  """
  Parse htseq-counts output into a df

  Parameters
  ----------
  counts_path : str
    Absolute path to the htseq-counts output file (in any format)
  
  Returns
  -------
  counts : pd.DataFrame
    A df where each row is a single gene and the columns contain counts, gene name (if available) and gene ID
  """
  rrna_dict = {'mgr01': 'RNR1', 'mgr02': 'RNR2'}
  try:
    counts = pd.read_csv(counts_path, names = ['ID','gene', 'counts'], delimiter = '\t')
    if any(counts['counts'].isnull()):
      counts.columns = ['ID','counts']
  except ValueError:
    counts = pd.read_csv(counts_path, names = ['ID','counts'], delimiter = '\t')
  no_feature = counts.loc[(counts.ID.str.contains('__')) & (counts.ID != '__not_aligned'), 'counts'].sum()
  counts = counts.loc[~counts.ID.str.contains('__'), :]
  for i, row in counts.iterrows():
    if 'mgr' in row.ID:
      try:
        rrna = row.ID.split('_')[1]
        counts.at[i, 'gene'] = rrna_dict[rrna]
      except (IndexError, KeyError) : continue

  return counts, no_feature

def extract_ID(att : str, format : str) -> str:
    if format == 'gff3': att = att.split(';')[0].replace('ID=','')
    else: att = att.split(' ')[1].replace('\"','').replace(';','')
    return att

def import_gff(gff_path : str, format : str) -> pd.DataFrame:
  """
  Parse a gff into a df, where each row is a gene and the cols are the sample name, ID, type of gene, start location and end location.

  Parameters
  ----------
  gff_path : str
    The absolute path to the gff3 file.
  
  Returns
  -------
  gff : pd.DataFrame
    a dataframe of the parsed gff3 file, columns: sample, type, start, end, ID. Only gene types are kept.
  """
  if format == 'gff3':
    with open(gff_path, 'r') as gff_file:
      lines = gff_file.readlines()
      skip = 3 if 'region' in lines[2] else 2
    gff = pd.read_csv(gff_path, delimiter = '\t', usecols = [0,2,3,4,8], skiprows = skip, names = ['sample', 'type', 'start', 'end', 'ID'])
    gff = gff.loc[gff.type == 'gene', :]
    gff.ID = gff.ID.apply(extract_ID, format = format)
    gff['glengths'] = (gff.end - gff.start) + 1

    try: mtlength = int(lines[0].split(' ')[-1].rstrip('\n'))
    except ValueError:
      print('The first line of the gff3 file must end with the gene length!\n')
      mtlength = 0
  else:
      with open(gff_path, 'r') as gff_file:
        skip = 0
        for line in gff_file:
          if '\t' not in line:
            skip +=1
          else: break
      gff = pd.read_csv(gff_path, delimiter = '\t', usecols = [0,2,3,4,8], skiprows = skip, names = ['sample', 'type', 'start', 'end', 'ID'])
      gff = gff.loc[gff.type == 'gene', :]
      gff['ID'] = gff.ID.apply(extract_ID, format = format)
      gff['glengths'] = (gff.end - gff.start) + 1

      mtlength = 0
  return gff, mtlength
  
def combine_lengths_counts(counts : pd.DataFrame, gff : pd.DataFrame) -> pd.DataFrame:
  """
  Merge the two dataframes (gff and counts)
  """
  return pd.merge(counts, gff, on = 'ID').drop(columns = ['type', 'start', 'end'])

def TPM(combined_df : pd.DataFrame, no_feature : int, mtlength : int) -> pd.DataFrame:
  """
  Calculate TPM for each gene on the merged dataframe

  Parameters
  ----------
  combined_df : pd.DataFrame
    The combined dataframe which contains data from both counts and gff files.
  
  Returns
  -------
  combined_df : pd.DataFrame
    The same dataframe with new columns for kbp, rpk and tpm
  """
  combined_df['kbp'] = combined_df['glengths'] / 1000
  combined_df['rpk'] = combined_df['counts'] / combined_df['kbp']
  no_feature_kbp = (mtlength - combined_df.glengths.sum())/1000
  total_rpk = combined_df.rpk.sum() + no_feature / no_feature_kbp
  if total_rpk == 0:
    print(f'There is an issue with this file! {combined_df.head()}')
    total_rpk = 1
  
  combined_df['tpm'] = combined_df['rpk'] / (total_rpk / 1000000)
  return combined_df

def RPKM(combined_df : pd.DataFrame, no_feature : int) -> pd.DataFrame:
  """
  Calculate RPKM for each gene on the merged dataframe

  Parameters
  ----------
  combined_df : pd.DataFrame
    The combined dataframe which contains data from both counts and gff files.
  
  Returns
  -------
  combined_df : pd.DataFrame
    The same dataframe with new columns for rpm, kbp and rpkm
  """
  total_reads = combined_df['counts'].sum() + no_feature
  if total_reads == 0:
    print(f'There is an issue with this file! {combined_df.loc[0, "ID"]}')
    total_reads = 1
  combined_df['rpm'] = combined_df['counts'] / (total_reads / 1000000)
  combined_df['kbp'] = combined_df['glengths'] / 1000
  combined_df['rpkm'] = combined_df['rpm'] / combined_df['kbp']

  return combined_df

def main(args : argparse.ArgumentParser):
  """
  Parse htseq-counts output, gff3 file and combine the, into a unified dataframe, calculate either tpm or rpkm and create a csv file in the chosen location.

  Parameters
  ----------
  args : argparse.ArguementParser
    An arguementparser object with all the arguements the user inputs in the command
  """
  format = args.format
  counts, no_feature = import_htseq(args.counts)
  gff, mtlengths = import_gff(args.gtf, format)
  if args.mtlength != 0: # If mtlength is defined, use that instead of the gff file
    mtlengths = args.mtlength
  elif mtlengths == 0: # If mtlength is not defined and the gff file is empty, print error and default to human mtDNA length
    print('Unable to parse the ref size from gtf and mtlength is 0. Please provide a ref size (defualting to human mtDNA).')
    mtlengths = 16569
  combined = combine_lengths_counts(counts, gff)
  if args.mode == 'tpm':
    final = TPM(combined, no_feature, mtlengths)
  elif args.mode == 'rpkm':
    final = RPKM(combined, no_feature)
  else:
    final = TPM(combined, no_feature, mtlengths)
    final = RPKM(combined, no_feature)

  final.to_csv(args.output)

def main_for_import(counts_path, gtf_path):
  """
  """
  counts, no_feature = import_htseq(counts_path)
  gff, mtlength = import_gff(gtf_path)
  combined = combine_lengths_counts(counts, gff)
  final = TPM(combined, no_feature, mtlength)
  return final

def path_check(path : str) -> str:
  """
  Check whether the path is valid
  """
  if os.path.isfile(path):
    return path
  else:
    raise argparse.ArgumentTypeError(f'{path} is not a valid path!')

def tpm_rpkm_check(string : str) -> str:
  """
  Check whether the str input for --mode is rpkm or tpm
  """
  if string not in ['tpm', 'rpkm']:
    raise argparse.ArgumentTypeError(f'{string} is not a valid option!')
  else:
    return string

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description = __doc__, prog = 'HTSeq-countsNormalizer')
  parser.add_argument('--counts','-c', help = 'Path to htseq-counts output file', type = path_check, metavar='counts_path')
  parser.add_argument('--gtf', '--gff', '-g', help = 'Path to gff3 output file', type = path_check,metavar='gff_path')
  parser.add_argument('--format', '-f', help = 'Format of the gff3 file. Default is gff3', default = 'gff3', choices = ['gff3', 'gtf'])
  parser.add_argument('--output', '-o',  help = 'Output path', metavar = 'output', default = './output.csv')
  parser.add_argument('--mode', help = 'Input mode (rpkm, tpm or both)', default = 'both', choices = ['tpm','rpkm', 'both'], metavar=('calc_type'))
  parser.add_argument('--mtlength', help = 'The length of the mitochondrial genome (default = try to parse from gff)', metavar = 'mtlength', default = 0, type = int)
  parser.add_argument('--version', action='version', version = '%(prog)s 1.0.0')

  args = parser.parse_args()
  main(args)