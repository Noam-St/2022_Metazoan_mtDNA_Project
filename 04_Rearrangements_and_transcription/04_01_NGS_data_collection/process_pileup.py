"""
Process a pileup format file into an annotated .csv file where each row is a single position, the columns are: position, coverage, z_normalized coverage, read_ends, feature, gene, length, strand
"""

import pandas as pd
import numpy as np
import os, sys, argparse
from Bio import SeqIO
from Bio import SeqFeature
from normalize_from_htseq import import_htseq

def dir_checker(dirname : str) -> None:
  """
  check if a directory exists, if it does not, create it.
  Parameters
  ----------
  dirname : str
    Name of the directory to be checked
  """
  if not os.path.isdir(os.path.join(PATH, dirname)):
    os.mkdir(os.path.join(PATH, dirname))

def dup_assigner(gene : str, genelist : list) -> str:
  """
  DEPRECATED
  """
  return gene + '*' * sum([gene in i for i in genelist])    


def annotate(org_name : str, record) -> pd.DataFrame:
  """
  Create a df of position and current feature and product annotation for each position
  
  Parameters
  ----------
  org_name : str
    Name of organism to work on
  
  Returns
  -------
  return_df : pd.DataFrame
    A df with the following columns: Position, Feature and Gene.
  
  """
  rrna = {'mgr01' : 'RNR1', 'mgr02' : 'RNR2'}
  if os.path.exists(os.path.join(PATH, 'annotated_df',f'{org_name}.csv')):
    return pd.read_csv(os.path.join(PATH, 'annotated_df',f'{org_name}.csv'), index_col = 0)
  features = record.features
  seq = record.seq
  return_df = pd.DataFrame(columns = ['Position','Feature','Gene','Length']) #The organism specific df
  return_df['Position'] = list(range(1,len(seq)+1))
  for feature in features: #Iterate over organism features
    if feature.type == 'gene': 
      gene = feature.qualifiers.get('gene')
      if gene:
        gene = gene[0]
      else:
        gene = feature.qualifiers.get('locus_tag')[0]
        try: gene = rrna[gene.split('_')[1]]
        except (IndexError, ValueError, KeyError) : continue
      for location in feature.location.parts:
        start = location.start
        end = location.end
        strand = location.strand
        #gene = dup_assigner(gene, genelist)
        #genelist.append(gene)
        return_df.loc[start:end, 'Gene'] = gene
        return_df.loc[start:end, 'Feature'] = feature.type
        return_df.loc[start:end, 'Strand'] = True if strand == 1 else False
        return_df.loc[start:end, 'Length'] = end - start
  dir_checker('annotated_df')
  return_df.to_csv(os.path.join(PATH, 'annotated_df', f'{org_name}.csv'))
  return return_df

def combine_with_sample(sample_df : pd.DataFrame, annotated_df : pd.DataFrame) -> pd.DataFrame:
  """
  Combine the annotation_df with sample_df, include all columns from both dfs and fill missing values from sample_df with 0
  Parameters
  ----------
  sample_df : pd.DataFrame
    The specific sample df with read coverage data.
  annotated_df : pd.DataFrame
    The df with all the gene locations for the specific organism.
  """
  merged = pd.merge(left = sample_df, right = annotated_df, on = 'Position', how = 'outer', suffixes = (False, False)).sort_values(by='Position', ascending = True)
  merged.coverage.fillna(0, inplace = True)
  return merged

def z_score(cov :int, mean : int, std : int) -> int:
  """
  Calculate Z score for given coverage
  """
  return (cov - mean)/std

def coverage(sample : pd.DataFrame, total_reads : int, flip_strands : bool) -> pd.DataFrame:
  """
  Measure normalized versions of the coverage (available from pileup file) and the percentage of reads that end in position

  Parameters
  ----------
  sample : sample DF.

  Returns
  -------
  sample : sample DF modified with the new cols.
  """
  #Count exact amount of exact ref matches, indels and snps as coverage.
  #IMPORTANT: CURRENTLY SET AS: FORWARD = Negative, REVERSE = Positive
  sample['ends_ratio'] = 0
  sample['neg_ends_ratio'] = 0
  sample['pos_ends_ratio'] = 0

  #Calculates the amount of reads that end in location in each strand
  if flip_strands: #If the strand is flipped, the reads are counted in the opposite direction
    pos_count_sign = ','
    pos_snp_sign = 'actgn'
    neg_count_sign = '.'
    neg_snp_sign = 'ACTGN'
  else: # This is in normal mode
    pos_count_sign = '.'
    pos_snp_sign = 'ACTGN'
    neg_count_sign = ','
    neg_snp_sign = 'actgn'
  
  sample['neg_read_end'] = sample.Location.str.count(rf'\$\{neg_count_sign}|\$[{neg_snp_sign}]') # Counts the number of reads that end in the negative strand for each position
  sample['pos_read_end'] = sample.Location.str.count(rf'\$\{pos_count_sign}|\$[{pos_snp_sign}]') # Counts the number of reads that end in the positive strand for each position
  sample['read_end'] = sample['pos_read_end'] + sample['neg_read_end'] # Sums the number of reads that end in each strand for the total number of reads that end in the location
  sample['read_end'].fillna(0) # Fills NaN with 0
  
  #Calculate the number of reads that either end or start in location in both strands
  sample['pos_end_start_counts'] = sample.Location.str.count(rf'(\$\{pos_count_sign}|\$[{pos_snp_sign}])|(\^\W\{pos_count_sign}|\$\W[{pos_snp_sign}])')
  sample['neg_end_start_counts'] = sample.Location.str.count(rf'(\$\{neg_count_sign}|\$[{neg_snp_sign}])|(\^\W\{neg_count_sign}|\$\W[{neg_snp_sign}])')
  sample['end_start_counts'] = sample['pos_end_start_counts'] + sample['neg_end_start_counts']
 
  sample.coverage = sample.coverage.fillna(0) # Fills NaN with 0
  sample['neg_coverage'] = sample.Location.str.count(rf'\{neg_count_sign}|[{neg_snp_sign}]|[\-\+][0-9]+[{neg_snp_sign}]+') # Counts the number of reads assigned to the negative strand
  sample['pos_coverage'] = sample.Location.str.count(rf'\{pos_count_sign}|[{pos_snp_sign}]|[\-\+][0-9]+[{pos_snp_sign}]+') # Counts the number of reads assigned to the positive strand
 
  #Count the number of reads that end in position and normalized to total read_ends
  pos_total = np.sum(sample['pos_read_end']) # Total number of reads that end in the positive strand
  neg_total = np.sum(sample['neg_read_end']) # Total number of reads that end in the negative strand
  sample['read_end'] = sample['read_end']/(pos_total + neg_total) # Normalize the read_end to the total number of reads that end in the sample

  #Calculate RPM based on the total_reads attained by htseq-count
  sample['RPM'] = sample['coverage']*1000000/total_reads
  sample['pos_RPM'] = sample['pos_coverage']*1000000/(pos_total)
  sample['neg_RPM'] = sample['neg_coverage']*1000000/(neg_total)
  
  #Calculate Z-scores
  tot_mean = np.mean(sample['coverage'])
  tot_std = np.std(sample['coverage'])
  pos_mean = np.mean(sample['pos_coverage'])
  neg_mean = np.mean(sample['neg_coverage'])
  pos_std = np.std(sample['pos_coverage'])
  neg_std = np.std(sample['neg_coverage'])
  sample['pos_z'] = sample.pos_coverage.apply(z_score, args = [pos_mean, pos_std])
  sample['neg_z'] = sample.neg_coverage.apply(z_score, args = [neg_mean, neg_std])
  sample['z'] = sample.coverage.apply(z_score, args = [tot_mean, tot_std])

  #Calculate ends and starts ratio only for reads that have more than mean - 1.5*std coverage
  sample.loc[sample['coverage'] > (tot_mean - tot_std*1), 'end_start_ratio'] = (sample['end_start_counts'] / sample['coverage'])
  sample.loc[sample['pos_coverage'] > (pos_mean - pos_std*1), 'pos_end_start_ratio'] = (sample['pos_end_start_counts'] / sample['pos_coverage'])
  sample.loc[sample['neg_coverage'] > (neg_mean - neg_std*1), 'neg_end_start_ratio'] = (sample['neg_end_start_counts'] / sample['neg_coverage'])

  return sample

def record_check(org : str, org_df : pd.DataFrame) -> None:
  """
  Check whether org genbank exists, if not download it. using a wget command (WORKS ONLY ON LINUX)

  Parameters
  ----------
  org : str
    The organism to check record for
  org_df : pd.DataFrame
    A database of all organisms used to retrieve RefSeq ID
  """

  ID = org_df.loc[org, 'RefSeq']
  if type(ID) != str: ID = ID.iloc[0]
  filename = os.path.join(PATH,'genbank_DB', f'{org}.gbk')
  if not os.path.isfile(filename):
    os.system(f'wget -O {filename} \"https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=genbank&retmode=txt&id={ID}\"\n')

def n_score(row, sample, window):
  """
  FOS = (C+1)/L + (C+1)/R
  C = average number of reads over current window.
  L = average number of reads over left-side window
  R = average number of reads over right-side window 
  """
  cur_pos = row.Position

def main(args : argparse.ArgumentParser) -> None:
  """
  Combine all above functions to create an annotated .csv based on .pileup file path and org_name
  """
  global PATH
  PATH = args.wd
  try:
    org_df = pd.read_csv(os.path.join(PATH, 'final.csv'))
    org_df.organism = org_df.organism.str.replace(' ','_')
    org_df = org_df.set_index('organism')
  except FileNotFoundError:
    raise FileNotFoundError('final.csv is required in the src folder for the program to operate.')

  org_name = args.org
  dir_checker('genbank_DB')
  dir_checker('annotated_pileup')
  dir_checker(os.path.join('annotated_pileup', org_name))
  filename = os.path.join(PATH, 'genbank_DB', f'{org_name}.gbk')
  record_check(org_name, org_df)

  record = SeqIO.read(filename, 'genbank')
  sample = pd.read_csv(args.pileup, delimiter = '\t',names = ['Chromosome','Position','Base','coverage','Location','Quallity'], index_col = None, usecols = list(range(0,6)))
  counts, no_feature = import_htseq(args.counts)
  total_reads = counts.counts.sum() + no_feature
  sample = coverage(sample, total_reads, args.flip_strands)
  sample = sample.drop(columns = ['Location','Quallity'])
  annotated = annotate(org_name, record)
  final = combine_with_sample(sample, annotated)
  final['dataset'] = os.path.basename(os.path.dirname(os.path.dirname(args.pileup)))
  accs = os.path.basename(args.pileup).replace(".pileup","")
  if args.output == '':
    final.to_csv(os.path.join(PATH, 'annotated_pileup', org_name, f'{accs}.csv'))
  else:
    final.to_csv(args.output)

def path_check(path : str) -> str:
  """
  Check whether the path is valid
  """
  if os.path.exists(path):
    return path
  else:
    raise argparse.ArgumentTypeError(f'{path} is not a valid path!')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description = __doc__, prog = 'PileupReader', epilog = 'This program automatically outputs the files into annotated_pileup folder in cwd')
  parser.add_argument('--pileup', help='Path to the pileup file', metavar='pileup_path', type = path_check)
  parser.add_argument('--org', help = 'Organism name', type = str, metavar='organism')
  parser.add_argument('--version', help = 'Version information', action = 'version', version = '%(prog)s version 1.0.0')
  parser.add_argument('--wd', help = 'Set the working dir to create directories in', default = os.getcwd(), type = path_check, metavar = 'WorkingDirectory')
  parser.add_argument('--counts', help = 'Path to htseq output', type = path_check, metavar = 'Counts')
  parser.add_argument('--output', help = 'output csv location', metavar = 'output', default = '')
  parser.add_argument('--flip-strands', help = 'Flip the heavy and light strands (mostly for PRO-seq data)', default = False, action = 'store_true')
  args = parser.parse_args()
  main(args)