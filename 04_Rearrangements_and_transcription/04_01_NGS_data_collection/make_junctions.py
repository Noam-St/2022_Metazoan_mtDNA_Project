
"""
FLANKING REGIONS CODE
This code filters the output of annotation.py with the desired flanking regions
TODO: Make this code a stand-alone unix pipeline step.
"""

import numpy as np
import pandas as pd
import os, re, argparse, utils
from importlib import reload
reload(utils)
def num_sim(n1 : int, n2 : int) -> int:
  """ calculates a similarity score between 2 numbers """
  return 1 - abs(n1 - n2) / (n1 + n2)

def genomic_ranges(sample : pd.DataFrame, no_trna : bool) -> list:
    """
    Based on a sample df made by the annotation.py script, return a list in the following fromat:
    [start, end, gene_name]
    Parameters
    ----------
    sample : pd.DataFrame
      A dataframe made using the annotation.py code. MUST HAVE COLUMNS:
      Gene, Position
    not_trna : bool
      True if junctions ignore tRNAs 
    Returns
    -------
    granges : list
      List of all the granges in the following format:
      [start, end, gene_name]
    """
    no_na = sample.Gene.unique().tolist()
    try: no_na.remove(np.nan)
    except ValueError: pass
    mtdna_len = len(sample)
    if no_trna:
      genelist = [gene for gene in no_na if 'trn' not in gene.lower()]     
    else:
      genelist = no_na
    granges = []
    for gene in genelist:
        cur_gene = sample.loc[sample.Gene == gene, :]
        if len(cur_gene) < 50: continue
        cur_pos = cur_gene.Position.to_list()
        zero_ind = cur_pos[0]
        last_ind = cur_pos[-1]
        size = abs(last_ind - zero_ind)
        if  size > (mtdna_len* 0.9): #This indicates that the gene happens to be 
          print(size)
          last_ind =  max([i for i in cur_pos if i < (mtdna_len* 0.9)])
          print(last_ind)
        granges.append([zero_ind, last_ind, gene])
    return granges

def wrap(start, end, pos, flank, size):
    """
    If the current position + the flanking region or the current position - the flanking region exceeds the 
    boundaries of the mtDNA, wrap around (because the mtDNA is circular)

    Parameters
    ----------
    start : int
      Gene start position
    end : int
      Gene ending position
    pos : int
      Current position (in the loop)
    flank : int
      The size of the flanking region we want
    size : int
      The mtDNA genome size
    
    Returns 
    -------
    bool
      A boolean that is True if the current position wraps around the mtDNA and is flanking.
    """
    if start - flank < 0:
        return pos >= start - flank + size
    elif start + flank > size:
        return pos < start + flank - size
    elif end - flank < 0:
        return pos >= end - flank + size
    elif end + flank > size:
        return pos < end + flank - size
    else: return False
  
def downstream_wrap(end : int, pos : int, flank : int, size : int) -> bool:
  """
  Wrap to the start of the genome incase end + flank is bigger than genome size (mtDNA is circular)
  return True if pos within range and False if not
  """
  if end + flank > size:
    return pos < end + flank - size   
  else: return False

def flanking(pos : int, genelist : list, flank : int, gsize : int, direction = 'downstream') -> bool:
    """
    Receive a genelist, return True if the current position is within flanking range (including wraps)
    
    Parameters
    ----------
    pos : int
      The current position to check
    genelist : list
      The list to compare the position with, MUST be in the following format: [start, end, gene_name]
    flank : int
      The flanking region size we are interested in.
    gsize : int
      The size of the mtDNA genome.
    
    Returns
    -------
    bool
      True if the pos is within flanking range.
    """
    if direction == 'both':
      return any([(grange[0] - flank <= pos < grange[0] + flank) or (grange[1] - flank <= pos < grange[1] + flank) or (wrap(grange[0], grange[1], pos, flank, gsize)) for grange in genelist])
    else:
      return any([(grange[1] < pos < grange[1] + flank) or (downstream_wrap(grange[1], pos, flank, gsize)) for grange in genelist])


def neigh_freq_and_expression(sample : pd.DataFrame, tpm : pd.DataFrame, prob_matrix : pd.DataFrame, granges : list, flank : int, org : str, read_size : int) -> pd.DataFrame:
  """
  Calculate junction-specific parameters such as junction size, gene expression on both sides, coexpression and the frequency of the neighbouring genes

  Parameters
  ----------
  sample : pd.DataFrame
    The input sample df
  tpm : pd.DataFrame
    The expression measurements df for the current sample
  prob_matrix : pd.DataFrame
    A dataframe of the contact frequencies for current phylum
  granges : list
    A list of granges created by created by genomic_ranges function
  flank : int 
    Flank size
  org : str 
    Current org name
  
  Returns
  -------
  sample : pd.DataFrame
    The current sample dataframe modified with the extra parameters
  """
  sample['neigh_freq'] = 0
  sample['strand_switch'] = np.nan
  sample['gpair'] = np.nan
  sample['coexpression'] = np.nan
  sample['left_gene'] = np.nan
  sample['right_gene'] = np.nan
  sample['left_tpm'] = np.nan
  sample['right_tpm'] = np.nan
  sample['org'] = org
  sample['junction_size'] = np.nan
  sample['sum_around_junction'] = np.nan

  for i, grange in enumerate(granges):
    cur_granges = []
    if grange[0] == granges[-1][0]:
      print(granges[-1][2])
      break
    grange2 = granges[i + 1]
    g1 = grange[2]
    g2 = grange2[2]
    st1 = sample.loc[sample.Gene == g1, 'Strand'].iloc[0]
    st2 = sample.loc[sample.Gene == g2, 'Strand'].iloc[0]
    try:
      TPM1 = tpm.loc[tpm.gene == g1, 'tpm'].iloc[0]
      TPM2 = tpm.loc[tpm.gene == g2, 'tpm'].iloc[0]
    except KeyError:
      print('TPM not found in normalize_from_htseq output')
    try:
      RPKM1 = tpm.loc[tpm.gene == g1, 'rpkm'].iloc[0]
      RPKM2 = tpm.loc[tpm.gene == g2, 'rpkm'].iloc[0]
    except KeyError:
      print('RPKM not found in normalize_from_htseq output')
    cur_granges.append(grange)
    cur_flank = sample.Position.apply(flanking, args = [cur_granges, flank, len(sample)])
    cur_flank_both = sample.Position.apply(flanking, args = [cur_granges, read_size, len(sample), 'both'])

    sum_around_junction = sample.loc[cur_flank_both, 'coverage'].rolling(read_size).mean().sum()

    sample.loc[cur_flank, 'neigh_freq'] = round(prob_matrix.loc[utils.old_to_new(g1), utils.old_to_new(g2)]*100,1)
    sample.loc[cur_flank, 'junction_size'] = grange2[0] - grange[1]
    sample.loc[cur_flank, 'left_strand'] = st1
    sample.loc[cur_flank, 'right_strand'] = st2
    sample.loc[cur_flank, 'strand_switch'] = (st1 != st2)
    sample.loc[cur_flank, 'gpair'] = '_'.join((g1, g2))
    sample.loc[cur_flank, 'coexpression'] = num_sim(TPM1, TPM2)
    sample.loc[cur_flank, 'left_gene'] = g1
    sample.loc[cur_flank, 'right_gene'] = g2
    sample.loc[cur_flank, 'sum_around_junction'] = sum_around_junction
    try:
      sample.loc[cur_flank, 'left_tpm'] = TPM1
      sample.loc[cur_flank, 'right_tpm'] = TPM2
    except NameError: pass 
    try:
      sample.loc[cur_flank, 'left_rpkm'] = RPKM1
      sample.loc[cur_flank, 'right_rpkm'] = RPKM2
    except NameError: pass 

  return sample

def junction_only(sample : pd.DataFrame, remove_trna : bool) -> pd.DataFrame:
  """
  """
  tdict = {k: 'first' for k in ['org', 'coexpression', 'left_tpm', 'right_tpm', 'left_rpkm','right_rpkm', 'right_gene', 'left_gene', 'left_strand', 'right_strand', 'neigh_freq', 'strand_switch', 'Chromosome', 'Feature', 'junction_size','sum_around_junction']}
  tdict.update({'Position':['first', 'last', 'count'], 'RPM':['mean','median'], 'coverage':['mean','median'], 'z':['mean', 'median'], 'ends_ratio' : ['mean', 'max']})
  #Remove tRNA lines
  grouped = sample.groupby('gpair').agg(tdict).sort_values(by = ('Position', 'first'), ascending=True)
  cols = grouped.columns.to_list()
  cols = [(i[1] + '_').replace('first', '') + i[0] for i in cols]
  grouped.columns = cols


  return grouped

def rolling_around_junction(sample, granges, flank):
  both_flanks = sample.Position.apply(flanking, args = [granges, flank, len(sample), 'both'])
  filtered = sample.loc[both_flanks, :]
  filtered.loc[:, 'coverage'] = filtered.loc[:, 'coverage'].rolling(flank).mean()

def main(args : argparse.ArgumentParser) -> None:
  """
  The main function, receives a sample and flank window, saves to output a modified dataframe that only contains rows that are within the flank range of all genes
  Parameters
  ----------
  argpase.ArgumentParser args:

  annotated_pileup : str
    Path to the annotated csv created by process_pileup.py
  tpm_from_htseq : str
    Path to the csv created by normalize_from_htseq.py
  output : str
    Path to write the csv into
  flank : int
    Flanking junction size

  Notes
  -----
  This program does not return anything, it parses the arguements from argparse and creates a junctions csv file in the desired output location.
  """
  global PATH
  PATH = args.wd
  #Read samples file
  sample = pd.read_csv(args.annotated_pileup, index_col=0)
  #Read gene counts file
  counts = pd.read_csv(args.tpm_from_htseq, index_col=0)
  #Path to output
  output = args.output
  #Junction size
  flank = args.junc
  #True if no trna
  no_trna = args.no_trna
  #True if stranded mode (with only same-strand junctions)
  junc_avg = args.junc_avg

  try:
    org_df = pd.read_csv(os.path.join(PATH, 'final.csv'), index_col = 'organism')
  except FileNotFoundError:
    raise FileNotFoundError('final.csv must be in the same folder as this code!\n')
  try:
    phylum = org_df.loc[args.org, 'phylum']
  except KeyError:
    phylum = org_df.loc[args.org.replace('_',' '), 'phylum']
  #Grab the appropriate phylum probabillity matrix
  ph_prob = phylum + '_prob.csv'
  #Create granges of genes
  if junc_avg:
    sample_pos = sample.loc[sample.Strand == True, :]
    sample_neg = sample.loc[sample.Strand == False, :]
    granges_pos = genomic_ranges(sample_pos, no_trna)
    granges_neg = genomic_ranges(sample_neg, no_trna)
    granges = granges_neg + granges_pos
  else:
    granges = genomic_ranges(sample, no_trna)
  flanks = sample.Position.apply(flanking, args = [granges, flank, len(sample)])
  try:
    prob_matrix = pd.read_csv(os.path.join(PATH, 'prob_matrices','no_trna' if no_trna else 'yes_trna', ph_prob), index_col = 0)
    sample = neigh_freq_and_expression(sample, counts, prob_matrix, granges, flank, args.org, args.read)
  except FileNotFoundError:
    print(f'Neighborhood frequency could not be calculated! Please make sure prob_matrices are in PATH - {PATH}')
  if not args.coding:
    sample = sample.loc[flanks & ((sample['Length'] < 100) | (sample['Gene'].isnull())), :]
  else:
    sample = sample.loc[flanks, :]
  #sample.to_csv(output)
  
  junctions = junction_only(sample, no_trna)
  junctions['dataset'] = os.path.basename(os.path.dirname(os.path.dirname(args.annotated_pileup)))
  junctions = junctions.reset_index()
  if args.junc_reads:
    norm = pd.read_csv(args.junc_reads, index_col = 0)
    agged = []
    for gpair in junctions.gpair.unique():
      cur = norm.loc[norm.gene.str.contains(gpair), :].agg('mean' if junc_avg else sum)

      cur['gene'] = gpair
      agged.append(cur)
    norm = pd.DataFrame(agged)
    norm = norm.drop(columns = ['sample', 'kbp', 'rpk', 'rpm'])
    norm = norm.rename(mapper = {'gene':'gpair', 'glenghts' : 'window', 'tpm':'junc_tpm','rpkm':'junc_rpkm','counts': 'junc_counts'}, axis = 1)
    junctions = junctions.merge(norm, on = 'gpair', how = 'inner')
  junctions = junctions.dropna(subset = ['junc_tpm', 'junc_rpkm', 'junc_counts'])
  junctions.to_csv(output)

def path_check(path : str) -> str:
  """
  Check whether the path is valid
  """
  if os.path.exists(path):
    return path
  else:
    raise argparse.ArgumentTypeError(f'{path} is not a valid path!')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description = __doc__, prog = 'JunctionMaker')

  parser.add_argument('--annotated-pileup', help = 'Sample file path', type = path_check, metavar = 'annotated_pileup')
  parser.add_argument('--tpm-from-htseq', help = 'Modified htseq-counts output csv file with TPM normalization', type = path_check, metavar = 'tpm_file')
  parser.add_argument('--output', '-o', help = 'Output path')
  parser.add_argument('--junc', help = 'Junction size', metavar = 'int', type = int, default = 50)
  parser.add_argument('--coding', help = 'Add flag if coding junctions are wanted aswell', action = 'store_true', default = False)
  parser.add_argument('--version', action = 'version', version = '%(prog)s 1.0.1')
  parser.add_argument('--wd', help = 'Set the working dir to create directories in', default = os.getcwd(), type = path_check, metavar = 'WorkingDirectory')
  parser.add_argument('--org', help = 'Organism to work on', type = str, metavar='org_name')
  parser.add_argument('--no-trna', '-n', help = 'Ignore tRNA genes while making the junctions', default = False, action = 'store_true')
  parser.add_argument('--read', help = 'Read size for the sum around junction', default = 100, metavar = 'ReadSize')
  parser.add_argument('--junc-reads', help = 'Normalized reads for the flanking regions around junctions', default = None, metavar = 'JuncReadsCSV')
  parser.add_argument('--junc-avg', help = 'Adding this flag will cause the left and right side of a junction to be averaged instead of summed', default = False, action = 'store_true')
  args = parser.parse_args()
  main(args)