#%%
import pandas as pd
import numpy as np
import seaborn as sns
from ast import literal_eval
from Bio import SeqIO
from Bio import Entrez
from Bio import SeqFeature
import os, sys
#%%


if sys.platform != 'win32':
    slash = '/'
else:
    slash = '\\'

Entrez.email = 'shtolz@post.bgu.ac.il'
Entrez.API = '8f3dd96e08ee6735d6cc2abc9c94beaf3709'
PATH = os.getcwd()

def dir_checker(dirname):
  """
  check if a directory exists, if it does not, create it.

  Parameters
  ----------
  dirname : str
    Name of the directory to be checked
  """
  if not os.path.isdir(f'{PATH}{slash}{dirname}'):
    os.mkdir(f'{PATH}{slash}{dirname}')
def dup_assigner(gene, genelist):
  """
  """
  return gene + '*' * sum([gene in i for i in genelist])
def annotate(org_name, record):
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
  if os.path.exists(f'{PATH}{slash}annotated_df{slash}{org_name}.csv'):
    return pd.read_csv(f'{PATH}{slash}annotated_df{slash}{org_name}.csv', index_col = 0)
  ignore = ['CDS','source'] #Remove these annotations (focusing only on genes)
  genelist = []
  features = record.features
  seq = record.seq
  return_df = pd.DataFrame(columns = ['Position','Feature','Gene','Length']) #The organism specific df
  return_df['Position'] = list(range(1,len(seq)+1))
  for feature in features: #Iterate over organism features
    if feature.type not in ignore: 
      product = feature.qualifiers.get('product') #Isolate the gene product name
      if product: product = product[0]
      gene = feature.qualifiers.get('gene')
      if gene: gene = gene[0]
      for location in feature.location.parts:
        start = location.start
        end = location.end+1
        strand = location.strand
        if product:
          product = dup_assigner(product, genelist)
          genelist.append(product)
          return_df.loc[start:end, 'Gene'] = product
        elif gene:
          gene = dup_assigner(gene, genelist)
          genelist.append(gene)
          return_df.loc[start:end, 'Gene'] = gene
        return_df.loc[start:end, 'Feature'] = feature.type
        return_df.loc[start:end, 'Strand'] = True if strand == 1 else False
        return_df.loc[start:end, 'Length'] = end - start
  dir_checker('annotated_df')
  return_df.to_csv(f'{PATH}{slash}annotated_df{slash}{org_name}.csv')
  return return_df

def combine_with_sample(sample_df, annotated_df):
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
def z_score(cov, mean, std):
  """
  Calculate Z score for given coverage
  """
  return (cov - mean)/std

def overlap_ratio(sample):
  """
  Calculate the ratio between reads that overlap a given position to total reads

  Parameters
  ----------
  sample : pd.DataFrame
    The sample df.
  
  Returns
  -------
  sample : pd.DataFrame
    The sample df with the new columns.
  """

  sample['pos_overlap_index'] = 0
  sample['neg_overlap_index'] = 0
  sample['overlap_index'] = 0

  sample['pos_end_counts'] = sample.Location.str.count('(\$\,|\$[actg])|(\^\W\,|\$\W[actg])')
  sample['neg_end_counts'] = sample.Location.str.count('(\$\.|\$[ACTG])|(\^\W\.|\$\W[ACTG])')
  sample['end_counts'] = sample['pos_end_counts'] + sample['neg_end_counts']

  pos_quant = sample.loc[sample['pos_coverage'] > 0, 'pos_coverage'].quantile(0.15)
  neg_quant = sample.loc[sample['neg_coverage'] > 0, 'neg_coverage'].quantile(0.15)
  quant = sample.loc[sample['coverage'] > 0, 'coverage'].quantile(0.15)


  sample.loc[sample['pos_coverage'] > pos_quant, 'pos_ends_ratio'] = ((sample.Location.str.count('(\$\,|\$[actg])|(\^\W\,|\$\W[actg])')) / sample['pos_coverage'])
  sample.loc[sample['neg_coverage'] > neg_quant, 'neg_ends_ratio'] = ((sample.Location.str.count('(\$\.|\$[ACTG])|(\^\W\.|\$\W[ACTG])')) / sample['neg_coverage'])
  sample.loc[sample['coverage'] > quant, 'ends_ratio'] = (sample['end_counts'] / sample['coverage']
  )
  return sample

def strand_coverage(sample):
  """
  Measure per-strand coverage, fraction of reads that end in each position and normalized Z score.

  Parameters
  ----------
  sample : sample DF.

  Returns
  -------
  sample : sample DF modified with the new cols.
  """
  #Count exact amount of exact ref matches, indels and snps as coverage.

  #IMPORTANT: CURRENTLY SET AS: FORWARD = Negative, REVERSE = Positive
  sample['neg_coverage'] = sample.Location.str.count('\.|[ACTGN]|[\-\+][0-9]+[ACTGN]+') 
  sample['pos_coverage'] = sample.Location.str.count('\,|[actgn]|[\-\+][0-9]+[actgn]+')
  sample['neg_read_end'] = sample.Location.str.count('\$\.|\$[ACTGN]')
  sample['pos_read_end'] = sample.Location.str.count('\$\,|\$[actgn]')
  
  sample['read_end'] = sample['pos_read_end'] + sample['neg_read_end']
  
  sample['pos_read_end'] = sample.pos_read_end.fillna(0)
  sample['neg_read_end'] = sample.neg_read_end.fillna(0)
  sample['read_end'].fillna(0)

  sample.pos_coverage = sample.pos_coverage.fillna(0)
  sample.neg_coverage = sample.neg_coverage.fillna(0)
  sample.coverage = sample.coverage.fillna(0)

  pos_total = np.sum(sample['pos_read_end'])
  neg_total = np.sum(sample['neg_read_end'])
  sample['read_end'] = sample['read_end']/(pos_total + neg_total)
  print(f'Negative values: {neg_total}\nPositive values: {pos_total}\nTotal values: {pos_total + neg_total}')
  sample['pos_read_end'] = sample.pos_read_end/pos_total
  sample['neg_read_end'] = sample.neg_read_end/neg_total

  sample['RPM'] = sample['coverage']*1000000/(pos_total + neg_total)
  sample['pos_RPM'] = sample['pos_coverage']*1000000/(pos_total)
  sample['neg_RPM'] = sample['neg_coverage']*1000000/(neg_total)



  pos_mean = np.mean(sample['pos_coverage'])
  neg_mean = np.mean(sample['neg_coverage'])
  tot_mean = np.mean(sample['coverage'])
  
  pos_std = np.std(sample['pos_coverage'])
  neg_std = np.std(sample['neg_coverage'])
  tot_std = np.std(sample['coverage'])

  sample['pos_z'] = sample.pos_coverage.apply(z_score, args = [pos_mean, pos_std])
  sample['neg_z'] = sample.neg_coverage.apply(z_score, args = [neg_mean, neg_std])
  sample['z'] = sample.coverage.apply(z_score, args = [tot_mean, tot_std])

  sample = overlap_ratio(sample)
  return sample

def TPM_calc(sample):
  RPKS = []
  sample['TPM'] = 0
  for gene in sample.Gene.dropna().unique():
    cur_gene = sample.Gene == gene
    glength = sample.loc[cur_gene, 'Length'].max()
    kb = glength/1000
    RPK = sample.loc[cur_gene, 'coverage'].mean()/kb #TODO - this needs to be replaced with proper coverage measurement!
    RPKS.append(RPK)
  total_RPK = sum(RPKS)/1000000
  for i,gene in enumerate(sample.Gene.dropna().unique()):
    cur_gene = sample.Gene == gene
    sample.loc[cur_gene, 'TPM'] = RPKS[i]/total_RPK
  
  return sample



  

def record_check(org):
  """
  Check whether org genbank exists, if not download it.
  """
  try: org_df = pd.read_csv(os.path.join('data', 'final.csv'), index_col='organism')
  except FileNotFoundError: raise FileNotFoundError('final.csv must be in the data folder!')
  ID = org_df.loc[org, 'RefSeq']
  if type(ID) != str: ID = ID.iloc[0]
  filename = f'{PATH}{slash}genbank_DB{slash}{org}.gbk'
  if not os.path.isfile(filename):
    net_handle = Entrez.efetch(db = 'nucleotide', id = ID, rettype = 'gb', retmode = 'text')
    out_handle = open(filename, 'w')
    out_handle.write(net_handle.read())
    net_handle.close()
    out_handle.close()

def main(path_to_accs, org_name):
  """
  Combine all above functions to create an annotated .csv based on .pileup file path and org_name
  """
  
  dir_checker('genbank_DB')
  dir_checker('annotated_pileup')
  dir_checker(f'annotated_pileup{slash}{org_name}')
  filename = f'{PATH}{slash}genbank_DB{slash}{org_name}.gbk'
  record_check(org_name)

  record = SeqIO.read(filename, 'genbank')
  sample = pd.read_csv(path_to_accs, delimiter = '\t',names = ['Chromosome','Position','Base','coverage','Location','Quallity'], index_col = None, usecols = list(range(0,6)))
  sample = strand_coverage(sample)
  sample = sample.drop(columns = ['Location','Quallity'])
  annotated = annotate(org_name, record)
  final = combine_with_sample(sample, annotated)
  final[['pos_coverage','neg_coverage']] = final[['pos_coverage','neg_coverage']].fillna(0)
  final = TPM_calc(final.copy())

  accs = path_to_accs.split("\\")[-1].replace(".pileup","")
  final.to_csv(f'{PATH}{slash}annotated_pileup{slash}{org_name}{slash}{accs}.csv')
    
def get_length(org):
  """
  Return the mtDNA sequence length of an organism
  """
  record_check(org)
  filename = f'{PATH}{slash}genbank_DB{slash}{org}.gbk'
  record = SeqIO.read(filename, 'genbank')
  return len(record.seq)


    

#%%
