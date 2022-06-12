import re
old_to_new_dict = {
    'tRNA-Ala':'trnA','tRNA-Cys':'trnC','tRNA-Asp':'trnD','tRNA-Glu':'trnE',
    'tRNA-Phe':'trnF','tRNA-Gly':'trnG','tRNA-His':'trnH','tRNA-Ile':'trnI',
    'tRNA-Lys':'trnK','tRNA-Leu':'trnL','tRNA-Met':'trnM','tRNA-Asn':'trnN',
    'tRNA-Pro':'trnP','tRNA-Gln':'trnQ','tRNA-Arg':'trnR','tRNA-Ser':'trnS',
    'tRNA-Thr':'trnT','tRNA-Val':'trnV','tRNA-Trp':'trnW','tRNA-Tyr':'trnY',
    'ATP6':'atp6', 'ATP8':'atp8', 'ND1': 'nad1', 'ND2': 'nad2', 'ND3': 'nad3', 'ND4': 'nad4',
    'ND4L':'nad4L', 'ND5':'nad5','ND6':'nad6','COX1': 'cox1',
    'COX2': 'cox2', 'COX3': 'cox3', '12S ribosomal RNA': 'rrnS',
    '16S ribosomal RNA': 'rrnL' , 'CYTB':'cob', 'l-rRNA': 'rrnL', 's-rRNA' : 'rrnS', 'lrRNA':'rrnL', 'srRNA' : 'rrnS', 'RNR1':'rrnS', 'RNR2':'rrnL', '12S rRNA':'rrnS', '16S rRNA':'rrnL', 'rnl': 'rrnL', 'rns': 'rrnS'}

new = list(old_to_new_dict.values())
new_to_old_dict = {v:k for k,v in old_to_new_dict.items()}
def old_to_new(gene):
  """
  Replace old format "tRNA-Ala" with shorter new format "trnA"
  """
  if 'trn' in gene.lower():
    gene = re.sub(r'\d', '', gene)
  try:
    gene = old_to_new_dict[gene]
  except KeyError:
    for i in old_to_new_dict.values():
      if i.lower() == gene.lower():
        return i
    print(f'{gene} not found!')
  return gene 

def new_to_old(gene):
  """
  Replace new format "trnA" with longer old format "tRNA-Ala"
  """
  try:
    gene = new_to_old(gene)[gene]
  except KeyError:
    print(f'{gene} not found!')
  return gene

