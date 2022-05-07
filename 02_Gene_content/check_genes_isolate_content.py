TRNA_DICT = {
    'tRNA-Ala':'trnA','tRNA-Cys':'trnC','tRNA-Asp':'trnD','tRNA-Glu':'trnE',
    'tRNA-Phe':'trnF','tRNA-Gly':'trnG','tRNA-His':'trnH','tRNA-Ile':'trnI',
    'tRNA-Lys':'trnK','tRNA-Leu':'trnL','tRNA-Met':'trnM','tRNA-Asn':'trnN',
    'tRNA-Pro':'trnP','tRNA-Gln':'trnQ','tRNA-Arg':'trnR','tRNA-Ser':'trnS',
    'tRNA-Thr':'trnT','tRNA-Val':'trnV','tRNA-Trp':'trnW','tRNA-Tyr':'trnY'}
TRNA_SET = set([i for i in TRNA_DICT.values()])


def full_trna_package(gorder, return_missing = False):
    """
    Recieve gene order and return True if it contains trnas for all 20 amino acids.
    """
    gorder = set([i.replace('-', '').replace('*', '')[0:4] for i in gorder if 'trn' in i])
    if return_missing:
        return list(TRNA_SET.difference(gorder))
    return TRNA_SET.difference(gorder) == set()

def check_gene(gorder, gene):
    """
    Recieve a gene order and return True if it contains gene
    """
    gorder = [i.replace('-', '').replace('*', '')[0:4] for i in gorder]
    return gene in gorder 

def check_gene_length(row, gene):
    """
    Recieve row from org_df and a gene, return its length
    """
    if type(row.Gene_locations) != list:
        raise ValueError('Please run literal_eval on dataframe before running this function!') 
    else:
        lengths = row.Gene_locations
        gorder = row.Gene_order
        try: ind = gorder.index(gene)
        except ValueError: return 0
        formula = lengths[ind].split(':')
        if len(formula) < 3: return 0
        glength = int(formula[1].replace('<', '').replace('>', '')) - int(formula[0].replace('<', '').replace('>', ''))
        return glength


def to_set(gorder, keep_dups = False):
    """
    """
    gorder = [i.replace('-', '').replace('5_', '').replace('_3', '') for i in gorder]
    gorder = [i.replace('*', '') if not keep_dups else i for i in gorder]
    gset = frozenset([(i[0:4] if 'trn' in i else i) + ('*' if '*' in i else '') for i in gorder])
    return gset

