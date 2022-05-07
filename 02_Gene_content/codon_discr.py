#%%
import numpy as np
import re,os, importlib
from Bio import Entrez
from Bio import SeqIO
from sys import platform
from textwrap import wrap
import utils
from functools import wraps
from consts import GENE_TO_ABB, SYNS, CODON_TABLES
import CAI as cai 
importlib.reload(utils)

if platform != 'win32':
    DASH = '/'
else:
    DASH = '\\'
path = os.getcwd()

"""
Function steps:
1. Retrieve sequence based on ID.
2. Retrieve CDR based on sequence.
3. Count codons of each gene within CDR dictionary.
4. For each organism, compare the tRNA availibillity to existing codons.
"""

def logger(func):
    """
    Simple logging decorator:
    Prior function call: log function name and arguements.
    Post funciton call: log returned variable values.
    """
    if not os.path.exists(f'{path}{DASH}logs'):
        os.mkdir(f'{path}{DASH}logs')
    @wraps(func)
    def wrapper(*args, **kwargs):
        with open(f'{path}{DASH}logs{DASH}codon_disc_dump.txt','a') as log:
            log.write(f'Running function {func.__name__} with parameters:\t {args} and keywords: {kwargs}\n\n\n')
            returned = func(*args,**kwargs)
            log.write(f'{func.__name__} returned:\t {returned}\n\n\n')
        return returned
    return wrapper

@logger
def retriever(row):
    """
    Isolate organism DNA sequence as a string based on recieved ID and return a dict of codon start locations

    Parameters
    ----------
    row : pd.Series
        A row of the database
    
    Returns
    -------
    Seq : str
        the organism's entire DNA sequence
    Codon_start : dict
        dict of CDS gene names and their codon start location
    """
    if not os.path.exists(f'{path}{DASH}genbank_DB'):
        os.mkdir(f'{path}{DASH}genbank_DB')
    genbank = f'{path}{DASH}genbank_DB{DASH}{row["RefSeq"]}.gbk' #based on row's RefSeq, create path to gbk file.
    codon_starts={}
    try:
        record = SeqIO.read(genbank, 'genbank') # see if the genbank ID exists
    except (FileNotFoundError, ValueError): # If it doesnt exist, retrieve from NCBI and open record
        net_handle = Entrez.efetch(db = 'nucleotide', id = row['RefSeq'], rettype = 'gb', retmode = 'text')
        out_handle = open(genbank, 'w')
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print('Saved {row["RefSeq"]}\n')
        record = SeqIO.read(genbank, 'genbank')
    for feature in record.features: # Iterate over features to grab gene names and codon_start locations
        if feature.type == 'CDS': #If the feature is a coding sequence (CDS)
            try: # TODO : change to a simple type check
                gene = feature.qualifiers.get('product')[0] # Take product (as a proxy for gene name)
            except TypeError: #If the gene does not exist (gene would be a None object raising typeerror if i index it)
                continue
            gene = utils.Which_tRNA(gene) # Use Create_database.py's Which_tRNA function to unify gene names
            if not gene: continue # If Which tRNA returned nothing
            if gene in codon_starts.keys(): gene+='*'
            try:
                codon_starts[gene] = feature.qualifiers.get('codon_start')[0] #grab the codon start site (where count starts from)
            except TypeError:
                print(f'Issue on {row["organism"]}')
                codon_starts[gene] = 1
    return record.seq, codon_starts #Returns the full mt sequence and the codon start site

@logger
def CDS(row, seq, codon_starts):
    """
    Based on an organism's org_df row, sequence and codon start locations,
    split an organism's coding sequences into a dictionary, taking into account incomplete stop codons.

    Parameters
    ----------
    row : pd.Series
        A row of the org_db database
    seq : pd.Str
        A string of the entire mtDNA sequence of the organism
    codon_starts :  list
        A list of all CDS translation start locations
    
    Returns
    -------
    genes : dict
        A dictionary of Gene_Name:sequence for every CDS in organism.
    """
    genes = {} #initialize gene_name:sequence dict.
    gorder = row['Gene_order']
    gloc = row['Gene_locations']
    notes = row['Notes']
    gcode = int(row['genetic_code'])
    if (len(gorder) * len(gloc) * len(notes)) != len(gorder) ** 3:
        raise ValueError(f"""Notes,
     Gene_order and Gene_locations do not have the exact same length!
     \ngene order length is: {len(gorder)} 
     \ngene location length is: {len(gloc)}
     \nnotes length is: {len(notes)}""")
    for gene,loc,note in zip(gorder, gloc, notes): # iterate over gene names, start:end:strand and notes from the org_db
        gene = gene.replace('-','')
        if gene in codon_starts.keys(): #filter for protein coding
            if len(note)>0: #ignore any notes that contain frameshift in notes
                if 'frameshift' in ''.join(note): continue
            try: start,end,strand = [int(i.replace('>','').replace('<','')) for i in loc.split(':')]
            except ValueError: continue
            codon_start = int(codon_starts[gene]) - 1
            start += codon_start
            if strand == 1:
                geneseq = seq[start:end]
            else:
                geneseq = seq[start:end].reverse_complement()
            stops = [i for i in CODON_TABLES[gcode].keys() if CODON_TABLES[gcode][i] == 'STOP']
            if geneseq[-3:] not in stops:
                while geneseq[-3:] not in stops:  # for i in range(1): if geneseq[-3] not in stops: geneseq+="A"
                    geneseq+='A'
                    if len(geneseq) > 10000: #can remove
                        print(f'It appears this organisms codon completion is done by something other than adenylylation\nID = {row.RefSeq}')
                        break
            if len(geneseq) > 10000: continue
            genes[gene] = str(geneseq)
    return genes 

def CDS2(row):
    """
    Based on an organism's org_df row, return a dictionary of gene names and their corresponding DNA sequence.

    Parameters
    ----------
    row : pd.Series
        A row of the org_db database
    
    Returns
    -------
    genes : dict
        A dictionary of Gene_Name:sequence for every CDS in organism.
    """
    genes = {}
    if not os.path.exists(f'{path}{DASH}genbank_DB'):
        os.mkdir(f'{path}{DASH}genbank_DB')
    genbank = f'{path}{DASH}genbank_DB{DASH}{row["RefSeq"]}.gbk' #based on row's RefSeq, create path to gbk file.
    codon_starts={}
    try:
        record = SeqIO.read(genbank, 'genbank') # see if the genbank ID exists
    except (FileNotFoundError, ValueError): # If it doesnt exist, retrieve from NCBI and open record
        net_handle = Entrez.efetch(db = 'nucleotide', id = row['RefSeq'], rettype = 'gb', retmode = 'text')
        out_handle = open(genbank, 'w')
        out_handle.write(net_handle.read())
        out_handle.close()
        net_handle.close()
        print('Saved {row["RefSeq"]}\n')
        record = SeqIO.read(genbank, 'genbank')
    for feature in record.features: # Iterate over features to grab gene names and codon_start locations
        if feature.type == 'CDS': #If the feature is a coding sequence (CDS)
            try:
                gene = feature.qualifiers.get('product')[0] # Take product (as a proxy for gene name)
            except TypeError: #If the gene does not exist (gene would be a None object raising typeerror if i index it)
                continue
            gene = utils.Which_tRNA(gene) # Use Create_database.py's Which_tRNA function to unify gene names
            if not gene: continue # If Which tRNA returned nothing
            if gene in codon_starts.keys(): gene+='*'
            genes[gene] = str(feature.extract(record.seq))
    return genes

@logger
def codon_counter_func(row, genes):
    """
    Based on org_df row and genes dictionary, count occurences of each codon

    Parameters
    ----------
    row : pd.Series
        A single row of org_df dataframe
    genes : dict
        A dictionary in this form: CDS_gene_name:sequence
    
    Returns
    -------
    codon_counter : dict
        A dictionary in this form: codon:amount_of_times_this_codon_appears
    """ 
    gcode = int(row['genetic_code'])
    codon_counter = {k:[v,0] for k,v in CODON_TABLES[gcode].items()}
    for name,gene in genes.items():
        for codon in wrap(gene,3):
            if re.search(pattern = '[^ATGC]', string = codon):
                continue
            if len(codon) < 3:
                print(f"Frameshift error! ID: {row.RefSeq} on gene: {name} of length: {len(gene)}")  
            try:
                codon_counter[codon][1] += 1
            except KeyError:
                #print(f'Issue on {row.name}\nsequence is not a codon! {codon}')
                continue
    return {k:v for k,v in codon_counter.items() if CODON_TABLES[gcode][k] != 'STOP'}

def check_if_recognized(row, wanted_tRNA):
    """ A sub-function of compare_codon which only checks which codons are recognized by the organism's mtDNA codons"""
    dicr = {k:np.nan for k in wanted_tRNA.values()}
    tRNA_recognized = {k:[] for k in wanted_tRNA.values()}

    for gene in row['Gene_order']: #re.search(r'\(\w{3}\)')
        gene = gene.replace('-','')
        gene = gene.replace('*','')
        if gene[0:4] in wanted_tRNA.values():
            codon = gene[5:8]
            if len(codon) < 3: continue
            codon = codon.replace('U','T')
            tRNA_recognized[gene[0:4]]+=[codon[0:2] + i for i in SYNS[codon[2]]]

    print(f'{tRNA_recognized=}')
    return tRNA_recognized

def compare_codon(row, codon_counter, wanted_tRNA = {'S':'trnS', 'L':'trnL'}, only_trn_recognized = False):
    """
    based on a row of org_db and a dictionary of codon counts, compare existing codons
    within the organim's CDS to its tRNA repertoire.

    Parameters
    ----------
    row : df.Series
        A single row of org_db dataframe
    codon_counter : dict
        A dicitonary in this format: codon:[translation(AA), count_of_appearances_within_CDS]
    
    Returns
    -------
    dicr : dict
        A dictionary in this format: amino_acid: True if complete match between tRNA repertoire and existing codons within CDS
    """
    
    dicr = {k:np.nan for k in wanted_tRNA.values()}
    tRNA_recognized = {k:[] for k in wanted_tRNA.values()}
    existing_codons = {k:[] for k in wanted_tRNA.values()}
    codon_counter = {k:v for k,v in codon_counter.items() if v[0] in wanted_tRNA.keys()}
    for tRNA in wanted_tRNA.values():
        existing_codons[tRNA] = [i[0] for i in codon_counter.items()\
                                    if wanted_tRNA[i[1][0]] == tRNA]
    for gene in row['Gene_order']: #re.search(r'\(\w{3}\)')
        gene = gene.replace('-','')
        gene = gene.replace('*','')
        if gene[0:4] in wanted_tRNA.values():
            codon = gene[5:8]
            if len(codon) < 3: continue
            codon = codon.replace('U','T')
            tRNA_recognized[gene[0:4]]+=[codon[0:2] + i for i in SYNS[codon[2]]]
    if only_trn_recognized:
        return tRNA_recognized
    for tRNA in wanted_tRNA.values():
        if len(tRNA_recognized[tRNA])<1: continue 
        dicr[tRNA] = set(tRNA_recognized[tRNA]) == set(existing_codons[tRNA])
    return dicr

def main_codon_dis(row):
    """
    Combine all above functions, recieve a single row of org_db dataframe, return a dis_dict 
    """
    seq,codon_starts = retriever(row)
    genes = CDS(row, seq, codon_starts)
    codon_counter = codon_counter_func(row, genes)
    if codon_counter != None:
        dis_dict = compare_codon(row, codon_counter)
        dis_dict.update({'organism':row.name})
        return dis_dict
    else: return {'organism':row.name, 'trnS':np.nan, 'trnL':np.nan}

def main_codon_dict(row):
    """
    """
    seq,codon_starts = retriever(row)
    genes = CDS(row, seq, codon_starts)
    codon_counter = codon_counter_func(row, genes)
    concordance = compare_codon(row, codon_counter)
    if concordance != None:
        return codon_counter

### RSCU CALCULATIONS BLOCK ###

def calculate_relative_adapt(row):
    """
    Calculate RSCU for each organism in a single row of org_db dataframe
    return a dictionary of relative adaptation values for each codon within the organism's CDS
    """
    seq,codon_starts = retriever(row)
    cds = CDS(row, seq, codon_starts)
    codon_counter = codon_counter_func(row, cds)
    if codon_counter == None:
        dicr = []
    else:
        dicr = compare_codon(row, codon_counter, only_trn_recognized = True)
    if dicr == None:
        dicr = []
    gcode = int(row['genetic_code'])
    try:
        raw_rscu = cai.RSCU(list(cds.values()), genetic_code = gcode)
        rscu = cai.relative_adaptiveness(list(cds.values()), genetic_code = gcode)

    except (ValueError, TypeError) as exc:
        print(f'Issue on {row.organism}!\n{exc}')
        return None
    return_dict = {
        'codon' : [i for i in rscu.keys()],
        'RSCU_adp' : [i for i in rscu.values()],
        'RSCU_raw' : [raw_rscu[i] for i in rscu],
        'codon_count': [codon_counter[i] for i in rscu.keys()],
        'organism' : [row.organism for _ in rscu.keys()],
        'phylum' : [row.phylum for _ in rscu.keys()],
        'class' : [row['class'] for _ in rscu.keys()],
        'AA' : [CODON_TABLES[gcode][i] for i in rscu.keys()],
        'agreement' : [i in dicr[GENE_TO_ABB[CODON_TABLES[gcode][i]]] if GENE_TO_ABB[CODON_TABLES[gcode][i]] in dicr.keys() else '' for i in rscu.keys()]
        }
    return return_dict


def calculate_cai(row, reference = 'all'):
    """
    For a given row, calculate a per-gene codon adaptation index
    """
    seq, codon_starts = retriever(row)
    cds = CDS(row, seq, codon_starts)
    gcode = int(row['genetic_code'])
    return_dict = {
        'cai' : [],
        'gene' : [],
        'organism' : [],
        'phylum' : [],
        'class' : []
    }
    for gene in cds.keys():
        try:
            # TODO(Noam): Think whether it is correct to use all genes as reference.
            cindex = cai.CAI(cds[gene], genetic_code = gcode)
        except (ValueError, TypeError, KeyError) as err:
            print(f'Could not calculate CAI for {gene} in {row.organism}!\n{err}')
            continue
        return_dict['cai'] += [cindex]
        return_dict['gene'] += [gene]
        return_dict['organism'] += [row.organism]
        return_dict['phylum'] += [row.phylum]
        return_dict['class'] += [row['class']]
    return return_dict

    