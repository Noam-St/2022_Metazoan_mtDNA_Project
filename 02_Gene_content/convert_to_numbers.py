#%%
from ast import literal_eval
import os, re
import utils as cd
import importlib
importlib.reload(cd)
path = os.getcwd()
def winPath2ubuPath(winpath):
    # d,p = os.path.splitdrive(winpath) # NG only works on windows!
    try:
        d,p = winpath.split(':')
        ubupath = '/mnt/'+d.lower()+p.replace('\\','/')   
        print (ubupath)
        return ubupath
    except ValueError:
        print('Unpack error! path recieved is most likely already unix')
        return winpath

REP_DICT={
          '1':'tRNA-Phe',
          '17':'tRNA-His',
          '14':'tRNA-Lys',
          '16':'tRNA-Arg',
          '13':'tRNA-Asp',
          '19':'tRNA-Glu',
          '12':'tRNA-Ser',
          '20':'tRNA-Thr',
          '9':'tRNA-Asn',
          '5':'tRNA-Gln',
          '8':'tRNA-Ala',
          '2':'tRNA-Val',
          '18':'tRNA-Leu',
          '11':'tRNA-Tyr',
          '7':'tRNA-Trp',
          '21':'tRNA-Pro',
          '15':'tRNA-Gly',
          '10':'tRNA-Cys',
          '6':'tRNA-Met',
          '4':'tRNA-Ile',
          '22':'12S ribosomal RNA',
          '23':'16S ribosomal RNA',
          '24':'ATP6',
          '25':'ATP8',
          '26':'COX1',
          '27':'COX2',
          '28':'COX3',
          '29':'CYTB',
          '30':'ND1',
          '31':'ND2',
          '32':'ND3',
          '33':'ND4',
          '34':'ND4L',
          '35':'ND5',
          '36':'ND6',
          '37':'ATP9',
          '38':'mutS',
          '39':'heg',
          '40':'secY',
          '41':'Reph',
          '42':'ORFX',
          '43':'RNAX'}
for k,v in REP_DICT.items():
    REP_DICT[k] = cd.Which_tRNA(v)

def remove_codon(gene_order):
    """ Remove codon annotation from gene name for Distance matrix creation"""
    gene_list2 = []
    gene_list = []
    for gene in gene_order:
        gene = str(gene)
        gene=gene.replace('*','')
        gene = re.sub(pattern='\(\w{3}\)',repl = '',string = gene)
        gene_list2.append(gene.replace('-',''))
        counter = gene_list2.count(gene.replace('-',''))
        gene = gene+'*'*(counter-1)
        gene_list.append(gene)
    return gene_list

def replace_func(to_rep,by_dict):
    """ recives a str and dict, replaces str with any matching values in dict and returns dict reversed to replace string annotations
    with numerical annotations"""
    counter=0
    bools=False
    if '*' in to_rep: #checks if this gene is a duplicate marked by *
        counter=to_rep.count('*')
    to_rep=to_rep.replace('*','')
    if '-' in to_rep[0]: #checks if the gene is in the complementary strand
        bools=True
        to_rep=to_rep.replace('-','',1)
    for k,v in by_dict.items(): #iterates over the replacement dictionary and looks for matching gene strings
        if to_rep==v:
            to_rep=k
            break
    to_rep=to_rep+'*'*counter
    if bools: to_rep='-'+to_rep #If gene is in complementry, if the answer is yes 
    return to_rep

def replace_func_on_list(listy):
    """ Apply replace function on each item in list"""
    if type(listy)==str:
        listy=literal_eval(listy)
    if type(listy)==list:
        for i in range(len(listy)):
            listy[i]=replace_func(listy[i],REP_DICT)
    return listy
