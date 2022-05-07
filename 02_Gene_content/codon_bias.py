# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 16:21:44 2020

@author: noam1
"""

#%%
from Bio.SeqUtils import CodonUsage
import numpy as np
import Bio
from Bio import Entrez
from Bio import SeqIO
import pandas as pd
import re,os
import matplotlib.pyplot as plt
import seaborn as sns
from ast import literal_eval
#%%
codon_tables={'Chordata':{ "CYS": ["TGT", "TGC"], 
                           "ASP": ["GAT", "GAC"],
                           "SER": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"],
                           "GLN": ["CAA", "CAG"],
                           "MET": ["ATG", "ATA"],
                           "ASN": ["AAC", "AAT"], 
                           "PRO": ["CCT", "CCG", "CCA", "CCC"],
                           "LYS": ["AAG", "AAA"],
                           "STOP": ["TAG", "TAA", "AGA", "AGG"], 
                           "THR": ["ACC", "ACA", "ACG", "ACT"],
                           "PHE": ["TTT", "TTC"],
                           "ALA": ["GCA", "GCC", "GCG", "GCT"],
                           "GLY": ["GGT", "GGG", "GGA", "GGC"],
                           "ILE": ["ATC", "ATT"],
                           "LEU": ["TTA", "TTG", "CTC", "CTT","CTG", "CTA"],
                           "HIS": ["CAT", "CAC"],
                           "ARG": ["CGA", "CGC", "CGG", "CGT" ],
                           "TRP": ["TGG", "TGA"], 
                           "VAL": ["GTA", "GTC", "GTG", "GTT"],
                           "GLU": ["GAG", "GAA"], 
                           "TYR": ["TAT", "TAC"] },
                'Others':{ "CYS": ["TGT", "TGC"],
                           "ASP": ["GAT", "GAC"],
                           "SER": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT", "AGA", "AGG"],
                           "GLN": ["CAA", "CAG"],
                           "MET": ["ATG", "ATA"],
                           "ASN": ["AAC", "AAT"], 
                           "PRO": ["CCT", "CCG", "CCA", "CCC"],
                           "LYS": ["AAG", "AAA"],
                           "STOP": ["TAG", "TAA"], 
                           "THR": ["ACC", "ACA", "ACG", "ACT"],
                           "PHE": ["TTT", "TTC"],
                           "ALA": ["GCA", "GCC", "GCG", "GCT"],
                           "GLY": ["GGT", "GGG", "GGA", "GGC"],
                           "ILE": ["ATC", "ATT"],
                           "LEU": ["TTA", "TTG", "CTC", "CTT","CTG", "CTA"],
                           "HIS": ["CAT", "CAC"],
                           "ARG": ["CGA", "CGC", "CGG", "CGT" ],
                           "TRP": ["TGG", "TGA"], 
                           "VAL": ["GTA", "GTC", "GTG", "GTT"],
                           "GLU": ["GAG", "GAA"], 
                           "TYR": ["TAT", "TAC"]},
    'starfish_flatworms':{ "CYS": ["TGT", "TGC"], 
                           "ASP": ["GAT", "GAC"],
                           "SER": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT", "AGA", "AGG"],
                           "GLN": ["CAA", "CAG"],
                           "MET": ["ATG", "ATA"],
                           "ASN": ["AAC", "AAT", "AAA"], 
                           "PRO": ["CCT", "CCG", "CCA", "CCC"],
                           "LYS": ["AAG"],
                           "STOP": ["TAG", "TAA"], 
                           "THR": ["ACC", "ACA", "ACG", "ACT"],
                           "PHE": ["TTT", "TTC"],
                           "ALA": ["GCA", "GCC", "GCG", "GCT"],
                           "GLY": ["GGT", "GGG", "GGA", "GGC"],
                           "ILE": ["ATC", "ATT"],
                           "LEU": ["TTA", "TTG", "CTC", "CTT","CTG", "CTA"],
                           "HIS": ["CAT", "CAC"],
                           "ARG": ["CGA", "CGC", "CGG", "CGT" ],
                           "TRP": ["TGG", "TGA"], 
                           "VAL": ["GTA", "GTC", "GTG", "GTT"],
                           "GLU": ["GAG", "GAA"], 
                           "TYR": ["TAT", "TAC"] }}
Entrez.email = "shtolz@post.bgu.ac.il"
Entrez.api_key="8f3dd96e08ee6735d6cc2abc9c94beaf3709"

#%%handle = Entrez.efetch(db='nucleotide',id=Id,rettype='gb',retmode='gbk')
def handle_fetcher(Id):
    """"Recieves ID, fetches records from Entrez and saves a record txt file that can be read as gb
    returns full filepath"""
    with Entrez.efetch(db='nucleotide',id=Id,rettype='gb',retmode='gbk') as handle: #Fetch record.
        fullpath = f'C:\\Users\\noam1\\Desktop\\Mitochondrial Gene comparision project\\Redoing the work with Bio\\Records\\{Id}.gbk'
        file = open(f'{fullpath}','w') #Open record file
        file.write(handle.read()) #write in the handle
        file.close()
        return fullpath #returns the path to handle file that was created
        print(f'Saved {Id}')



#for feature in record.features:
#    if feature.type == 'gene':
#        
#        print(feature.qualifiers.get('gene'))
def unite_seq(record,seq):
    """" Recieves record obj and sequence for x org, returns the entire coding sequence of that org with
    correct strands"""
    try:
        logfile = open(r'C:\Users\noam1\Desktop\Mitochondrial Gene comparision project\Redoing the work with Bio\Logs\unite_seq_log.txt','w')
        finalseq = {}
        for feature in record.features: #iterate over features
            if (feature.type=='gene'): #Filter only genes
                if feature.qualifiers.get('gene')==None : continue
                genecode = str(feature.qualifiers.get('gene')[0])
                if not ('TR' in genecode or 'RNR' in genecode or 'tr' in genecode): #only protein coding genes (not rRNA or tRNA)
                    loc = feature.location
                    logfile.write(f'This is the start of {genecode}: {loc.start}\n')
                    logfile.write(f'This is the end of {genecode}: {loc.end}\n')
                    if loc.strand == 1:
                        finalseq[genecode]=str(seq[loc.start:loc.end])  # If + strand
                    else:
                        finalseq[genecode]=str(seq[loc.start:loc.end].reverse_complement())# if - strand take reverse complement
                    logfile.write(f'Length of {genecode}: {len(finalseq[genecode])}\n')
        return finalseq
    except TypeError as err:
        print(f'{err}: Stuck on {feature}')
def count_codons(seq,Id):
    """" Recieves dict of all coding sequences annotated, and current org Id, counts codons and 
    returns a dict with the observed number of each codon"""
    try:
        
        table = full_taxa[full_taxa.Replicons==Id]['phylum'].to_string( index = False, header = False) # Define only the phylum name as table
        if table == 'Echinodermata' or table == 'Platyhelminthes': #Detect the correct table
            table = 'starfish_flatworms'
        elif table!='Chordata':
            table = 'Others'
        codon_counts={k:{i:0 for i in codon_tables[table][k]} for k in codon_tables[table].keys() } #Define a dict to fill, with A.A as keys and empty dict values.
        for k in codon_tables[table]:
            for gene in seq.values():
                codons = (gene[i:i+3] for i in range(0,len(gene),3))
                for codon in codons:
                    if codon in codon_tables[table][k]: codon_counts[k][codon] += 1
        return codon_counts
    except KeyError as err:
        print(f'{err}\n Please make sure the codon_tables dict contains all amino acids as keys and codons as values')
    except TypeError as err:
        print(f'{err}\n Please enter the following:\n seq : A dictionary with gene names as keys and sequence as values\n id : the replicon ID of the current organism')

def calculate_rscu(codon_counts,Id):
    """" Recieves dict of codon counts and current organism Id, calculates RSCU values for each codon and
    returns a nested dict, A.A as keys with values being another dict with codons as keys and rscu as values"""
    logfile = open(r'C:\Users\noam1\Desktop\Mitochondrial Gene comparision project\Redoing the work with Bio\Logs\calculate_rscu.txt','w')
    logfile.write('Current Organism:' + full_taxa[full_taxa.Replicons==Id].index[0] + ' phylum: ' + full_taxa[full_taxa.Replicons==Id]['phylum'].to_string( index = False, header = False))
    codon_rscu = {k:{i:0 for i in codon_counts[k]} for k in codon_counts.keys()}
    for k,v in codon_counts.items():
        ka = len(v) #Amount of synonymous codons for this A.A
        Oa = sum(v.values()) #Sum of observed codons for this specific A.A
        if Oa==0:
            continue
        logfile.write(f'\nCurrent A.A: {k}\t amount of synonym codons: {ka}\t Sum of all codons {Oa}\n')
        for i,Oac in v.items():
            codon_rscu[k][i] = Oac/((1/ka)*Oa) #RSCU calculation for each codon
            logfile.write(f'Current Codon: {i}\t amount of codons found: {Oac}\t RSCU: {codon_rscu[k][i]}\n')
    return codon_rscu

def main(Id):
    """" Unites all functions together, recieves Id and returns codon rscu dict"""
    filepath=f'C:\\Users\\noam1\\Desktop\\Mitochondrial Gene comparision project\\Redoing the work with Bio\\Records\\{Id}.gbk'
    if not os.path.exists(filepath): # If the file doesnt exist yet, creates it using handle fetcher
        filepath = handle_fetcher(Id)
        print('handle_fetcher was invoked')
    record = SeqIO.read(filepath, "genbank") # saves the record for the current organism
    sequence = record.seq #saves the sequence for the current organism
    sequence = unite_seq(record,sequence) #unites all coding sequences for codon bias calcs
    if type(sequence) == 'NoneType':  # checks that the unite_seq worked 
        print('No genes on this organism')
        return None
    codon_counts = count_codons(sequence, Id)  # counts how many codons of each type are in the CDS
    codon_rscu = calculate_rscu(codon_counts,Id) #calculates RSCU
    return codon_rscu

if __name__ == '__main__':
    x = handle_fetcher(full_taxa.loc['Coatitermes_kartaboensis','Replicons'])
    record = SeqIO.read(x,'genbank')
    sequence = record.seq
    y = unite_seq(record,sequence)
    #%%
    full_taxa['RSCU']=''
    full_taxa['RSCU']=full_taxa.Replicons.apply(main)
    #%%
    full_taxa.to_csv(r'C:\Users\noam1\Desktop\Mitochondrial Gene comparision project\Redoing the work with Bio\Codons\organisms_with_codon_usage.csv')

    #%% plotting the data
    plot_dict = {}
    phyl = 'Echinodermata'
    for aa in codon_tables['starfish_flatworms']:
        plot_dict[aa]={}
        for codon in codon_tables['starfish_flatworms'][aa]:
            plot_dict[aa][codon]=[]
            for i in full_taxa[full_taxa.phylum == phyl].index:
                plot_dict[aa][codon].append(full_taxa.RSCU[i][aa][codon])
    #%%
    for aa in codon_tables['starfish_flatworms']:
        temp = pd.DataFrame.from_dict(data = plot_dict[aa], orient = 'columns')
        plt.style.use('seaborn-darkgrid')
        fig1, ax1 = plt.subplots(figsize=(15,10))
        ax1 = sns.boxplot(data=temp, orient = 'v')
        ax1.set_xlabel('Codon')
        ax1.set_ylabel('RSCU')
        ax1.set_title(f'{phyl} codon usage for {aa}')
        ax1.set_yticks(np.arange(0,6,0.5))
        plt.savefig(rf'C:\Users\noam1\Desktop\Mitochondrial Gene comparision project\Redoing the work with Bio\Codons\Codon_usage\{phyl}\{phyl}_{aa}.pdf')
        
    #%%
    Id = full_taxa.loc['Homo_sapiens','Replicons']
    path = handle_fetcher(Id)
    record = SeqIO.read(path,'genbank')
    sequence=record.seq
    sequence = unite_seq(record, sequence)
    full_taxa=pd.read_csv(r'C:\Users\noam1\Desktop\Mitochondrial Gene comparision project\Redoing the work with Bio\Codons\organisms_with_codon_usage.csv')