import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import Entrez
from Bio import Seq
from ast import literal_eval
import re
import glob
from utils import Which_tRNA
from consts import TRNA_DICT


PATH = os.getcwd()

def generate_trna_fasta(ID, file_name, folder):
    """
    Appends the trna sequences of a given organism to a fasta file

    Parameters
    ----------
    ID : str
        The NCBI ID of the organism
    file_name : str
        The name of the fasta file to be created
    
    Returns
    -------
    None
    """
    with open(os.path.join(PATH, folder, file_name), 'a') as f:
        record = SeqIO.read(os.path.join(PATH, 'genbank_DB', f'{ID}.gbk'), 'genbank')
        try: seq = record.seq
        except Seq.UndefinedSequenceError:
            print(f'Sequence content undefined for {ID}! ...Skipping\n')
            return
        features = record.features # list of SeqFeature objects
        for feature in features:
            if feature.type == 'tRNA':
                try:
                    feat_seq = feature.extract(seq)
                    f.write(f'>{ID}_{feature.qualifiers["product"][0]}_{feature.location.start}_{feature.location.end}_{feature.strand}\n')
                    f.write(str(feat_seq) + '\n')
                except Seq.UndefinedSequenceError:
                    print(f'Sequence content undefined for {ID}! ...Skipping\n')
                    return
                
def extract_anticodon(ID):
    """
    """
    with open(os.path.join(PATH, 'genbank_DB', f'{ID}.gbk'), 'r') as f:
        record = SeqIO.read(f, 'genbank')
        anticodons = {}
        features = record.features # list of SeqFeature objects
        for feature in features:
            if feature.type == 'tRNA':
                try: gene = (f'{TRNA_DICT[(feature.qualifiers["product"][0].replace(" ", ""))]}')
                except KeyError:
                    pass
                try:
                    # Isolate only the actual anticodon sequence
                    anticodon = feature.qualifiers['anticodon'][0][-4:-1]
                    anticodons[gene] = (re.sub('[^ATCGatcg]', '', str(anticodon)))
                except KeyError:
                    anticodons[gene] = 'None'
                
        return anticodons

def generate_trna_fasta_all(IDs, file_name, folder, n_seq_split = 100):
    """
    Appends the trna sequences of all organisms to a fasta file

    Parameters
    ----------
    file_name : str
        The name of the fasta file to be created
    
    Returns
    -------
    None
    """
    last_split = 0
    for i, ID in enumerate(IDs):
        if i % n_seq_split == 0:
            last_split = i
        generate_trna_fasta(ID, f'{file_name}.{last_split}.fasta', folder)
    
def make_run_files(fasta_folder, runfiles_folder, trnascan_outfolder):
    """
    Creates a run file for each fasta file in the fasta folder

    Parameters
    ----------
    fasta_folder : str  
        The folder containing the fasta files
    runfiles_folder : str
        The folder to store the run files
    trnascan_outfolder : str    
        The folder to store the trnascan output files
    
    Returns
    -------
    None 
    """
    fastas = glob.glob(os.path.join(PATH, fasta_folder, '*.fasta'))
    for fasta in fastas:
        fasta_ind = fasta.split("/")[-1].split(".")[1]
        with open(os.path.join(PATH, runfiles_folder, f'run.{fasta_ind}.sh'), 'w') as f:
            
            COMMAND = f'tRNAscan-SE -M vert -o {os.path.join(PATH, trnascan_outfolder, f"{fasta_ind}.txt")} {fasta}\n'
            f.write('#!/bin/env conda run -n noam_env\n')
            f.write(f'#$ -cwd\n')
            f.write(f'#$ -V\n')
            f.write(f'#$ -N trnascan_{fasta_ind}\n')
            f.write(f'#$ -q bioinfo.q\n')
            f.write(f'#$ -o {os.path.join(PATH, "stdout")}\n')
            f.write(f'#$ -e {os.path.join(PATH, "stderr")}\n')
            f.write(f'#$ -S /bin/bash\n')
            f.write(COMMAND)


if __name__ == '__main__':
    org_df = pd.read_csv(os.path.join(PATH, 'DB_csvs', 'final.csv'))
    org_df.Gene_order = org_df.Gene_order.apply(literal_eval)
    org_df.Gene_locations = org_df.Gene_locations.apply(literal_eval)
    # Convert org_df into a dataframe where each row is a tRNA gene
    # and each column is a codon
    trna_df = pd.DataFrame(columns=['trna', 'organism', 'RefSeq', 'family', 'order', 'class', 'phylum', 'AA', 'codon', 'anticodon'])
    trna_dict = {i:[] for i in ['trna', 'id', 'organism', 'RefSeq', 'family', 'order', 'class', 'phylum', 'AA', 'codon', 'anticodon']}
    for _, row in org_df.iterrows():
        anticodons = extract_anticodon(row.RefSeq)
        for i, gene in enumerate(row.Gene_order):
            gene = re.sub(pattern=r'\(\w{3}\)',repl = '',string = gene)
            gene = gene.replace('-','')
            gene = gene.replace('*','')  
            if 'trn' in gene:
                loc = row.Gene_locations[i].split(':')
                if loc[0] == 'V': continue
                trna_dict['trna'].append(gene)
                trna_dict['organism'].append(row.organism)
                trna_dict['RefSeq'].append(row.RefSeq)
                trna_dict['family'].append(row.family)
                trna_dict['order'].append(row.order)
                trna_dict['class'].append(row['class'])
                trna_dict['phylum'].append(row.phylum)
                trna_dict['AA'].append(gene[-1])
                trna_dict['id'].append(f'{loc[0]}_{loc[1]}_{loc[2]}')
                try:
                    trna_dict['anticodon'].append(anticodons[gene])
                except (KeyError, TypeError):
                    trna_dict['anticodon'].append('None')
                if trna_dict['anticodon'][-1] != 'None':
                    trna_dict['codon'].append(str(Seq.Seq(trna_dict['anticodon'][-1]).reverse_complement()))
                else:
                    trna_dict['codon'].append('None')
    trna_df = pd.DataFrame(trna_dict)
    trna_df.to_csv(os.path.join(PATH, 'trna_df.csv'), index=False)       
    generate_trna_fasta_all(trna_df[trna_df.anticodon == 'None'].RefSeq.unique(), 'trna_seqs', folder = 'trna_fastas')
    make_run_files('trna_fastas', 'runfiles', 'trnascan_out')
