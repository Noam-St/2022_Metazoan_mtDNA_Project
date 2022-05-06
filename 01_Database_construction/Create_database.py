""" 
This code is the basis of my project, creating a Metazoan database from NCBI Organelle database 
using the Bio package in python.
This code contains the main Retreive_details function which requires a RefSeq ID, it constructs
a dictionary with several key features of the organism:
scientific_name; kingdom; phylum; class; order; family; genus; genetic_codon_table; haplotype; topology; completeness; gene_order; Introns.
Additionally, the main function takes advantage of four other functions: Which_tRNA, gene_fasta, check_tRNA and Intron_handler,
Which_tRNA unifies gene annotations, gene_fasta creates fasta files with the sequence for a specific gene,
check_tRNA validates a chosen tRNA using up to two tRNA detecting algrorithms, tRNAscan-SE and/or ARWEN and 
Intron_handler is activated in case an intron exists in order to annotate it appropriately.
Dependencies: Bio; re; os; sys; numpy; pandas; urlib.
"""
#%%
import pandas as pd
import re,os
from consts import TRNA_DICT, REPL_DICT
from Bio import Entrez
from Bio import SeqIO
from ast import literal_eval
from Bio import Medline
from Bio.Seq import UndefinedSequenceError
from sys import platform
import multiprocessing as mp
import random
import urllib.request
path = os.getcwd()
if platform != 'win32':
    dash = '/'
else:
    dash = '\\'


def dir_checker(dirname : str) -> None:
  """
  check if a directory exists, if it does not, create it.

  Parameters
  ----------
  dirname : str
    Name of the directory to be checked
  """
  dirname = dirname.replace('/', dash)
  if not os.path.isdir(f'{path}{dash}{dirname}'):
    os.mkdir(f'{path}{dash}{dirname}')

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

#%% tRNA annotation unification code
def Which_tRNA(name):
    """
    This code recieves a gene name, tries to match it with common tRNA motifs, and then with any matching values in the replacement dict, if nothing matches, returns None
    which is later removed
    """
    if name in REPL_DICT.keys(): return name 
    missing_log = open(f'{path}{dash}logs{dash}missing_log.txt', 'a')
    c=0
    name=str(name)
    if 'transfer' in name:
        name = name.replace('transfer ','t')
    codon = re.search(r'\(\w{3}\)',name) # Grab codon if it exists
    if codon: codon = codon.group(0)
    else: codon=''
    name = re.sub(pattern=r'\(\w{3}\)',repl = '',string = name)
    if name in TRNA_DICT.keys():  # if the tRNA is in tRNA-X format, replace with my format and return
        name = TRNA_DICT[name]
        return name+codon

    #Following lines try to match with each possible tRNA.
    if re.match(r't[r,R]\D{1,2}(H|His)',name):
        name='tRNA-His'
        
    if re.match(r't[r,R]\D{1,2}(K|Lys)',name):
        name='tRNA-Lys'
        
    if re.match(r't[r,R]\D{1,2}(R|Arg)',name):
        name='tRNA-Arg'
        
    if re.match(r't[r,R]\D{1,2}(D|Asp)',name):
        name='tRNA-Asp'
        
    if re.match(r't[r,R]\D{1,2}(E|Glu)',name):
        name='tRNA-Glu'
        
    if re.match(r't[r,R]\D{1,2}S',name):
        name='tRNA-Ser'
        
    if re.match(r't[r,R]\D{1,2}T',name):
        name='tRNA-Thr'
        
    if re.match(r't[r,R]\D{1,2}(N|Asn)',name):
        name='tRNA-Asn'
        
    if re.match(r't[r,R]\D{1,2}(Q|Gln)',name):
        name='tRNA-Gln'
        
    if re.match(r't[r,R]\D{1,2}V',name):
        name='tRNA-Val'
        
    if re.match(r't[r,R]\D{1,2}L',name):
        name='tRNA-Leu'
        
    if re.match(r't[r,R]\D{1,2}I',name):
        name='tRNA-Ile'
        
    if re.match(r't[r,R]\D{1,2}M',name):
        name='tRNA-Met'
        
    if re.match(r't[r,R]\D{1,2}(F|Phe)',name):
        name='tRNA-Phe'
        
    if re.match(r't[r,R]\D{1,2}(Y|Tyr)',name):
        name='tRNA-Tyr'
        
    if re.match(r't[r,R]\D{1,2}(W|Trp)',name):
        name='tRNA-Trp'
        
    if re.match(r't[r,R]\D{1,2}P',name):
        name='tRNA-Pro'
        
    if re.match(r't[r,R]\D{1,2}G',name):
        name='tRNA-Gly'
        
    if re.match(r't[r,R]\D{1,2}C',name):
        name='tRNA-Cys'
    
    if re.match(r't[r,R]\D{1,2}A',name):
        name='tRNA-Ala'
        
    for k,v in REPL_DICT.items(): # if not tRNA, this loop matches according to replacement dict
        if name in v:
            name=k
            return name
        if c==len(REPL_DICT)-1:
            name=None
            return name
        c=+1
    try:
        if name in TRNA_DICT.keys(): # if matched with any of the if's above (for tRNAs), replace with my annotation accordingly.
            name = TRNA_DICT[name]
            return name+codon
    except KeyError:
        pass
    print(name)
    missing_log.write(f'%%%%%%%%%%\n{name}\n%%%%%%%%%')
    if 'RNA' in name or 'rna' in name:
        return None
    elif 'orf' in name or 'ORF' in name:
        return 'ORFX'
    else:
        return None


#%%
""" handle = Entrez.efetch(db ='nucleotide', id = 'NC_006886', rettype = 'gb', retmode = 'text') # grab gene_features for urrent organism based on ID
record = SeqIO.read(handle, 'genbank')
org = ' '.join(record.annotations['organism'].split(' ')[0:2])

pre_tax = Entrez.esearch(db = 'Taxonomy', term = org) #Grab TAX_ID based on RefSeq ID
pre_tax_record = Entrez.read(pre_tax)
TAX_ID = pre_tax_record['IdList'][0]
lineage = pre_tax_record[0]['LineageEx'] """

#%%
def gene_fasta(gene,seq_range,seq):
    """"" Recieves gene seq range and full org seq file, creates a fasta file with the gene sequence
    and returns the path to that file """
    gene = re.sub(pattern=r'\(\w{3}\)',repl = '',string = gene)
    gene = gene.replace('-','')
    gene = gene.replace('*','')
    d = random.random()
    genepath = f'{path}{dash}temp{dash}{gene}{d}.fasta'
    start,end,strand = seq_range.split(':') # Split the range file i created in the following format : (stard:end:strand) to its sub-components
    start = int(str(start).replace('>','').replace('<',''))
    end = int(str(end).replace('>','').replace('<',''))
    strand = int(strand)
    if start > 50: addition = 50
    else: addition = 0
    if (end - start) > 100:
        print(f'Abnormally long tRNA! gene: {gene}')
    try:
        if strand == 1: #if positive strand
            gene_seq = str(seq[start-addition:end+50])
        else:
            gene_seq = str(seq[start-addition:end+50].reverse_complement()) #if negative strand, take reverse complement instead
    except UndefinedSequenceError:
        return None
    print(gene_seq)
    with open(genepath, 'w') as fasta: # creates fasta file
        fasta.write(f'>{gene}\n')
        fasta.write(gene_seq)
    return genepath # returns the path to the created fasta file.


def check_tRNA(fasta, trnsl, gene):
    """" fasta with tRNA sequence, uses algorithms to validate it, returns alternative tRNA if annotation is wrong
    or True if annotation is correct MUST DELETE FASTA.

    tRNAscan-SE form: tRNAscan-SE -O {input_path} -a {output_path}
    ARWEN form: ./{compiled_code_path} -gc{translation_table_code} -l -s {output_path} {input_path}"""
    d = random.random()
    gene = re.sub(pattern=r'\(\w{3}\)',repl = '',string = gene)
    gene = gene.replace('-','')
    gene = gene.replace('*','')
    try: os.remove('/tmp/tscan14053.fpass')
    except OSError: pass
    logfile = open(f'{path}{dash}logs{dash}check_tRNA_log.txt', 'a')
    command_tRNAscan = f'tRNAscan-SE -q -O {fasta} -a {path}{dash}outputs{dash}tRNA{d}.txt' #tRNAscan-SE command structure
    os.system(command_tRNAscan) #send tRNAscan-SE command to bash
    with open(f'{path}{dash}outputs{dash}tRNA{d}.txt', 'r') as result: #open the resultant file
        output_tscan = result.read() #read file into string
        if output_tscan != '': #if string not empty, meaning tRNA was found
            output_tscan = output_tscan.split('\n')[0].split(' ')[1].split(':')[0].replace('-','')
            output_tscan = Which_tRNA(output_tscan) #grab just the tRNA name
            logfile.write(f'tRNAscan-SE: {output_tscan}\n')
            print(output_tscan)
        os.remove(f'{path}{dash}outputs{dash}tRNA{d}.txt') #remove result file created
    if gene != output_tscan or output_tscan == '': #if no result found or different result than gene
        command_ARWEN = f'.//{ARWEN_FILE} -gc{trnsl} -l -s -fo -o {path}{dash}outputs{dash}ARWEN{d}.txt {fasta}' #ARWEN command
        os.system(command_ARWEN) #send ARWEN command to bash
        with open(f'{path}{dash}outputs{dash}ARWEN{d}.txt', 'r') as result: #open ARWEN result file
            output_arwen = result.read()
            if output_arwen != '': #if string not empty, meaning tRNA was found
                TV = 0
                output_arwen = output_arwen.split('\n')[0].split(' ')
                if '?' in output_arwen[0]: #If ? in output, that means ARWEN is considering more than one tRNA, split the resultant tRNAs and find which one equals the recieved tRNA.
                    try:   
                        trna1 = Which_tRNA('tRNA-'+output_arwen[0].split('|')[0].replace('(','').replace(')',''))
                        trna2 = Which_tRNA('tRNA-'+output_arwen[0].split('|')[1].replace('(','').replace(')',''))
                        if trna1 == gene:
                            output_arwen = trna1
                        elif trna2 == gene:
                            output_arwen = trna2
                        else:
                            output_arwen = ''
                    except IndexError:
                        print('Unknown tRNA')
                        output_arwen = ''
                else:
                    if len(output_arwen) > 2: TV = 1
                    output_arwen = output_arwen[TV].replace('>','')[1:9]
                    output_arwen = Which_tRNA(output_arwen) #grab just the tRNA name
                    logfile.write('ARWEN: ' + output_arwen+'\n')
                    print(output_arwen)
            os.remove(f'{path}{dash}outputs{dash}ARWEN{d}.txt') #delete the ARWEN result file            
        if output_arwen == '' and output_tscan == '': #if both algorithms didnt find anything, return None
            print(f'Error! Could not validate {gene} in either algorithm! skip organism!')
            logfile.write(f'Error! Could not validate {gene} in either algorithm! skip organism!\n')
            return None
        if output_arwen == output_tscan: #if both algorithms cofirmed the same correct result, return the correct tRNA
            print(f'{gene} is mis-annotated! returning true annotation confirmed by both algorithms')
            logfile.write(f'{gene} is mis-annotated! returning true annotation confirmed by both algorithms\n')
            return output_arwen
        elif output_arwen == gene: #if only ARWEN is correct, return gene, annotation is assumed to be correct
            print(f'Only ARWEN confirmed gene, annotation is correct {gene}')
            logfile.write(f'Only ARWEN confirmed gene, annotation is correct {gene}\n')
            return gene
    print(f'Gene validated {gene}')
    logfile.write(f'Gene validated {gene}\n')
    logfile.close()
    return gene #normal result, tRNAscan-SE found tRNA and ARWEN wasnt needed.

def intron_handler(record,genetic_code,topology,seq):
    """ Special function to generate gene_order incase the organism has introns, returns 
    gene order and gene locations as a list, introns with genes are split into 5' and 3'"""
    gene_list=[]
    gene_list2=[]
    loc_list=[]
    note_list=[]
    tRNA_validation=[]
    checker = False
    intron_gene = 0
    gene_start = 0
    for feature in record.features: #iterate over gene features
        feat_type = feature.type
        if feat_type == 'gene': #grab all ORF genes (protein coding) locations
            gene_start = feature.location.start
            gene_end = feature.location.end
            gene_strand = feature.location.strand
        if feat_type == 'intron': #look for intron marks
            intron_start = feature.location.start 
            intron_end = feature.location.end
            if intron_start >= gene_start and intron_end <= gene_end: #if intron, and if within latest gene range, split gene into 3' end and 5' end (opposite for complementry strand gene)
                print(f'Intron is inside {gene}!')
                intron_gene_start = gene_start
                intron_gene_end = gene_end
                intron_gene_strand = gene_strand
                tRNA_validation.append('N')
                checker = True
                if gene_strand == 1:
                    intron_gene = gene
                    gene_list.append(f'5_{gene}')  
                    note_list.append(note)
                    loc_list.append(f'{gene_start}:{gene_end}:{gene_strand}')
                else:
                    intron_gene = gene
                    gene_list.append(f'-3_{gene}')
                    note_list.append(note)
                    loc_list.append(f'{gene_start}:{gene_end}:{gene_strand}')
        if feat_type == 'CDS' or feat_type == 'tRNA' or feat_type == 'rRNA': #if feature is a gene/tRNA/rRNA
            gene = str(feature.qualifiers.get('product')).replace('[','').replace(']','').replace('\\','').replace('\'','') #Grab product
            gene = Which_tRNA((gene)) #Send to annotation unification function
            codon =  str(feature.qualifiers.get('codon_recognized')).replace('[','').replace(']','').replace('\\','') #grab codon
            note = str(feature.qualifiers.get('note')) #grab note
            CDS_start = feature.location.start #grab feature start location
            CDS_end = feature.location.end #grab feature end location
            CDS_strand = feature.location.strand #grab feature strand
            loc = f'{CDS_start}:{CDS_end}:{CDS_strand}'
            if checker == True and CDS_start > intron_gene_end: # if current feature is outside of intron_gene range that means we reached the end of intron_gene
                tRNA_validation.append('N')
                loc_list.append('')
                note_list.append('None')
                checker = False
                if intron_gene_strand ==1: #if intron gene is in positive strand
                    gene_list.append(f'{intron_gene}_3')
                else:
                    gene_list.append(f'-{intron_gene}_5')
            if gene==None: #if Which_tRNA returned none, remove it from gene-list
                print('Which_tRNA returned None!')
                continue
            if gene in TRNA_DICT.values(): #if gene is a tRNA 
                fasta = gene_fasta(gene=gene, seq_range=loc, seq=seq) #first make fasta file for checking algorithms
                if fasta == None:
                    print('Invalid fasta! ignoring organism!')
                    return None
                gene_check = check_tRNA(fasta=fasta, trnsl=genetic_code, gene=gene[0:4])
                os.remove(fasta)
                if gene_check != None: #if validated
                    gene = gene_check
                    tRNA_validation.append('V')
                else:
                    tRNA_validation.append('X') #if not validated
            else : tRNA_validation.append('N')
            if codon != 'None' and not re.match(r'\(\w{3}\)',gene):
                gene+=f'({codon})'.replace('\'','')
            gene_list2.append(gene) #appends to secondary gene list without strand or duplicate annotations
            counter=gene_list2.count(gene) # counts how many copies of a certain gene exist in the un-annotated list
            gene=gene+'*'*(counter-1) #adds duplicate signs based on amount of copies
            if CDS_strand == -1:
                gene_list.append(f'-{gene}')
            else:
                gene_list.append(gene)
            note_list.append(note)
            loc_list.append(loc)
        if intron_gene in gene_list: #the code written caputes the intron_gene name aswell (without _5 or _3) so this line removes it form gene_list and derivative lists.
            intron_ind = gene_list.index(intron_gene)
            del gene_list[intron_ind]
            del note_list[intron_ind]
            del loc_list[intron_ind]
            del tRNA_validation[intron_ind]

    return gene_list,note_list,loc_list, tRNA_validation
try:
    org_df = pd.read_csv(f'{path}{dash}DB_csvs{dash}final_2021.csv',index_col = 0)
    org_list = list(org_df.index)
except FileNotFoundError:
    org_list = []



def retrieve_details(ID):
    """" Recieves organism NCBI ID and returns gene order, and codon usage
    converts into universal format,returns a dict of gene_order,codon,taxonomy(?)"""
    ranks = ['domain','kingdom','phylum','class','order','family','genus']
    filename = f'{path}{dash}genbank_DB{dash}{ID}.gbk'
    missing_log = open(f'{path}{dash}logs{dash}missing_log.txt', 'a')
    error_log = open(f'{path}{dash}logs{dash}error_log.txt', 'a')

    if 'CM' in ID:
        error_log.write(f'NOT REFSEQ ID: {ID}\n') 
        return None
    try:
        if not os.path.isfile(filename): # Check if genebank file exists on system, if not creates it
            print(f'Downloading {ID}\n')
            net_handle = Entrez.efetch(db ='nucleotide', id = ID, rettype = 'gb', retmode = 'text') # grab gene_features for current organism based on ID
            out_handle = open(filename, 'w')
            out_handle.write(net_handle.read())
            out_handle.close()
            net_handle.close()
            print(f'Saved {ID}\n')
    except:
        print(f'File does not exist! ID:{ID}')
        error_log.write(f'Genbank file not found! ID: {ID}\n')
        return None
    print(f'Parsing {ID}\n')
    return_dict = {}
    record = SeqIO.read(filename, 'genbank') # Read genbank file record into variable
    org = record.annotations['organism']
    master_log = open(f'{path}{dash}logs{dash}master_log.txt','a+')
    if not os.path.isfile(f'{path}{dash}tax_DB{dash}{ID}.txt'):
        master_log.write(f'Creating taxfile: {ID}\n')
        pre_tax = Entrez.esearch(db = 'Taxonomy', term = org) #Grab TAX_ID based on RefSeq ID
        try: pre_tax_record = Entrez.read(pre_tax)
        except RuntimeError:
            print(f'Search backend failed! {org}')
            return None
        try:
            TAX_ID = pre_tax_record['IdList'][0]
        except KeyError: #If tax_record doesn't exist, report to error_log and skip organism
            print('Tax file not found')
            error_log.write(f'Tax file not found! ID: {ID}\t org: {org}\n')
            return None
        except IndexError:
            print('No info in taxfile')
            error_log.write(f'Tax file doesnt have info! ID: {ID}\t org: {org}\n')
            return None
        master_log.write(f'Tax parsing {TAX_ID}')
        tax_handle = Entrez.efetch(db = 'Taxonomy', id = TAX_ID, retmode = 'xml')
        tax_record = Entrez.read(tax_handle) # read taxonomy file into variable.
        try:
            lineage = tax_record[0]['LineageEx'] # Grab the dictionary of lineage
        except KeyError: #If organism does not have a lineage report, report to error_log and skip organism
            print('No lineage on taxfile! (KeyError)')
            error_log.write(f'No lineage on organism! (KeyError) ID: {ID}\t org: {org}\n')
            return None
        except IndexError:
            print('No lineage on taxfile (IndexError)')
            error_log.write(f'No lineage of organism! (IndexError) ID: {ID}\t org: {org}\n')
            return None
        return_dict['genetic_code'] = tax_record[0]['MitoGeneticCode']['MGCId'] # grab genetic code table ID
        with open(f'{path}{dash}tax_DB{dash}{ID}.txt', 'w') as taxfile: #if tax file not created, creates it now.
            taxfile.write(str(lineage).replace('[','').replace(']','')+'\n')
            taxfile.write(str(return_dict['genetic_code']))
    else:
        master_log.write(f'Taxfile exists: {ID}\n')
        with open(f'{path}{dash}tax_DB{dash}{ID}.txt', 'r') as taxfile: #If taxfile exists, parses it into the return dict.
            tax = taxfile.read()
            lineage = literal_eval(tax.split('\n')[0])
            return_dict['genetic_code'] = int(tax.split('\n')[1])
    for rank in lineage:  # If lineage is one of dear king phillip came over for good _ then add it to return dict
        if rank['Rank'] in ranks:
            return_dict[rank['Rank']] = rank['ScientificName']
    org = ' '.join(org.split(' ')[0:2])
    if org in org_list:
        error_log.write(f'Organism already in .csv! ID: {ID}\t org: {org}\n')
        return None
    else:
        org_list.append(org)
    missing_log.write(f'### {org}\n')
    return_dict['organism'] = org
    return_dict['topology'] = record.annotations['topology'] #grab mtDNA topology
    if return_dict['topology'] == 'fragmented': return_dict['Full'] = False #if topology is fragmented mark partial genome
    temp = record.features
    seq = record.seq
    gene_list = []
    gene_list2 = []
    note_list = []
    loc_list = []
    tRNA_validation = []
    intron = False
    for seqfeat in temp: #Iterate over gene_features
        product = str(seqfeat.qualifiers.get('product'))
        codon_temp = str(seqfeat.qualifiers.get('codon_recognized'))
        if seqfeat.type == 'intron':
            intron = True
            break
        if product != 'None':
            gene = product.replace('[','').replace(']','').replace('\\','').replace('\'','')
            codon = codon_temp.replace('[','').replace(']','').replace('\\','')
            note = str(seqfeat.qualifiers.get('note'))
            loc = f'{seqfeat.location.start}:{seqfeat.location.end}:{seqfeat.location.strand}' #grab the range and strand of current gene
            gene=Which_tRNA(gene) # Replace with common annotation
            if gene == None: #if Which_tRNA returned none, remove it from gene-list
                error_log.write(f'Which_tRNA returned None! ID: {ID}\t gene: {product}\n')
                continue
            if gene in TRNA_DICT.values():
                fasta = gene_fasta(gene=gene, seq_range=loc, seq=seq)
                if fasta == None:
                    print(f'Invalid fasta! ignoring organism! {org}')
                    return None
                gene_check = check_tRNA(fasta=fasta, trnsl=return_dict['genetic_code'], gene=gene[0:4])
                os.remove(fasta)
                if gene_check != None:
                    gene = gene_check
                    tRNA_validation.append('V')
                else:
                    tRNA_validation.append('X')
            else:
                tRNA_validation.append('N')
            if codon != 'None' and not re.match(r'\(\w{3}\)',gene):
                gene+=f'({codon})'.replace('\'','')
            gene_list2.append(gene) #appends to secondary gene list without strand or duplicate annotations
            counter=gene_list2.count(gene) # counts how many copies of a certain gene exist in the un-annotated list
            gene=gene+'*'*(counter-1) #adds duplicate signs based on amount of copies
            if seqfeat.strand == -1: #checks if negative strand, if it is adds - sign
                gene='-'+str(gene)
            gene_list.append(gene)
            loc_list.append(loc)
            note_list.append(note)
    # add all collected variables into the return_dict
    if intron:
        intron_lists = intron_handler(record,return_dict['genetic_code'], return_dict['topology'],seq)
        gene_list = intron_lists[0]
        note_list = intron_lists[1]
        loc_list = intron_lists[2]
        tRNA_validation = intron_lists[3]
    if return_dict['topology'] == 'circular': #If dna topology is circular, anchor gene order around nad1, to prevent incorrect gene orders
        if 'nad1' in gene_list: minus = ''
        elif '-nad1' in gene_list: minus = '-'
        else: minus = False
        if minus != False:
            print('We reached wrapper')
            nad1 = gene_list.index(f'{minus}nad1')
            gene_list = gene_list[nad1:] + gene_list[:nad1]
            note_list = note_list[nad1:] + note_list[:nad1]
            tRNA_validation = tRNA_validation[nad1:] + tRNA_validation[:nad1]
            loc_list = loc_list[nad1:] + loc_list[:nad1]
            
    return_dict['Gene_locations'] = loc_list
    return_dict['Gene_order'] = gene_list
    return_dict['Notes'] = note_list
    return_dict['tRNA_validation'] = tRNA_validation
    return_dict['Haplotype']=''
    return_dict['Intron'] = intron
    for line in record.format('genbank').split('\n')[0:40]: # This loop searches for the PUBID
        if 'COMPLETENESS' in line: # look for completeness report
            if 'full' in line:
                print('Full length indeed!')
                return_dict['Full'] = True
            else:
                return_dict['Full'] = False
        if 'PUBMED' in line:
            pub_ID = line.split(' ')[-1]
    try:
        if return_dict['class'] == 'Bivalvia': 
            try:
                medline_handle = Entrez.efetch(db = 'pubmed', id = pub_ID, rettype = 'medline', retmode = 'text') # grab the pubmed record
                medline_record = Medline.parse(medline_handle)
                for medline in medline_record: #within the title and abstract, searches whether the haplotype is male or female.
                    TI = medline['TI']
                    AB = medline['AB'] #ISSUE HERE, female contains all letters of male.
                    if re.search(r'female',TI,flags = re.IGNORECASE) or re.search(r'female',AB,flags = re.IGNORECASE):
                        return_dict['Haplotype'] = 'female'
                    elif re.search(r'male',TI,flags = re.IGNORECASE) or re.search(r'male',AB,flags = re.IGNORECASE):
                        return_dict['Haplotype'] = 'male'
                    else:
                        return_dict['Haplotype'] = 'unknown'
            except NameError:
                print(f'No PUBMED ID in this organism: {org}')
                error_log.write(f'No PUBMED ID on this organism features: {org}')
                return_dict['Haplotype'] = 'unknown'
            except KeyError:
                print(f'PubMed paper missing!')
                return_dict['Haplotype'] = 'unknown'
        else: return_dict['Haplotype'] = ''
    except KeyError:
        return_dict['Haplotype'] = ''
    return_dict['RefSeq'] = ID
    haplo = record.features[0].qualifiers.get('haplotype')
    sex = record.features[0].qualifiers.get('sex')
    group = record.features[0].qualifiers.get('haplogroup')
    if sex != '' or haplo != '':
        if sex == 'F' or sex == 'female' or sex == 'Female' or haplo == 'F' or haplo == 'female' or haplo == 'Female' or group == 'F' or group == 'female' or group == 'Female':
            return_dict['Haplotype'] = 'female'
        elif sex == 'M' or sex == 'male' or sex == 'Male' or haplo == 'M' or haplo == 'male' or haplo == 'Male' or group == 'M' or group == 'male' or group == 'Male':
            return_dict['Haplotype'] = 'male'  
    master_log.write(str(return_dict.items()) + '\n')
    return return_dict

def make_org_df(file):
    """" Recieves file path, imports into a df, grabs gene orders from ncbi returns
    a df in this format: org|Domain|Kingdom|Subgroup|Replicons|Genes"""
    full_org = pd.read_csv(file, index_col=0)
    full_org.drop(columns = ['Organism Groups','Strain','BioSample','Size(Mb)','GC%','Type','BioProject','Release Date'],inplace= True)
    full_org.Replicons = full_org.Replicons.apply(lambda x: x.split('/')[0].replace('MT:',''))
    full_org.replace(to_replace={'MT:':'','\d{1}:':'','mt:':'','mtDNA:':'','M:':'','I:':'','II:':''},regex=True,inplace=True)

    return full_org

#%%
# tRNAscan-SE format : "tRNAscan-SE -O (fasta file path) -o (output path)"
""" handle = Entrez.efetch(db ='nucleotide', id = 'NC_012920.1', rettype = 'gb', retmode = 'text') # grab gene_features for urrent organism based on ID
record = SeqIO.read(handle, 'genbank')
for feature in record.features:
    name = str(feature.qualifiers.get('product'))
    if name == None: continue
    if 'tRNA' in name or 'trn' in name:
        loc = feature.location
        print(loc)
        print(name)
        gene_fasta(name.replace('[','').replace(']','').replace('\'',''), f'{loc.start}:{loc.end}:{loc.strand}', record.seq)

 """
if __name__ == '__main__':
#%% This block runs the retrieve_details function on all organisms in organelle.csv

    ARWEN_FILE = f'output'
    # CHECK ARWEN FILE

    try:
        open(f'{path}{dash}output')
    except FileNotFoundError:
        try:
            open(f'{path}{dash}arwen1.2.3.c')
        except FileNotFoundError:
            print('ARWEN .c file not found, Downloading...')
            URL = 'http://130.235.244.92/ARWEN/arwen1.2.3.c'
            urllib.request.urlretrieve(URL, f'{path}{dash}arwen1.2.3.c')
        print(f'Compiling ARWEN...')
        os.system(f'gcc {path}{dash}arwen1.2.3.c -o output')
    full_org = make_org_df(f'{path}{dash}organelles.csv')
    c= len(org_list) -1
    for i in ['tax_DB', 'genbank_DB', 'logs', 'DB_csvs','temp','outputs']:
        dir_checker(i)

    for ID in full_org.iloc[:, 0]:
        dict_list = []
        print(f'current position: {c}')
        return_dict = retrieve_details(ID)
        if return_dict != None:
            dict_list.append(return_dict)
            dummy = pd.DataFrame(dict_list)
            dummy.set_index(keys = 'organism',inplace = True)
            dummy.to_csv(f'{path}{dash}DB_csvs{dash}final_2021.csv',header = False,mode = 'a')
        c+=1

