## Repository for "The metazoan landscape of mitochondrial DNA gene order and content is shaped by selection and affects mitochondrial transcription" by Noam Shtolz and Dan Mishmar
This code contains 4 main parts:

- 01 Database construction
- 02 Gene content analysis 
- 03 Gene order analysis
    - 03.01 Gene order-based distance clustering
    - 03.02 Gene order variance analysis (AR analysis)
    - 03.03 Gene neighbor prevalence calculations
- 04 Poly-cistrone detection using transcriptomic data (RNA-seq and PRO-seq)
    - 04.01 Next generation sequencing (NGS) data collection and preparation
    - 04.02 RNA-seq data analysis
    - 04.03 Precision run-on sequencing (PRO-seq) and global run-on sequencing (GRO-seq) analyses

### Background
This repository contains all the code used to calculate and generate the results presented in the paper `link here`.
In the paper, we study and compare the mitochondrial gene order and contents of over 8,000 available metazoan oragnisms, downloaded from NCBI Organelle database. We then analyzed available RNA-seq data of 56 different metazoans to disover the polycistronic organization of a variaty of mitochondrial DNAs (mtDNAs) based on intergenic junction expression patterns.
### Installation steps
- Clone this repository into any folder.
- Use `pip install -r requirements.txt` to install all dependencies used.
- The code should be run in the same order listed: 01 > 02 > 03 > 04. In general, `.py` files contain functions and should not be run directly, the functions within `.py` files can be run by running the `.ipynb` files.