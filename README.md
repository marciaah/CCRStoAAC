# CCRStoAAC

R package for mapping the Constrained Coding Regions (CCRs) in the Human genome to proteins


## Overview

Constrained Coding Regions (CCRs) are focal regions in the Human coding genome depleted of protein changing variants (i.e. missense, stop gain/loss, frameshift indels) [Havrilla et al., 2019, Nature Genetics](https://doi.org/10.1038%2Fs41588-018-0294-6) ( [GitHub](https://github.com/quinlan-lab/ccr) ). These regions were originally identified using the whole exome sequencing data from large cohorts of healthy control populations aggregated in [gnomAD](https://gnomad.broadinstitute.org/) (The Genome Aggregation Database, version 2.0, GRCh37/hg19 reference genome, including 125.748 human exomes). 
Here, we extended this by calculating the CCRs using gnomAD3.0 (GRCh38, 76.156  Human genomes) and mapping these regions to the amino acids in Human protein sequences of UniProtKB, via reference transcripts of Ensembl which are part of the [GENCODE basic transcripts](https://www.gencodegenes.org/human/) (version 35).


 

## Citation

If you find this package useful for your work, please mention us: 

**Manuscript in Press:**

*'Mapping the Constrained Coding Regions in the human genome to their corresponding proteins'*,
Marcia A. Hasenahuer, Alba Sanchis-Juan, Roman A. Laskowski, James A. Baker, James D. Stephenson, Christine A. Orengo,F Lucy Raymond, Janet M. Thornton.
Available online 21 November 2022, 167892
https://doi.org/10.1016/j.jmb.2022.167892

## Installation

### Dependencies

Before installing this package, ensure you have the following: 
  - A Unix based OS: by now this package will only run in this systems
  - [R](https://cran.r-project.org/) version >= 4.0.X. Full list of R packages is in the DESCRIPTION file ('Imports') 
  - [Blastp](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download). On Ubuntu/Debian systems you can get it as part of the ncbi-blast+ suite:
  ``sudo apt-get install ncbi-blast+``
  
  No need to create any local database of sequences! as the specific Human protein sequences will be downloaded on-demand by this package.

### How to install

From R, run the following commands:
	
``install.packages("devtools")``

``devtools::install_github("marciaah/CCRStoAAC")``


## How to use

#### Load the package
 ``library(CCRStoAAC)``

#### Some examples on how to use this package:

- Mapping CCRs for only one gene
  
  ``CCRStoAAC( gene="VPS4B" )``
 
- Mapping CCRs for a list of three genes
 
 ``CCRStoAAC( gene="WDR7,VAPA,VPS4B", gnomad_version="gnomAD3_0")``

- Mapping CCRs for a complete chromosome, specifying 10 processors to run in parallel and starting from a previous
 output file (i.e. ./out/gnomAD3_0/vep_101/aac_weightedresiduals-cpg-synonymous-novariant_18.tsv). Useful when a previous run had stopped (e.g. memory overload)
 
 ``CCRStoAAC( chromosome=18, nproc=10, keep=1 )``

- Mapping of all genes in autosomes and X chromosome that have CCRs calculated
 
 ``CCRStoAAC()`` 


## Full list of parameters

- **gene**: Gene(s) of interest, comma separated. Not specifying triggers the mapping of all genes
- **chromosome**: Only autosomes (1 to 22) and X chromosome are accepted. If not provided, will be obtained from the CCRs output file 
- **nproc**: Number of cores for parallel processing. By default, it will be max numbr of cores of your computer minus two
- **gnomad_version**: gnomAD database version. Default is "gnomad3_0", as our analysis was based on it. It will be used for naming output files
- **ens_version**: Ensembl database version employed for annotating gnomAD VCF files. Will define the set of transcripts for the mapping. Default is "101", as our analysis was based on it
- **out_path**: Path to output folder
- **ccr_file**: Path to the raw CCRs file, obtained after running the CCRs pipeline
- **mapping_file**: Path to the file with the mapping of IDs bewteen UniProtKB, Ensembl, MANE and gnomAD constraint scores
- **gtf_file**: Path to folder where Ensembl GTF files are, splitted by chromosome for an easier lookup
- **fastas**: Path to the folder where protein fasta files are/will be organized by chromosome (will be extended with 'chr')
- **keep**: keep=1 keeps previous output for same chromosome and resumes after last mapped gene, keep=0 overwrites the file. 



## Input files

The following files are essential, are all provided in the ``data/`` folder and automatically loaded by the package with default parameters

  - The CCRs raw file, this is the output from running the CCRs model pipeline  [GitHub](https://github.com/quinlan-lab/ccr). We provide this file, obtained using gnomAD3.0 (GRCh38) variants and annotations of VEP101. Here you will find in one file the constraint for autosomes and X chromosome, although they were obtained separatelly ( please, refer to [Havrilla et al., 2019, Nature Genetics](https://doi.org/10.1038%2Fs41588-018-0294-6) and [GitHub](https://github.com/quinlan-lab/ccr) for further details) 
  ``data/rawCCRs/gnomad3_0/vep101/sort_weightedresiduals-cpg-synonymous-novariant.txt.gz``

  - A file with the mapping table that provides the correspondence of identifiers between [UniProtKB](https://www.uniprot.org/) (version 10-2020) and [ENSEMBL](https://www.ensembl.org/index.html) (version 101):
  ``data/mapping_tables/ensembl_uniprot_MANE_metrics_07102020.tsv.gz``
  Here you can also find the identification of [MANE transcripts](https://www.ncbi.nlm.nih.gov/refseq/MANE/) and the [gene constraint scores](https://gnomad.broadinstitute.org/help/constraint) of gnomAD2.1 for genes and transcript, when available.

  - Pre-filtered and simplified GTF files from ENSEMBL (v101):
  ``data/GTF101/Homo_sapiens.GRCh38.101.chr##.gtf.gz``

  This is used to identify and select the protein coding [GENCODE basic transcripts](https://www.gencodegenes.org/human/) (version 35), a subset of representative transcripts for each human gene.

## Output files

The resulting dataframe with the mapping will not be returned by the function, but instead it will be dumped in one file by chromosome in ``./out/gnomAD3_0/vep_101/aac_weightedresiduals-cpg-synonymous-novariant_##.tsv``, where ##: chromosome name. This is to prevent R to crash, because it can take a significant amount of RAM memory to load the mapping all at once, depending on the number and length of genes that were requested to be mapped.

You can then read-in the output files by-chromosome, for that using ``data.table::fread()`` is highly recommended, but watch memory usage!.

Columns in the output file(s) correspond to:
  
- **chr** : chromosome 
- **pos** : genomic position
- **ensembl_gene_name** : Ensembl gene name
- **strand** : strand (+ or -)
- **ensembl_gene_id** : Ensembl gene ID (e.g. ENSG00000######)
- **weighted_pct** : the CCRs percentile
- **simple_uniprot_id_SP_C** : if the genomic position encodes for an amino acid in a UniProtKB/SwissProt Canonical protein (SP_C), it will be in this column (e.g. Q9P0L0-1,247,F Phe in position 247 of protein Q9P0L0-1) 
- **ensembl_transcript_id__ensembl_protein_id__uniprot_id__pos__aac_SPC** : if the genomic position encodes for an amino acid in a UniProtKB/SwissProt Canonical protein (SP_C), this will be here and also the Ensembl transcript and protein identifiers (e.g. ENST00000400000,ENSP00000382880,247,F,Q9P0L0-1,247,F)
- **ensembl_transcript_id__ensembl_protein_id__uniprot_id__pos__aac_SP** : if the genomic position encodes for an amino acid in a UniProtKB/SwissProt non-canonical protein (SP), this will be here and also the Ensembl transcript and protein identifiers (e.g. ENST00000340541,ENSP00000345656,292,F,Q9P0L0-2,292,F)
- **ensembl_transcript_id__ensembl_protein_id__uniprot_id__pos__aac_TR** : if the genomic position encodes for an amino acid in a UniProtKB/trEMBL protein (TR), this will be here and also the Ensembl transcript and protein identifiers (e.g. ENST00000624697,ENSP00000485643,157,S,A0A096LPJ4,157,S)
- **ensembl_transcript_id__ensembl_protein_id__uniprot_id__pos__aac_NM** : if the genomic position encodes for an amino acid in a Ensembl protein that does not match in sequence (NM) with a UniProtKB protein, showing only Ensembl transcript and protein identifiers (e.g. ENST00000602528,ENSP00000501156,34,G,no_match,NA,NA)
- **chr__start__end__ensembl_gene_name** : region ID, to unequivocally identify genomic regions (e.g. 1-1014475-1014475-ISG15) 

When it is not possible to download the protein sequences, either or both from Ensembl or UniProtKB, or if for a fragment of a CCR region there are no Ensembl/UniProtKB encoded, then identifiers, amino acids and positions will have "NA".  



