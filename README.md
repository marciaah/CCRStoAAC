### CCRStoAAC
R package for mapping the Constrained Coding Regions (CCRs) in the Human genome to proteins


### Overview
Constrained Coding Regions (CCRs) are focal regions in the Human coding genome depleted of protein changing variations (i.e. missense, stop gain/loss, frameshift indels) [Havrilla et al., 2019, Nature Genetics](https://doi.org/10.1038%2Fs41588-018-0294-6) [GitHub](https://github.com/quinlan-lab/ccr). These regions were originally identified using the whole exome sequencing data from large cohorts of healthy control populations aggregated in gnomAD2.0 (Genome Aggregation Database 2.0, GRCh37/hg19, 125.748 human exomes). 
Here, we extended this by calculating the CCRs using gnomAD3.0 (GRCh38, 76.156  Human genomes) and mapping these regions to the amino acids in Human protein sequences of UniProtKB

- Citation: Manuscript in preparation 
- Dependencies:
  - A linux based OS
  - R packages:
  - Blastp
  - AWK 

### Input files:
In order to do the mapping, the package needs the following files. They are all provided in the data/ folder 
  - The CCRs raw file, as it comes as output from the CCRs model pipeline  [GitHub](https://github.com/quinlan-lab/ccr). We provide this file, obtained using gnomAD3.0 (GRCh38) variants and annotations of VEP101, in the folder data/rawCCRs/gnomad3_0/vep101/sort_weightedresiduals-cpg-synonymous-novariant.txt.gz
  - A mapping table that provides the identifiers correspondences between [UniProtKB](https://www.uniprot.org/) (version 10-2020) and [ENSEMBL](https://www.ensembl.org/index.html) (version 101). You can find this table in data/mapping_tables/ensembl_uniprot_MANE_metrics_07102020.tsv.gz. It also includes the identification of MANE transcripts and the constraint scores of gnomAD2.1 for genes and transcript, when available
  - Pre-filtered and simplified GTF files from ENSEMBL (v101). This is used to identify and select the [GENCODE basic transcripts](https://www.gencodegenes.org/human/) (version 35), a subset of representative transcripts for each gene. This subset prioritises full-length protein coding transcripts over partial or non-protein coding transcripts within the same gene, and intends to highlight those transcripts that will be useful to the majority of users
### Output files
  - The output with the mapping will be dumped in /out/gnomAD3_0/vep_101/aac_weightedresiduals-cpg-synonymous-novariant_##.tsv, where ##= the chromosome name
  - Columns in this file:
     - chr: chromosome 
     - pos: genomic position
     - ensembl_gene_name: Ensembl gene name
     - strand: strand (+ or -)
     - ensembl_gene_id: Ensembl gene ID (i.e. ENSG00000######)
     - weighted_pct: the CCRs percentile
     - simple_uniprot_id_SP_C: if the genomic position encodes for an amino acid in a UniProtKB/SwissProt Canonical protein (SP_C), it will be in this column (e.g. Q9P0L0-1,247,F Phe in position 247 of protein Q9P0L0-1) 
     - ensembl_transcript_id-ensembl_protein_id-uniprot_id-pos-aac_SP_C: if the genomic position encodes for an amino acid in a UniProtKB/SwissProt Canonical protein (SP_C), this will be here and also the Ensembl transcript and protein identifiers (e.g. ENST00000400000,ENSP00000382880,247,F,Q9P0L0-1,247,F)
     - ensembl_transcript_id-ensembl_protein_id-uniprot_id-pos-aac_SP: if the genomic position encodes for an amino acid in a UniProtKB/SwissProt non-canonical protein (SP), this will be here and also the Ensembl transcript and protein identifiers (e.g. ENST00000340541,ENSP00000345656,292,F,Q9P0L0-2,292,F)
     - ensembl_transcript_id-ensembl_protein_id-uniprot_id-pos-aac_TR: if the genomic position encodes for an amino acid in a UniProtKB/trEMBL protein (TR), this will be here and also the Ensembl transcript and protein identifiers (e.g. ENST00000624697,ENSP00000485643,157,S,A0A096LPJ4,157,S)
     - ensembl_transcript_id-ensembl_protein_id-uniprot-id-pos-aac_NM: if the genomic position encodes for an amino acid in a Ensembl protein that does not match in sequence (NM) with a UniProtKB protein, showing only Ensembl transcript and protein identifiers (e.g. ENST00000602528,ENSP00000501156,34,G,no_match,NA,NA)
     - chr-start-end-ensembl_gene_name: identifier, to unequivocally identify genomic regions (e.g. 1-1014475-1014475-ISG15) 
### Input parameters:
- gene: Genes of interest, comma separated. Passing nothing triggers the mapping of all genes
- chromosome: Only autosomes and X chromosome. Mandatory parameter.
- nproc: Number of cores for parallel processing
- gnom_version: gnomAD database version. Default is "gnomad3_0". It will be used for naming output files
- ens_version: Ensembl database version employed for annotating gnomAD VCF files. Will define the set of transcripts for the mapping
- out_path: Path to output folder
- ccr_file: Path to the raw CCRs file, obtained after running the CCRs pipeline
- mapping_file: Path to the file with the mapping of IDs bewteen UniProtKB, Ensembl, MANE and gnomAD constraint scores
- gtf_file: Path to folder where Ensembl GTF files are, splitted by chromosome for an easier lookup
- fastas: Path to the folder where protein fasta files are/will be organized by chromosome (will be extended with 'chr')
- keep: keep=1 keeps previous output for same chromosome and resumes after last mapped gene, keep=0 overwrites the file. 


### Use examples:
- CCRStoAAC(gene="VPS4B,VAPA,WDR7"), will map to proteins the CCRs from genes "VPS4B,VAPA,WDR7"  
- CCRStoAAC(chromosome=18,out_path="out/",nproc=15,keep=1), triggers the mapping from all genes in chromosome 18 that have CCRs percentiles calculated, but starts from a previous mapping file for the same chromosome, if it exists. Useful when a previous run had stopped (e.g. memory overload). It will use 15 cores of your computer.
- CCRStoAAC(), triggers the mapping of all genes in autosomes and X chromosome that have CCRs calculated. By default it will use 10 cores of your computer.

