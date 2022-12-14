% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CCRStoAAC.R
\name{CCRStoAAC}
\alias{CCRStoAAC}
\title{CCRStoAAC}
\usage{
CCRStoAAC(
  gene = NULL,
  chromosome = NULL,
  nproc = NULL,
  gnomad_version = "gnomAD3_0",
  ens_version = 101,
  out_path = "out/",
  ccr_file =
    "data/rawCCRs/gnomad3_0/vep101/sort_weightedresiduals-cpg-synonymous-novariant.txt.gz",
  mapping_file = "data/mapping_tables/ensembl_uniprot_MANE_metrics_07102020.tsv.gz",
  gtf = "data/GTF101/",
  fastas = "data/fastas",
  keep = 0,
  logfile = "/dev/null"
)
}
\arguments{
\item{gene}{Genes of interest, comma separated. Not specifying triggers by default the mapping of all genes.}

\item{chromosome}{Only autosomes and X chromosome. Mandatory parameter.}

\item{nproc}{Number of cores for parallel processing}

\item{gnomad_version}{gnomAD database version. Default is "gnomad3_0". It will be used for naming output files}

\item{ens_version}{Ensembl database version employed for annotating gnomAD VCF files. Will define the set of transcripts for the #' #' mapping}

\item{out_path}{Path to output folder}

\item{ccr_file}{Path to the raw CCRs file, obtained after running the CCRs pipeline}

\item{mapping_file}{Path to the file with the mapping of IDs bewteen UniProtKB, Ensembl, MANE and gnomAD constraint scores}

\item{fastas}{Path to the folder where protein fasta files are/will be organized by chromosome (will be extended with 'chr')}

\item{keep}{keep=1 keeps previous output and resume after last mapped gene, keep=0 overwrites}

\item{logfile}{Path to a log file where the log of the run will be dumped. Defaul logfile="/dev/null" the output will be discarded (only works in UNIX systems)}

\item{gtf_file}{Path to folder where Ensembl GTF files are, splitted by chromosome for an easier lookup}
}
\value{
The function itself calls to do the mapping and saves the results to a file (e.g. \emph{"./out/gnomAD3_0/vep_101/aac_weightedresiduals-cpg-synonymous-novariant_##.tsv"}), returning nothing at it's end.
This is because the mapping tables can be quite large depending on the amount/length of the genes you requested.
You can then read-in the output files by-chromosome, for that using \pkg{data.table::fread()} is highly recommended
}
\description{
This code is for mapping the Human Constrained Coding Regions (CCRs) from genomic
coordinates (chromosomes 1-22 and X) to amino acid sites in protein sequences
of UniProtKB, via reference annotated Ensembl transcripts of the GENCODE basic Project.

Mandatory input files are set by default and can be found in different sub-folders in the \emph{'./data/'} sub-directory of this package

External dependencies (non R libraries):
\pkg{blastp, awk, grep, tr, sort}
}
\examples{
 Mapping CCRs for only one gene
 \code{CCRStoAAC( gene="VPS4B" )}

 Mapping CCRs for complete chromosome, specifying 10 processors to run in parallel and starting from a previous
 output file (i.e. \emph{./out/gnomAD3_0/vep_101/aac_weightedresiduals-cpg-synonymous-novariant_18.tsv})
 \code{CCRStoAAC( chromosome=18, nproc=10, keep=1 )}

 Mapping CCRs for a list of three genes
 \code{CCRStoAAC( gene="WDR7,VAPA,VPS4B", gnomad_version="gnomAD3_0")}

}
