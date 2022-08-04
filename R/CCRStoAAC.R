#Installing manually as roxygen2 fails to install them
suppressMessages({
  if (!requireNamespace("AnnotationHub", quietly = TRUE))
    BiocManager::install("AnnotationHub",force = T,quiet=T)
  if (!requireNamespace("Gviz", quietly = TRUE))
    BiocManager::install("Gviz",force = T,quiet=T)
})


#CCRStoAAC----
#' @title CCRStoAAC
#'
#' @description This code is for mapping the Human Constrained Coding Regions (CCRs) from genomic
#' coordinates (chromosomes 1-22 and X) to amino acid sites in protein sequences
#' of UniProtKB, via reference annotated Ensembl transcripts of the GENCODE basic Project.
#'
#' The parameters are set to default in order to demonstrate the mapping of CCRs for 3 genes of chromosome 18,
#' using as input the files organised in different sub-folders in the 'data/' directory of this package
#'
#' External dependencies (non R libraries):
#'    blastp
#'    awk
#'    grep
#'    tr
#'    sort
#'
#'
#' @param  gene Genes of interest, comma separated. Passing "all" (default) triggers the mapping of all genes in a chromosome
#' @param  chromosome Only autosomes and X chromosome. Mandatory parameter.
#' @param  nproc Number of cores for parallel processing
#' @param  gnom_version gnomAD database version. Default is "gnomad3_0". It will be used for naming output files
#' @param  ens_version Ensembl database version employed for annotating gnomAD VCF files. Will define the set of transcripts for the #' #' mapping
#' @param  out_path Path to output folder
#' @param  ccr_file Path to the raw CCRs file, obtained after running the CCRs pipeline
#' @param  mapping_file Path to the file with the mapping of IDs bewteen UniProtKB, Ensembl, MANE and gnomAD constraint scores
#' @param  gtf_file Path to folder where Ensembl GTF files are, splitted by chromosome for an easier lookup
#' @param  fastas Path to the folder where protein fasta files are/will be organized by chromosome (will be extended with 'chr')
#' @param  keep keep=1 keeps previous output and resume after last mapped gene, keep=0 overwrites
#'
#' @return The function itself calls to do the mapping and saves the results to a file, returning nothing at it's end,
#' this is because the mapping tables can be quite large for a whole chromosome for example
#'
#' @examples
#'
#'  CCRStoAAC(gene="VPS4B",chromosome=18)
#'  CCRStoAAC( gene="all",chromosome=18,out_path="out/",nproc=10)
#'  CCRStoAAC( gene="WDR7,VAPA,VPS4B", chromosome=18,nproc=NULL, gnom_version="gnomAD3_0",
#'  ens_version=101,out_path="out/",ccr_file="data/rawCCRs/gnomad3_0/vep101/sort_weightedresiduals-cpg-synonymous-novariant.txt.gz",
#'  mapping_file="data/mapping_tables/ensembl_uniprot_MANE_metrics_07102020.tsv.gz",gtf="data/GTF101/",fastas="data/fastas",keep=0)
#'
#' @export
CCRStoAAC<-function(gene=NULL, chromosome=NULL, nproc=NULL, gnom_version="gnomAD3_0", ens_version=101,out_path="out/",
                    ccr_file="data/rawCCRs/gnomad3_0/vep101/sort_weightedresiduals-cpg-synonymous-novariant.txt.gz",
                    mapping_file="data/mapping_tables/ensembl_uniprot_MANE_metrics_07102020.tsv.gz",
                    gtf="data/GTF101/",
                    fastas="data/fastas",
                    keep=0)
{

  #Avoid scientific notation, otherwise genomic positions with a big number might be printed out in scientific notation
  options(scipen=999)
  #Suppressing warnings
  options(warn=-1)
  #Suppressing dplyr summaries
  options(dplyr.summarise.inform = FALSE)

  #DEBUG
  # gene="WDR7,VAPA,VPS4B"
  # gene<-"VAPA"
  # chromosome=18
  # nproc=10
  # gnom_version="gnomAD3_0"
  # ens_version=101
  # out_path="out/"
  # ccr_file="data/rawCCRs/gnomad3_0/vep101/sort_weightedresiduals-cpg-synonymous-novariant.txt.gz"
  # mapping_file="data/mapping_tables/ensembl_uniprot_MANE_metrics_07102020.tsv.gz"
  # gtf="data/GTF101/"
  # fastas="data/fastas"
  # keep=0


  ##Checking the variables that were provided
  #Number of cores for processing
  np<- ifelse(is.null(nproc), parallel::detectCores()-2,
              ifelse( nproc > parallel::detectCores()-2,
                      stop(paste("STOP. Your computer has in total cores=" ,parallel::detectCores(),
                                 ".\nThe amount of cores you provide might causes a hung of your system if nproc > (total cores - 2)\nPlease, specify as maximum nproc=",parallel::detectCores()-2,sep="")),nproc))

  k <- keep

  gene_chr<-NULL
  ext<-NULL

  #If gene and chromosome were not provided, all of them will be pulled from the CCRs file
  if(is.null(gene) & is.null(chromosome)){
    cat("|-----------WARNING: 'gene' and 'chromosome' were not provided.\nWill map to proteins all the genes in autosomes and X chromosome which have CCRs assigned.\nThis will consume more time and resources, watch memory usage as this might makes the code to crash!\n\n")
    ccr<-getRawCCRs(ccr_file)
    #ext<-"chr"
   }else{
    #If gene was provided, but chromosome was not, it will be pulled from the CCRs file
    if(!is.null(gene) & is.null(chromosome)){

      mapgene<-unlist(strsplit(gene, ","))
      cat("|-----------'gene' were provided.\nWill try to obtain chromosomes from the CCRs file\n\n")
      ccr<-getRawCCRs(ccr_file,gene = mapgene)
      #ext<-"gene"
      }else{
        #If gene was not provided, but chromosome was, it will be pulled from the CCRs file
        if(is.null(gene) & !is.null(chromosome)){
        chr<-unlist(strsplit(as.character(chromosome), ","))
        cat("|-----------'chromosome' were provided.\nWill try to obtain 'gene' from the CCRs file\n\n")
        ccr<-getRawCCRs(ccr_file,chr = chr)
        #ext<-"chr"
        }else{
          #If both gene and chromosome were provided, I'll use only the gene and get the chromosomes from the CCRs file
          if(!is.null(gene) & !is.null(chromosome)){

            mapgene<-unlist(strsplit(gene, ","))
            cat("|-----------'gene' and 'chromosome' were provided.\nWill only use 'gene' and obtain chr from the CCRs file\n\n")
            ccr<-getRawCCRs(ccr_file,gene = mapgene)
            #ext<-"chr"
          }
        }
      }
   }

  #Obtaining the pairs gene-chromosome, this will be use from not to get the list of genes-chr to process
  gene_chr<-unique(ccr[,c("gene","chrom")])


  # if(mapgene=="all"){
  #
  #   #all genes of the chr in the raw ccrs file will be mapped
  #   g <- unique(subset(ccr, chrom==chr,select = gene))
  #
  # }else{
  #   #only mapping a list of genes (1 or n)
  #   g <- unique(subset(ccr, chrom==chr & gene %in% mapgene,select = gene))
  # }



  ### Going by chromose and by gene and calling the mapping

  if(!is.null(gene_chr)){

    for (chr in unique(gene_chr$chr)){

    #protein fasta sequences will be downloaded and saved organized by chromosome,
    #check if folders exist and create them
    fastafolder<-paste(fastas,"/",chr,"/",sep="")
    tryCatch(if(!dir.exists(fastas)) { dir.create(fastas) },
             error=function(e) {
               message(paste("\nERROR: Failed creating the fastas folder:\n",fastas,"\ndir.create() returned:"))
               message(e)
               message("\n\nExecution halted")
             }
    )
    tryCatch(if(!dir.exists(fastafolder)) { dir.create(fastafolder) },
             error=function(e) {
               message(paste("\nERROR: Failed creating the fastas folder:\n",fastafolder,"\ndir.create() returned:"))
               message(e)
               message("\n\nExecution halted")
             }
    )


    #Path to GTF file
    gtf_file<-paste(gtf,"/Homo_sapiens.GRCh38.",ens_version,".chr",chr,".gtf.gz",sep="")

    #Setting the path where the output mapping will be saved
    outfile <- paste(out_path,gnom_version,"/vep_",ens_version,"/aac_weightedresiduals-cpg-synonymous-novariant_",chr,".tsv",sep="")

    ### LOADING THE EnsDb and filtering to the chromosome of interest
    # supply version of interest (ens_version=101 for replicating our analysis)
    edb<-getEnsDB(ens_version,chr)


    ### GET THE MAPPING OF IDS for this chromosome
    gencode_basic_uniprot<-getIDs( mapping_file,gtf_file,chr,ens_version)

    ## write the header for the output file
    # only if the request is to map only one gene, the name of such gene will be in
    # the file name.
    #if(length(mapgene)==1 & mapgene!="all"){
    # if(ext=="chr"){
    #   outfile <- paste(out_path,gnom_version,"/vep_",ens_version,"/aac_weightedresiduals-cpg-synonymous-novariant_",chr,".tsv",sep="")
    #   }else{
    #       if(ext=="gene"){
    #         outfile <- paste(out_path,gnom_version,"/vep_",ens_version,"/aac_weightedresiduals-cpg-synonymous-novariant_",chr,"_",g,".tsv",sep="")
    #        }
    #

    outsubfolder<-paste(out_path,gnom_version,"/vep_",ens_version,sep="")
    outfolder<-paste(out_path,gnom_version,sep="")
    #First create the out/gnomad_###/vep_###/ folder structure, if it doesn't exist already
    tryCatch(if(!dir.exists(outfolder)) { dir.create(outfolder) },
             error=function(e) {
               message(paste("\nERROR: Failed creating the output folder:\n",outfolder,"\ndir.create() returned:"))
               message(e)
               message("\n\nExecution halted")
             }
    )
    tryCatch(if(!dir.exists(outsubfolder)) { dir.create(outsubfolder) },
             error=function(e) {
               message(paste("\nERROR: Failed creating the output folder:\n",outsubfolder,"\ndir.create() returned:"))
               message(e)
               message("\n\nExecution halted")
             }
    )



    # k=0  parameter overwrites previous outputs, k=1 identifies genes that were already mapped
    # from a previous output and attempts to start from there
    if(k==0){
      header <- paste("chr","pos","ensembl_gene_name","strand","ensembl_gene_id","weighted_pct","simple_uniprot_id_SP_C","ensembl_transcript_id-ensembl_protein_id-uniprot_id-pos-aac_SP_C","ensembl_transcript_id-ensembl_protein_id-uniprot_id-pos-aac_SP","ensembl_transcript_id-ensembl_protein_id-uniprot_id-pos-aac_TR","ensembl_transcript_id-ensembl_protein_id-uniprot-id-pos-aac_NM","chr-start-end-ensembl_gene_name",sep="\t")
      write(header, outfile,append = F)} else if(k==1) {
        message(paste("Resuming from previous ",outfile,sep=" "))
        previous<-data.table::fread(outfile)
        #to resume, these are the genes to avoid because they are already  in the previous output
        avoid <- unique(subset(previous, chr==chr,select = ensembl_gene_name))
        newgenes<-setdiff(unique(gene_chr$gene),avoid$ensembl_gene_name)
        message(paste("New genes to process: ",ifelse(length(unique(newgenes$gene))>0,unique(newgenes$gene),"none"),sep=""))
        g<-NULL

        gene_chr<-gene_chr[gene_chr$gene %in% newgenes,]
        #g$gene<-newgenes
      }



    ### CALL TO GO GENE-BY-GENE PERFORMING THE MAPPING OF GENOMIC POSITIONS ####

    #if any gene to process
    if(length(gene_chr$gene)>0){
      tryCatch( byGene(gene_chr[gene_chr$chrom==chr,],gencode_basic_uniprot,edb,ccr,outfile,chr,fastafolder,ens_version,np),
                error = function(e){
                  message(paste("\nERROR: Failed doing the mapping for your genes. \nbyGene() returned:\n"))
                  message(e)
                  message("\n\nExecution halted")
                  stop()
                }
               )
    cat(paste("Wrote mapping to ",outfile,"\n"))
    }else{
      message(paste("'",mapgene,"' IS/ARE NOT in the CCRs file. \nPlease use 'all'to specify all genes of a chromosome or a comma separated list of HGNC gene names",sep=""))

      }

    #Free up some memory space
    suppressMessages(gc())

    cat("\nDone!\n")

  }
  }
}

#getEnsDB----
#' @title getEnsDB
#' @description
#' Gets annotations of a specific version of Ensembl DB by using 'ensembldb'
#' Bioconductor R Package
#' Provide version of ensembl of your preference, but please use the same one used when
#' annotating gnomAD VCF file with VEP when running the CCRs pipeline
#' (ideally based on GRCh38/hg38) otherwise you might find discrepancies on the mapping!
#'
#' @param ens_version Version of Ensembl
#' @param chr Chromosome
#' @return A data frame
#' @export
getEnsDB<-function(ens_version,chr){


  #Avoid scientific notation, otherwise genomic positions with a big number might be printed out in scientific notation
  options(scipen=999)
  #Suppressing warnings
  options(warn=-1)

  cat(paste("|-----------Loading Ensembl v",ens_version," for Human\n",sep=""))

  # start an instance of AnnotationHub
  ah <- AnnotationHub::AnnotationHub()

  # This is how you query for all available EnsDb databases, to see what's available
  #query(ah, "EnsDb")

  # Query for Human available EnsDb database in a specific version
  ahDb <- AnnotationHub::query(ah, pattern = c("Homo sapiens", "EnsDb", ens_version))
  edb <- ahDb[[1]]

  # filter the database to the actual chromosome
  edb <- ensembldb::filter(edb, filter = ~ seq_name == chr)

  return(edb)
}


#getRawCCRs----
#' @title getRawCCRs
#' @description Reads in the CCRs file
#' Make sure the CCRs file is sorted by chrom and genomic position, this will further
#' simplify loading EnsDb.Hsapiens.v## in memory, doing it only once per chr
#'
#' @param ccr_file The location of the raw CCRs file, obtained from running the model
#' @param chr Chromosome
#' @param gene Gene
#' @return A data frame with CCRs with the regions extended to each genomic base inside its boundaries
#' @export
getRawCCRs<-function(ccr_file,chr=NULL,gene=NULL){

  #Avoid scientific notation, otherwise genomic positions with a big number might be printed out in scientific notation
  options(scipen=999)
  #Suppressing warnings
  options(warn=-1)

  cat("|-----------loading CCRs file\n")
  ccr <- as.data.frame(data.table::fread(ccr_file))
  names(ccr)[1]<-"chrom"

  #simplifies the comparisons in if()
  chr<-ifelse(is.null(chr),NA,chr)

  chromosomes<-c(seq(1,22,by=1),"X")


  #If gene/s and chromosome/s were provided, filter according to this
  if(!is.na(chr) & !is.null(gene) ) {
    cat("|-----------WARNING: both 'gene' and 'chromosome' were provided.\nFiltering CCRs file by 'gene' and obtaining 'chromosome' from it\n\n  ")
    ccr<-ccr[ccr$gene %in% gene,]
  } else{

    #If only chromosome/s provided, filter according to this
    if(!is.na(chr) & is.null(gene) & chr %in% chromosomes ) {
      ccr<-ccr[ccr$chrom %in% chr,]
      } else {

        #If only gene/s was/were provided, filter according to this
        if(is.na(chr) & !is.null(gene) ) {
          ccr<-ccr[ccr$gene %in% gene,]
        }
       }
    }



  if(nrow(ccr) == 0) {stop( paste("ERROR: ",gene," and/or ",chr,"IS/ARE NOT in the CCRs file" )) }

  return(ccr)
}


#getIDs----
#' @title getIDs
#' @description Gets the Ensembl-Uniprot IDs mapping, for those transcripts
#' that encode for proteins and that are part of the Gencode basic set of transcripts
#'
#' @param mapping_file The location of the gtf file for the actual chromosome
#' @param gtf_file The location of the gtf file for the actual chromosome
#' @param chr Chromosome
#' @param ens_version Ensembl version
#' @return A data frame with identifiers
#' @export
getIDs <- function( mapping_file,gtf_file,chr,ens_version) {

  cat("|-----------loading ID mappings file\n")


  #Loading The ENSEMBL - UniProtKb identifiers mapping file
  command <- paste("zcat ", mapping_file,sep="")
  map <- data.frame(data.table::fread(text = system(command, intern = T), sep = "\t"))
  map<-map[,c("ensembl_transcript_id","ensembl_protein_id","uniprot_blastp_id","ensp_uniprot_blastp_match","uniprot_blastp_type")]
  #names(map) <- c("transcript_id","uniprot_id")


  ## GET THE ENSEMBL GTF file ###########################
  # Make sure you splitted the original Homo_sapiens.GRCh##.###.gtf by chr
  # using the splitGTFbyChr.sh script

  #file <- paste(gtf_file, "Homo_sapiens.GRCh38.",ens_version,".chr",chr,".gtf.gz",sep="")
  gtf <- rtracklayer::import(con=gtf_file)
  gtf_df<-as.data.frame(gtf)

  # filter protein coding GENCODE_basic transcripts
  gencode_basic <- unique(subset(gtf_df,type=="CDS" & tag == "basic" & transcript_biotype == "protein_coding", select=c(seqnames,gene_name,gene_id,start,end,strand,exon_number,transcript_id,protein_id,tag)))
  gencode_basic_uniprot <- unique(subset(merge(x=gencode_basic,y=map,by.y="ensembl_transcript_id",by.x="transcript_id",all=T),seqnames==chr,select=c(seqnames,strand,gene_name,gene_id,transcript_id,protein_id,uniprot_blastp_id,tag,ensembl_protein_id,ensp_uniprot_blastp_match,uniprot_blastp_type)))
  gencode_basic_uniprot <- gencode_basic_uniprot[order(gencode_basic_uniprot$gene_id),]

  colnames(gencode_basic_uniprot)[names(gencode_basic_uniprot) == "transcript_id"] <- "ensembl_transcript_id"
  colnames(gencode_basic_uniprot)[names(gencode_basic_uniprot) == "protein_id"] <- "ensembl_protein_id"
  colnames(gencode_basic_uniprot)[names(gencode_basic_uniprot) == "gene_id"] <- "ensembl_gene_id"
  colnames(gencode_basic_uniprot)[names(gencode_basic_uniprot) == "gene_name"] <- "ensembl_gene_name"
  colnames(gencode_basic_uniprot)[names(gencode_basic_uniprot) == "seqnames"] <- "chr"

  return(gencode_basic_uniprot)
}


#getENSPseq----
#' @title getENSPseq
#' @description Downloads the Ensembl protein sequence for a ENSP identifier and
#' saves in a specified file
#'
#' @param ensembl_protein_id The Ensembl identifier (ENSP) of a protein sequence
#'  to download in fasta format
#' @param file TThe location of the file where to save the *.fasta sequence
#' @return Nothing
#' @export
getENSPseq<-function(ensembl_protein_id,file,gen_pad){

  cat(paste(gen_pad,"|-----------Obtaining FASTA sequence for ",ensembl_protein_id," from Ensembl\nSequence will be in ",file,"\n",sep=""))






  requestURL <- paste("https://rest.ensembl.org/sequence/id/",ensembl_protein_id,"?",sep="")

  r <- httr::GET(requestURL, httr::content_type("text/x-fasta"),charset = "UTF-8")
  seq<-NULL

  seq <- suppressMessages(httr::content(r,as = "text",charset = "UTF-8"))

  #only if properly download print out, otherwise the mapping to this protein will be skipped
  if(!(grepl("error|not found|exceed|Service Unavailable",seq)) && !purrr::is_empty(seq))
  {
    write(seq, file,append = T)
  }

}


#getUniPseq----
#' @title getUniPseq
#' @description Downloads the UniProtKB protein sequence for a given identifier and
#' saves in a specified file
#'
#' @param uniprot_blastp_ids A vector with UniProtKB identifier/s of a protein sequence/s.
#' to download in fasta format
#' @param file The location of the file where to save the *.fasta sequence
#' @param gencode_basic_uniprot Mapping table, related identifiers will be taken from here
#' @return Nothing
#' @export
getUniPseq<-function(uniprot_blastp_ids,gencode_basic_uniprot,file,gen_pad){

  #Avoid scientific notation, otherwise genomic positions with a big number might be printed out in scientific notation
  options(scipen=999)
  #Suppressing warnings
  options(warn=-1)

  for (uniprot_blastp_id in unique(uniprot_blastp_ids)){

  cat(paste(gen_pad,"|-----------Obtaining FASTA sequence for ",uniprot_blastp_id," from UniProtKB\nSequence will be in ",file,"\n",sep=""))

  # isoforms and canonical retrieval
  requestURL <- paste("https://www.ebi.ac.uk/proteins/api/proteins/",uniprot_blastp_id,"/isoforms",sep="")


  #DEBUG:
  #requestURL <- "https://www.ebi.ac.uk/proteins/api/proteins/K7EKL1/isoforms"

  r<-NULL
  r <- httr::GET(requestURL, httr::accept("text/x-fasta"))
  seq<-NULL
  seq <- suppressMessages(httr::content(r,as = "text"))

  # Usually works better for retrieving canonical or trEMBL sequences
  # if no sequences were retrieved when asking for isoforms, or if the code of the
  # canonical is not in the isoforms, try to make the request downloading only the canonical
  if(purrr::is_empty(seq) || grepl(paste("\\|",uniprot_blastp_id,"\\|",sep=""),seq) == F )# uniprot_blastp_id
  {
    #cat(">>>|-----------Downloading canonical or trEMBL sequence\n")
    requestURL <- paste("https://www.ebi.ac.uk/proteins/api/proteins/",uniprot_blastp_id,sep="")
    r<-NULL
    r <- httr::GET(requestURL, httr::accept("text/x-fasta"))
    seqaux<-suppressMessages(httr::content(r,as = "text"))

    type<-unique(gencode_basic_uniprot$uniprot_blastp_type[gencode_basic_uniprot$uniprot_blastp_id==uniprot_blastp_id])

    #get the correct uniprot id for the canonical sequence, as it is stated in the
    # uniprot plain text for the entry
    correct_unipid<-unique(gencode_basic_uniprot$uniprot_blastp_id[grepl(uniprot_blastp_id,gencode_basic_uniprot$uniprot_blastp_id) & gencode_basic_uniprot$uniprot_blastp_type %in% c("sp_c","tr")])


    correct_pattern<-paste(">",type,"\\|",correct_unipid,"\\|",sep="")
    toreplace<-paste(">",type,"\\|",uniprot_blastp_id,"\\|",sep="")

    #replace for the correct sp canonical uniprotid when it corresponds
    #e.g.:  Q96GD0 -> Q96GD0-1
    if(grepl(toreplace,seqaux) && !grepl(correct_pattern,seqaux))
    {seqaux<-stringr::str_replace(seqaux,toreplace,correct_pattern)}


    seq <- paste( seq, seqaux,sep="") #paste with content of previous retrieve
  }

  # write out sequences only if no Server error
  # and if at least a sequence was retrieved (either canonical, isoform or trEMBL),
  if(!purrr::is_empty(seq) && !(grepl("error|not found|exceed|Service Unavailable",seq))  )
  {
    write(seq, file,append = T)
  }
 }
}


#mapping----
#' @title mapping
#' @description Performs the mapping of genomic positions to ensembl protein
# positions, and adds the corresponding UniProtAcc identifier, positions and amino acids
#'
#' @param uniprot_blastp_id The UniProtKB identifier of a protein sequence
#'  to download in fasta format
#' @param fastafolder The location of the  *.fasta sequences
#' @param gencode_basic_uniprot Mapping table, related identifiers will be taken from here
#' @param gnm database connection, generated with with 'ensembldb'
#' @param bysite a dataframe with CCRs genomic regions transformed into one genomic
#' position by line
#' @param selectedtx the selected transcripts for the mapping
#' @param ens_version version of Ensembl
#' @param chr chromosome of interest
#' @return A dataframe with mapping of coordinates
#' @export
mapping <- function(gnm,bysite,selectedtx,gencode_basic_uniprot,blast_matches_seq,fastafolder,ens_version,chr){

  #Avoid scientific notation, otherwise genomic positions with a big number might be printed out in scientific notation
  options(scipen=999)
  #Suppressing warnings
  options(warn=-1)



  cat("|-----------Doing the mapping via the selected transcripts\n")

  # database connections need to be established and filtered in each child process,
  # they cannot be passed as a '.package' or 'function' argument when calling %dopar%
  #suppressPackageStartupMessages(library(EnsDb.Hsapiens.v101))
  #Calling functions defined here inside function that is running in child processes is not working,
  #It triggers : << task 1 failed - "could not find function "getEnsDB">>
  #edb<-getEnsDB(ens_version,chr)
  #This is the code of getEnsDB

  cat(paste("|-----------Loading Ensembl v",ens_version," for Human  (process ID: ",Sys.getpid(),"@",Sys.info()[['nodename']]," )\n",sep=""))

  # start an instance of AnnotationHub
  ah <- AnnotationHub::AnnotationHub()

  # This is how you query for all available EnsDb databases, to see what's available
  #query(ah, "EnsDb")

  # Query for Human available EnsDb database in a specific version
  ahDb <- AnnotationHub::query(ah, pattern = c("Homo sapiens", "EnsDb", ens_version))
  edb <- ahDb[[1]]

  # filter the database to the actual chromosome
  edb <- ensembldb::filter(edb, filter = ~ seq_name == chr)
  #End of the code of getEnsDB

  gnm_prt<-NULL
  mapped<-NULL
  result<-data.frame(uniprot_spc=NA,ensembl_uniprot_spc=NA,ensembl_uniprot_sp=NA,ensembl_uniprot_tr=NA, ensembl_uniprot_nm=NA)

  # select only the transcripts of interest in EnsDb, this will speedup searches/mappings
  edbxx <- ensembldb::filter(edb, filter = TxIdFilter(selectedtx$ensembl_transcript_id))

  # KEY ACTION: performs the mapping from genomic coordinates to amino acids
  # considering the selected transcripts (selectedtx)
  # that were filtered in the edbxx subset
  # TAKES A WHILE TO PROCESS, depending in the amount of genomic positions to map
  try(gnm_prt <- ensembldb::genomeToProtein(gnm, edbxx),silent=F)

  #DEBUG: WARNING only works when gnm_prt is empty, otherwise error:  "no method for coercing this S4 class to a vector"
  #message(paste("gnm_prt*",gnm_prt,sep=""))

  if(!is.null(gnm_prt)){
    unlisted_gnm_prt<-unlist(gnm_prt)

    mapped<-data.frame(cbind(unlisted_gnm_prt@elementMetadata@listData$seq_start,
                             ifelse(unlisted_gnm_prt@elementMetadata@listData$tx_id == "" | is.na(unlisted_gnm_prt@elementMetadata@listData$tx_id),NA,unlisted_gnm_prt@elementMetadata@listData$tx_id),
                             #ifelse(sjmisc::purrr::is_empty(unlist(unlisted_gnm_prt@NAMES)),NA,unlist(unlisted_gnm_prt@NAMES)),
                             ifelse(unlisted_gnm_prt@start==-1,NA,unlisted_gnm_prt@start) )
                       ,stringsAsFactors = FALSE)

    names(mapped)<-c("pos","ensembl_transcript_id","ensembl_aa_pos")


    #DEBUG: use this for for selecting smaller slices of 'mapped'
    small<-as.data.frame(mapped)
    names(small)<-c("pos","ensembl_transcript_id","ensembl_aa_pos")


    for(pos in unique(small$pos)){

      #obtain ENSG identifiers
      grf <- GRangesFilter(GRanges(unique(bysite$chrom), ranges = IRanges::IRanges(start=as.numeric(pos),width=1)), type = "any")
      gene_region<-genes(edbxx, filter = grf)


      if( !purrr::is_empty(small$ensembl_transcript_id[small$pos==pos]) )
      {

        # looping through each ENSG identifier
        for(gene_id in unique(gene_region$gene_id)) {
          #ENSG associated ENST
          txs <- as.data.frame(transcripts(edbxx, filter = GeneIdFilter(gene_id)))


          for(tx in txs$tx_id){
            if(tx %in% small$ensembl_transcript_id[small$pos==pos]) #if ENST associated to position
            {
              mapped$ensembl_gene_id[mapped$pos==pos & mapped$ensembl_transcript_id==tx]<-unique(txs$gene_id[txs$tx_id==tx])
              mapped$strand[mapped$pos==pos & mapped$ensembl_transcript_id==tx]<-as.character(unique(txs$strand[txs$tx_id==tx]))
            }
          }
        }
      } else if (is.na(unique(small$ensembl_transcript_id[small$pos==pos])))  { #end if the position has a ENS associated
        for(gene_id in unique(gene_region$gene_id)) {
          #region spanned by the ENSG gene
          gene_start<- gene_region@ranges@start[gene_region@ranges@NAMES==gene_id]
          gene_end<- gene_start + gene_region@ranges@width[gene_region@ranges@NAMES==gene_id] -1

          mapped$ensembl_gene_id[mapped$pos==pos & is.na(mapped$ensembl_transcript_id)]<-ifelse(pos>=gene_start && pos<=gene_end,gene_id,NA)
          mapped$strand[mapped$pos==pos & is.na(mapped$ensembl_transcript_id)]<-as.character(unique(gene_region@strand[gene_region@ranges@NAMES==gene_id]))

        }

      }
    }



    #creating a column of NAs for ensembl protein ids (ENSP)
    mapped$ensembl_protein_id <- NA

    #populating the ensembl protein ids (ENSP)
    aux<-unlisted_gnm_prt@NAMES
    aux[aux==0 | aux ==""] <- NA
    if(length(aux)>0)
    {
      mapped$ensembl_protein_id <- aux
    }


    #ensembl_protein_id and sequences
    ensp_seqs<-data.frame(ensembl_protein_id=NA,ensembl_seq=NA)

    for ( ensembl_protein_id in as.character(na.omit(unique(mapped$ensembl_protein_id))))
    {
      ensembl_seq<-blast_matches_seq$ensembl_seq[blast_matches_seq$ensembl_protein_id==ensembl_protein_id]

      # If it was not possible to obtain the sequence, the amino acids will be NA for the positions
      if(purrr::is_empty(ensembl_seq)) {
        ensembl_seq=NA
        cat(paste("Failed obtaining protein sequence for: ",ensembl_protein_id," it will not be in the mappings\t"))
      }

      line<-cbind(ensembl_protein_id,ensembl_seq)
      ensp_seqs<-rbind(ensp_seqs,line)
    }


    ensp_seqs<-unique(na.omit(ensp_seqs))

    # creating a column of NA for the amino acid letters
    mapped$ensembl_aa<-rep(NA,length(mapped$ensembl_aa_pos))

    #Populating the NA column with ENSP amino acids, obtaiend from the fasta sequences file
    mapped$ensembl_aa<-apply(mapped, MARGIN = 1,function(x) if(!is.na(x[3]) && x[3] != "NA" &&  !purrr::is_empty(substr(unique(ensp_seqs$ensembl_seq[ensp_seqs$ensembl_protein_id==x[6] ]), x[3], x[3])))
    {
      substr(unique(ensp_seqs$ensembl_seq[ensp_seqs$ensembl_protein_id==x[6] ]), x[3], x[3])
    }else{NA} #only when the ENSP was not found in the fasta file or amino acid positions where NA

    )

  }else{  #if the genomic positions don't map to any GENCODE_basic transcript
    mapped<-data.frame(pos=gnm@ranges@start)
    mapped$ensembl_transcript_id <-NA
    mapped$ensembl_protein_id <- NA
    mapped$ensembl_aa_pos <- NA
    mapped$ensembl_aa <- NA
    aux<-data.frame(pos=NA,ensembl_gene_id=NA,strand=NA)

    for(pos in unique(mapped$pos)){
      grf <- GRangesFilter(GRanges(unique(bysite$chrom), ranges = IRanges::IRanges(start=as.numeric(pos),width=1)), type = "any")
      gene_region<-genes(edbxx, filter = grf)

      if(!purrr::is_empty(unique(gene_region$gene_id))){

        for(gene_id in unique(gene_region$gene_id)) {
          # Region spanned by the ENSG
          gene_start<- gene_region@ranges@start[gene_region@ranges@NAMES==gene_id]
          gene_end<- gene_start + gene_region@ranges@width[gene_region@ranges@NAMES==gene_id] -1

          ensembl_gene_id<-ifelse(pos>=gene_start && pos<=gene_end,gene_id,NA)
          strand<-ifelse(!is.na(ensembl_gene_id),as.character(unique(gene_region@strand[gene_region@ranges@NAMES==gene_id])),NA)
          aux<-rbind(aux,cbind(pos,ensembl_gene_id,strand))
        }

      }else{
        ensembl_gene_id<-NA
        strand<-NA
        aux<-rbind(aux,cbind(pos,ensembl_gene_id,strand))
      }
    }

    mapped<-merge(x=mapped,y=aux,by=c("pos"),all.x=T)
  }


  mapped_tr<-data.frame(pos=NA,ensembl_gene_id.x=NA,strand.x=NA,ensembl_transcript_id=NA,ensembl_protein_id.x=NA,ensembl_protein_id.y=NA,ensembl_aa_pos=NA,ensembl_aa=NA,uniprot_blastp_id=NA)
  mapped_sp<-data.frame(pos=NA,ensembl_gene_id.x=NA,strand.x=NA,ensembl_transcript_id=NA,ensembl_protein_id.x=NA,ensembl_protein_id.y=NA,ensembl_aa_pos=NA,ensembl_aa=NA,uniprot_blastp_id=NA)
  mapped_spc<-data.frame(pos=NA,ensembl_gene_id.x=NA,strand.x=NA,ensembl_transcript_id=NA,ensembl_protein_id.x=NA,ensembl_protein_id.y=NA,ensembl_aa_pos=NA,ensembl_aa=NA,uniprot_blastp_id=NA)
  mapped_nm<-data.frame(pos=NA,ensembl_gene_id.x=NA,strand.x=NA,ensembl_transcript_id=NA,ensembl_protein_id.x=NA,ensembl_protein_id.y=NA,ensembl_aa_pos=NA,ensembl_aa=NA,uniprot_blastp_id=NA)


  #If there are informative rows out of merging the mapping and the gencode_basic_uniprot mapping, for each type of UniProt identifier

  #if(!purrr::is_empty(subset(merge(x=mapped,y=gencode_basic_uniprot,by=c("ensembl_transcript_id"),all.x=T),ensp_uniprot_blastp_match!="no_match" & uniprot_blastp_type=="tr" & ensembl_protein_id.x != "NA" & !is.na(ensembl_protein_id.x) ,select=c(chr,pos,ensembl_gene_id.x,strand.x,ensembl_transcript_id,ensembl_protein_id.x,ensembl_aa_pos,ensembl_aa,uniprot_blastp_id)))){
  if(nrow(subset(merge(x=mapped,y=gencode_basic_uniprot,by=c("ensembl_transcript_id"),all.x=T),ensp_uniprot_blastp_match!="no_match" & uniprot_blastp_type=="tr" & ensembl_protein_id.x != "NA" & !is.na(ensembl_protein_id.x) ,select=c(chr,pos,ensembl_gene_id.x,strand.x,ensembl_transcript_id,ensembl_protein_id.x,ensembl_aa_pos,ensembl_aa,uniprot_blastp_id)))){
      mapped_tr<-subset(merge(x=mapped,y=gencode_basic_uniprot,by=c("ensembl_transcript_id"),all.x=T),ensp_uniprot_blastp_match!="no_match" & uniprot_blastp_type=="tr" & ensembl_protein_id.x != "NA" & !is.na(ensembl_protein_id.x) ,select=c(chr,pos,ensembl_gene_id.x,strand.x,ensembl_transcript_id,ensembl_protein_id.x,ensembl_aa_pos,ensembl_aa,uniprot_blastp_id))}

  #if(!purrr::is_empty(subset(merge(x=mapped,y=gencode_basic_uniprot,by=c("ensembl_transcript_id"),all.x=T),ensp_uniprot_blastp_match!="no_match" & uniprot_blastp_type=="sp" & ensembl_protein_id.x != "NA" & !is.na(ensembl_protein_id.x) ,select=c(chr,pos,ensembl_gene_id.x,strand.x,ensembl_transcript_id,ensembl_protein_id.x,ensembl_aa_pos,ensembl_aa,uniprot_blastp_id)))){
  if(nrow(subset(merge(x=mapped,y=gencode_basic_uniprot,by=c("ensembl_transcript_id"),all.x=T),ensp_uniprot_blastp_match!="no_match" & uniprot_blastp_type=="sp" & ensembl_protein_id.x != "NA" & !is.na(ensembl_protein_id.x) ,select=c(chr,pos,ensembl_gene_id.x,strand.x,ensembl_transcript_id,ensembl_protein_id.x,ensembl_aa_pos,ensembl_aa,uniprot_blastp_id)))){
     mapped_sp<-subset(merge(x=mapped,y=gencode_basic_uniprot,by=c("ensembl_transcript_id"),all.x=T),ensp_uniprot_blastp_match!="no_match" & uniprot_blastp_type=="sp" & ensembl_protein_id.x != "NA" & !is.na(ensembl_protein_id.x) ,select=c(chr,pos,ensembl_gene_id.x,strand.x,ensembl_transcript_id,ensembl_protein_id.x,ensembl_aa_pos,ensembl_aa,uniprot_blastp_id))}

  #if(!purrr::is_empty(subset(merge(x=mapped,y=gencode_basic_uniprot,by=c("ensembl_transcript_id"),all.x=T),ensp_uniprot_blastp_match!="no_match" & uniprot_blastp_type=="sp_c" & ensembl_protein_id.x != "NA" & !is.na(ensembl_protein_id.x) ,select=c(chr,pos,ensembl_gene_id.x,strand.x,ensembl_transcript_id,ensembl_protein_id.x,ensembl_aa_pos,ensembl_aa,uniprot_blastp_id)))){
  if(nrow(subset(merge(x=mapped,y=gencode_basic_uniprot,by=c("ensembl_transcript_id"),all.x=T),ensp_uniprot_blastp_match!="no_match" & uniprot_blastp_type=="sp_c" & ensembl_protein_id.x != "NA" & !is.na(ensembl_protein_id.x) ,select=c(chr,pos,ensembl_gene_id.x,strand.x,ensembl_transcript_id,ensembl_protein_id.x,ensembl_aa_pos,ensembl_aa,uniprot_blastp_id)))){
    mapped_spc<-subset(merge(x=mapped,y=gencode_basic_uniprot,by=c("ensembl_transcript_id"),all.x=T),ensp_uniprot_blastp_match!="no_match" & uniprot_blastp_type=="sp_c" & ensembl_protein_id.x != "NA" & !is.na(ensembl_protein_id.x) ,select=c(chr,pos,ensembl_gene_id.x,strand.x,ensembl_transcript_id,ensembl_protein_id.x,ensembl_aa_pos,ensembl_aa,uniprot_blastp_id))}

  #if(!purrr::is_empty(subset(merge(x=mapped,y=gencode_basic_uniprot,by=c("ensembl_transcript_id"),all.x=T),ensp_uniprot_blastp_match=="no_match",select=c(chr,pos,ensembl_gene_id.x,strand.x,ensembl_transcript_id,ensembl_protein_id.x,ensembl_aa_pos,ensembl_aa,uniprot_blastp_id)))){
  if(nrow(subset(merge(x=mapped,y=gencode_basic_uniprot,by=c("ensembl_transcript_id"),all.x=T),ensp_uniprot_blastp_match=="no_match" ,select=c(chr,pos,ensembl_gene_id.x,strand.x,ensembl_transcript_id,ensembl_protein_id.x,ensembl_aa_pos,ensembl_aa,uniprot_blastp_id)))){
    mapped_nm<-subset(merge(x=mapped,y=gencode_basic_uniprot,by=c("ensembl_transcript_id"),all.x=T),ensp_uniprot_blastp_match=="no_match",select=c(chr,pos,ensembl_gene_id.x,strand.x,ensembl_transcript_id,ensembl_protein_id.x,ensembl_aa_pos,ensembl_aa,uniprot_blastp_id))}


  #For cases where more than one ENSG ID are assigned to one gene_name (HGNC Symbol) (e.g. TBCE in chr 1 has ENSG00000284770 and ENSG00000285053)
  #if(length(unique(gencode_basic_uniprot$gene_id[gencode_basic_uniprot$gene_name==unique(bysite$gene)])) == 1){
  #mapped$gene_id<-unique(gencode_basic_uniprot$gene_id[gencode_basic_uniprot$gene_name==unique(bysite$gene)] )
  #}else{
  # geneid_tx<-unique(subset(gencode_basic_uniprot, gene_name==bysite$gene,select = c(gene_id,transcript_id)))
  #mapped$gene_id<-paste(unique(geneid_tx$gene_id[geneid_tx$transcript_id==mapped$transcript_id]),collapse = ",")
  #mapped$gene_id<-paste(unique(gencode_basic_uniprot$gene_id[gencode_basic_uniprot$gene_name==unique(bysite$gene) & gencode_basic_uniprot$transcript_id[gencode_basic_uniprot$transcript_id==mapped$transcript_id] ]),collapse=",")

  #}




  #Renaming columns. Merge adds ".x", ".y" etc when the merged tables have columns
  #with same name
  mapped_tr <- mapped_tr %>% dplyr::rename (ensembl_protein_id=ensembl_protein_id.x)
  mapped_sp <- mapped_sp %>% dplyr::rename (ensembl_protein_id=ensembl_protein_id.x)
  mapped_spc <- mapped_spc %>% dplyr::rename (ensembl_protein_id=ensembl_protein_id.x)
  mapped_nm <- mapped_nm %>% dplyr::rename (ensembl_protein_id=ensembl_protein_id.x)


  mapped_tr <- mapped_tr %>% dplyr::rename (ensembl_gene_id=ensembl_gene_id.x)
  mapped_sp <- mapped_sp %>% dplyr::rename (ensembl_gene_id=ensembl_gene_id.x)
  mapped_spc <- mapped_spc %>% dplyr::rename (ensembl_gene_id=ensembl_gene_id.x)
  mapped_nm <- mapped_nm %>% dplyr::rename (ensembl_gene_id=ensembl_gene_id.x)


  mapped_tr <- mapped_tr %>% dplyr::rename (strand=strand.x)
  mapped_sp <- mapped_sp %>% dplyr::rename (strand=strand.x)
  mapped_spc <- mapped_spc %>% dplyr::rename (strand=strand.x)
  mapped_nm <- mapped_nm %>% dplyr::rename (strand=strand.x)


  #if the genomic position did not mapped to a ensembl protein position (NAs or empty), then the uniprot_id is also turned into NA
  mapped_tr$uniprot_blastp_id[is.na(mapped_tr$ensembl_protein_id) & (is.na(mapped_tr$ensembl_aa_pos) | length(mapped_tr$ensembl_aa_pos) == 0)] <- NA
  mapped_sp$uniprot_blastp_id[is.na(mapped_sp$ensembl_protein_id) & (is.na(mapped_sp$ensembl_aa_pos) | length(mapped_sp$ensembl_aa_pos) == 0)] <- NA
  mapped_spc$uniprot_blastp_id[is.na(mapped_spc$ensembl_protein_id) & (is.na(mapped_spc$ensembl_aa_pos) | length(mapped_spc$ensembl_aa_pos) == 0)] <- NA
  mapped_nm$uniprot_blastp_id[is.na(mapped_nm$ensembl_protein_id) & (is.na(mapped_nm$ensembl_aa_pos) | length(mapped_nm$ensembl_aa_pos) == 0)] <- NA

  mapped_tr <- unique(mapped_tr)
  mapped_spc <- unique(mapped_spc)
  mapped_sp <- unique(mapped_sp)
  mapped_nm <- unique(mapped_nm)


  cat("Adding amino acids from UniProtKB sequences...\n")

  mapped_tr$unipaa_pos<-NA
  mapped_tr$unipaa_pos<-apply(mapped_tr, MARGIN = 1,function(x) if(!is.na(x[7]) && x[7] != "NA" && !is.na(x[9]) && x[9] != "NA" ) {x[7]}else{NA})

  mapped_nm$unipaa_pos<-NA


  mapped_sp$unipaa_pos<-NA
  mapped_sp$unipaa_pos<-apply(mapped_sp, MARGIN = 1,function(x) if(!is.na(x[7]) && x[7] != "NA" && !is.na(x[9]) && x[9] != "NA" ) {x[7]}else{NA})


  mapped_spc$unipaa_pos<-NA
  mapped_spc$unipaa_pos<-apply(mapped_spc, MARGIN = 1,function(x) if(!is.na(x[7]) && x[7] != "NA" && !is.na(x[9]) && x[9] != "NA" ) {x[7]}else{NA})



  mapped_tr$unipaa<-NA
  mapped_tr$unipaa<-apply(mapped_tr, MARGIN = 1,function(x) if(!is.na(x[10]) && x[10] != "NA" &&  !purrr::is_empty(substr(unique(blast_matches_seq$uniprot_seq[blast_matches_seq$uniprot_blastp_id==x[9]]), x[10], x[10])))
  {unique(substr(unique((blast_matches_seq$uniprot_seq[blast_matches_seq$uniprot_blastp_id==x[9] ])), x[10], x[10]))}else{NA})


  mapped_nm$unipaa<-NA

  mapped_sp$unipaa<-NA
  mapped_sp$unipaa<-apply(mapped_sp, MARGIN = 1,function(x) if(!is.na(x[10]) && x[10] != "NA" &&  !purrr::is_empty(substr(unique(blast_matches_seq$uniprot_seq[blast_matches_seq$uniprot_blastp_id==x[9]]), x[10], x[10])))
  {unique(substr(unique((blast_matches_seq$uniprot_seq[blast_matches_seq$uniprot_blastp_id==x[9] ])), x[10], x[10]))}else{NA})


  mapped_spc$unipaa<-NA
  mapped_spc$unipaa<-apply(mapped_spc, MARGIN = 1,function(x) if(!is.na(x[10]) && x[10] != "NA" &&  !purrr::is_empty(substr(unique(blast_matches_seq$uniprot_seq[blast_matches_seq$uniprot_blastp_id==x[9] & blast_matches_seq$ensembl_protein_id==x[6]]), x[10], x[10])))
  {unique(substr(unique((blast_matches_seq$uniprot_seq[blast_matches_seq$uniprot_blastp_id==x[9] & blast_matches_seq$ensembl_protein_id==x[6]  ])), x[10], x[10]))}else{NA})


  #creating a single "map" column, colapsing ENST,ENSP,aa_pos,uniprot_id, "|" separated, for each position
  if(!purrr::is_empty(unique(mapped_tr$pos))){
    gnm_prt_pos_colapsed_tr<-data.table::setDT(mapped_tr)[,.(map_tr=paste(ensembl_transcript_id,ensembl_protein_id,ensembl_aa_pos,ensembl_aa,uniprot_blastp_id,unipaa_pos,unipaa,collapse ="|",sep=",")),by=.(chr,pos,strand,ensembl_gene_id)]}else
    {gnm_prt_pos_colapsed_tr<-data.table::data.table(chr=NA,pos=NA,strand=NA,ensembl_gene_id=NA,map_tr="NA,NA,NA,NA,NA,NA,NA")}


  if(!purrr::is_empty(unique(mapped_sp$pos))){
    gnm_prt_pos_colapsed_sp<-data.table::setDT(mapped_sp)[,.(map_sp=paste(ensembl_transcript_id,ensembl_protein_id,ensembl_aa_pos,ensembl_aa,uniprot_blastp_id,unipaa_pos,unipaa,collapse ="|",sep=",")),by=.(chr,pos,strand,ensembl_gene_id)]}else
    {gnm_prt_pos_colapsed_sp<-data.table::data.table(chr=NA,pos=NA,strand=NA,ensembl_gene_id=NA,map_sp="NA,NA,NA,NA,NA,NA,NA")}


  if(!purrr::is_empty(unique(mapped_spc$pos))){
    gnm_prt_pos_colapsed_spc<-data.table::setDT(mapped_spc)[,.(map_spc=paste(ensembl_transcript_id,ensembl_protein_id,ensembl_aa_pos,ensembl_aa,uniprot_blastp_id,unipaa_pos,unipaa,collapse ="|",sep=",")),by=.(chr,pos,strand,ensembl_gene_id)]}else
    {gnm_prt_pos_colapsed_spc<-data.table::data.table(chr=NA,pos=NA,strand=NA,ensembl_gene_id=NA,map_spc="NA,NA,NA,NA,NA,NA,NA")}


  if(!purrr::is_empty(unique(mapped_nm$pos))){
    gnm_prt_pos_colapsed_nm<-data.table::setDT(mapped_nm)[,.(map_nm=paste(ensembl_transcript_id,ensembl_protein_id,ensembl_aa_pos,ensembl_aa,uniprot_blastp_id,unipaa_pos,unipaa,collapse ="|",sep=",")),by=.(chr,pos,strand,ensembl_gene_id)]}else
    {gnm_prt_pos_colapsed_nm<-data.table::data.table(chr=NA,pos=NA,strand=NA,ensembl_gene_id=NA,map_nm="NA,NA,NA,NA,NA,NA,NA")}



  if(!purrr::is_empty(unique(mapped_spc$pos)) & length(na.omit(mapped_spc$pos))!=0){
    simple_spc<-unique((mapped_spc %>% dplyr::select(chr,pos,strand,ensembl_gene_id,uniprot_blastp_id,unipaa_pos,unipaa)))}else
    {simple_spc<-data.table::data.table(chr=NA,pos=NA,strand=NA,ensembl_gene_id=NA,uniprot_blastp_id=NA,unipaa_pos=NA,unipaa=NA)}



  if(!purrr::is_empty(unique(simple_spc$pos))){
    gnm_prt_pos_colapsed_simple_spc<-data.table::setDT(simple_spc)[,.(simple_spc=paste(uniprot_blastp_id,unipaa_pos,unipaa,collapse ="|",sep=",")),by=.(chr,pos,strand,ensembl_gene_id)]}else
    {gnm_prt_pos_colapsed_simple_spc<-data.table::data.table(chr=NA,pos=NA,strand=NA,ensembl_gene_id=NA,simple_spc="NA,NA,NA")}



  gnm_prt_pos_colapsed_tr$pos <- sapply(gnm_prt_pos_colapsed_tr$pos,as.numeric)
  gnm_prt_pos_colapsed_sp$pos <- sapply(gnm_prt_pos_colapsed_sp$pos,as.numeric)
  gnm_prt_pos_colapsed_spc$pos <- sapply(gnm_prt_pos_colapsed_spc$pos,as.numeric)
  gnm_prt_pos_colapsed_nm$pos <- sapply(gnm_prt_pos_colapsed_nm$pos,as.numeric)
  gnm_prt_pos_colapsed_simple_spc$pos <- sapply(gnm_prt_pos_colapsed_simple_spc$pos,as.numeric)


  mapped$pos <- sapply(mapped$pos,as.numeric)

  #DEBUG:
  # cat(paste("mapped",  Sys.getpid(),"@",Sys.info()[['nodename']], "\n"))
  # message(colnames(mapped))
  # message("A small sample...")
  # message(mapped[1:15,])


  bysite<-subset(merge(x=bysite,y=mapped,by=c("pos"),all=T,allow.cartesian=TRUE),select=c(chrom,pos,gene,strand,ensembl_gene_id,weighted_pct,chrom_start_end_gen))
  bysite$pos <- sapply(bysite$pos,as.numeric)

  #DEBUG:
  #cat(paste("Merging mapped sites with CCRs information (process ID: ",Sys.getpid(),"@",Sys.info()[['nodename']]," )\n"))

  #merge of the colapsed mapping and the bysite information from the extended ccrs regions
  if(!is.na(unique(gnm_prt_pos_colapsed_tr$chr))){
    result<-subset(merge(x=gnm_prt_pos_colapsed_tr,y=bysite,by=c("pos","ensembl_gene_id"),all.y=T),select=c(chrom,pos,gene,strand.y,ensembl_gene_id,map_tr,weighted_pct,chrom_start_end_gen))
    result$pos <- sapply(result$pos,as.numeric)}else
    {
      result<-bysite %>% dplyr::select(chrom,pos,strand,gene,ensembl_gene_id,weighted_pct,chrom_start_end_gen)
      result$strand<-rep(unique(gencode_basic_uniprot$strand[gencode_basic_uniprot$ensembl_gene_name==unique(bysite$gene)],length(result$pos)))
      result$map_tr<-rep(gnm_prt_pos_colapsed_tr$map_tr,length(result$pos))
    }


  if(!is.na(unique(gnm_prt_pos_colapsed_sp$chr))){
    if("ensembl_gene_id"  %in% colnames(result)) {result_sp<-subset(merge(x=gnm_prt_pos_colapsed_sp,y=result,by=c("pos","ensembl_gene_id"),all.y=T),select=c(chrom,pos,gene,strand.y,ensembl_gene_id,map_sp,map_tr,weighted_pct,chrom_start_end_gen))}else{
      result_sp<-subset(merge(x=gnm_prt_pos_colapsed_sp,y=result,by=c("pos"),all.y=T),select=c(chrom,pos,gene,strand.y,ensembl_gene_id,map_sp,map_tr,weighted_pct,chrom_start_end_gen))

    }

    result_sp$pos <- sapply(result_sp$pos,as.numeric)
  }else{
    result_sp<-result
    result_sp$map_sp<-rep(gnm_prt_pos_colapsed_sp$map_sp,length(result_sp$pos))
  }



  if(!is.na(unique(gnm_prt_pos_colapsed_spc$chr))){
    if("ensembl_gene_id"  %in% colnames(result_sp)) {result_spc<-subset(merge(x=gnm_prt_pos_colapsed_spc,y=result_sp,by=c("pos","ensembl_gene_id"),all.y=T),select=c(chrom,pos,gene,strand.y,ensembl_gene_id,map_spc,map_sp,map_tr,weighted_pct,chrom_start_end_gen)) }else{
      result_spc<-subset(merge(x=gnm_prt_pos_colapsed_spc,y=result_sp,by=c("pos"),all.y=T),select=c(chrom,pos,gene,strand.y,ensembl_gene_id,map_spc,map_sp,map_tr,weighted_pct,chrom_start_end_gen))
    }
    result_spc$pos <- sapply(result_spc$pos,as.numeric)}else
    {
      result_spc<-result_sp
      result_spc$map_spc<-rep(gnm_prt_pos_colapsed_spc$map_spc,length(result_spc$pos))
    }



  if(!is.na(unique(gnm_prt_pos_colapsed_nm$chr))){
    if("ensembl_gene_id"  %in% colnames(result_spc)) {result_nm<-subset(merge(x=gnm_prt_pos_colapsed_nm,y=result_spc,by=c("pos","ensembl_gene_id"),all.y=T),select=c(chrom,pos,gene,strand.y,ensembl_gene_id,map_spc,map_sp,map_tr,map_nm,weighted_pct,chrom_start_end_gen))}else{
      result_nm<-subset(merge(x=gnm_prt_pos_colapsed_nm,y=result_spc,by=c("pos"),all.y=T),select=c(chrom,pos,gene,strand.y,ensembl_gene_id,map_spc,map_sp,map_tr,map_nm,weighted_pct,chrom_start_end_gen))
    }
    result_nm$pos <- sapply(result_nm$pos,as.numeric)}else
    {
      result_nm<-result_spc
      result_nm$map_nm<-rep(gnm_prt_pos_colapsed_nm$map_nm,length(result_nm$pos))
    }



  if(!is.na(unique(gnm_prt_pos_colapsed_simple_spc$chr))){
    if("ensembl_gene_id"  %in% colnames(result_nm)) {
      result_simple_spc<-subset(merge(x=gnm_prt_pos_colapsed_simple_spc,y=result_nm,by=c("pos","ensembl_gene_id"),all.y=T,allow.cartesian=T),select=c(chrom,pos,gene,strand.y,ensembl_gene_id,simple_spc,map_spc,map_sp,map_tr,map_nm,weighted_pct,chrom_start_end_gen))
    }else{
      result_simple_spc<-subset(merge(x=gnm_prt_pos_colapsed_simple_spc,y=result_nm,by=c("pos"),all.y=T,allow.cartesian=T),select=c(chrom,pos,gene,strand.y,ensembl_gene_id,simple_spc,map_spc,map_sp,map_tr,map_nm,weighted_pct,chrom_start_end_gen))
    }
    result_simple_spc$pos <- sapply(result_simple_spc$pos,as.numeric)}else
    {
      result_simple_spc<-result_nm
      result_simple_spc$simple_spc<-rep(gnm_prt_pos_colapsed_simple_spc$simple_spc,length(result_simple_spc$pos))
    }



  final_result<-result_simple_spc

  final_result$pos <- sapply(final_result$pos,as.numeric)

  if (  "strand.y"  %in% names(final_result)) {final_result<-final_result  %>% dplyr::rename (strand=strand.y)  }
  if (  "ensembl_gene_id.y"  %in% names(final_result))  {final_result<-final_result  %>% dplyr::rename (ensembl_gene_id=ensembl_gene_id.y)}
  if (  !("ensembl_gene_id"  %in% names(final_result)))  {final_result$ensembl_gene_id<-rep(NA,length(final_result$pos))}

  final_result$simple_spc[is.na(final_result$simple_spc)]<-"NA,NA,NA"
  final_result$map_spc[is.na(final_result$map_spc)]<-"NA,NA,NA,NA,NA,NA,NA"
  final_result$map_nm[is.na(final_result$map_nm)]<-"NA,NA,NA,NA,NA,NA,NA"
  final_result$map_sp[is.na(final_result$map_sp)]<-"NA,NA,NA,NA,NA,NA,NA"
  final_result$map_tr[is.na(final_result$map_tr)]<-"NA,NA,NA,NA,NA,NA,NA"


  #order the columns
  final_result<-unique(final_result[,c("chrom","pos","gene","strand","ensembl_gene_id","weighted_pct","simple_spc","map_spc","map_sp","map_tr","map_nm","chrom_start_end_gen")])

  cat("Returning mapped CCRs (process ID: ",Sys.getpid(),"@",Sys.info()[['nodename']]," )")
  return(final_result)
  #return(gnm_prt_pos_colapsed)
}



#matchedProteinSeqs----
#' @title matchedProteinSeqs
#' @description Gets the amino acid sequences and correspondences between ENSP and UniProt
#' Ensembl and UniProtKB proteins may differ in sequence even when they are assigned
#' as corresponding in these databases
#' This function runs blastp for comparing Ensembl and UniProtKB protein sequences
#' and makes sure they are 100% identical
#'
#' @param ensembl_gene_ids The UniProtKB identifier of a protein sequence
#'  to download in fasta format
#' @param gencode_basic_uniprot Mapping table, related identifiers will be taken from here
#' @param fastafolder will check if the protein sequences are already here, if not donwload
#' @param ens_version version of Ensembl
#' @return A dataframe
#' @export
matchedProteinSeqs <- function(ensembl_gene_ids,gencode_basic_uniprot,fastafolder,remove,gen_pad){

  cat(paste(gen_pad,"|-----------Running Blastp to check Ensembl vs UniProtKB protein sequences identity\n",sep=""))



  blast_matches_seq<-data.frame(uniprot_blastp_id=NA,ensembl_protein_id=NA,uniprot_seq=NA,ensembl_seq=NA)


  ##TODO: foreach in dopar parallel, although there is a limited amount of downloads which can be done from
  # the the UniProt API under the same user IP

  #For each gene identifier
  for (ensembl_gene_id in ensembl_gene_ids){

    #DEBUG:
    #ensembl_gene_id<-ensembl_gene_ids[1]

    query_file<-paste(fastafolder,ensembl_gene_id,"_unip.fasta",sep="") #queries are the UniProtSP sequences
    subject_file<-paste(fastafolder,ensembl_gene_id,"_ENSP.fasta",sep="") #subjects are the Ensembl protein sequences


    for(ensembl_protein_id in gencode_basic_uniprot$ensembl_protein_id[gencode_basic_uniprot$ensembl_gene_id==ensembl_gene_id])
    {



      # Obtain the UniProtKB identifiers associated with the actual ENSP identifier
      uniprot_blastp_id<-unique(gencode_basic_uniprot$uniprot_blastp_id[gencode_basic_uniprot$ensembl_protein_id==ensembl_protein_id &
                                                                   gencode_basic_uniprot$uniprot_blastp_id != "no_match"])

      #If there was a match wit uniprot
      if(!purrr::is_empty(uniprot_blastp_id) & length(na.omit(uniprot_blastp_id)) != 0 ){
        # Check whether the fasta sequence for the ENSP and uniprot proteins are present in 'data/mapping_tables/fastas/(chromosome_number)/'
        # if not, attempt to download the sequence and save it to the corresponding file
        if ( !file.exists(subject_file)  || purrr::is_empty(unique(grep(pattern = ensembl_protein_id,readLines(subject_file))))){ getENSPseq(ensembl_protein_id,subject_file,gen_pad)  }
        if ( !file.exists(query_file)  || purrr::is_empty(unique(grep(pattern = uniprot_blastp_id,readLines(query_file))))){ getUniPseq(uniprot_blastp_id,gencode_basic_uniprot,query_file,gen_pad)  }



        # Obtaining the sequence from the multifasta file, using bash command lines
        #comm<-paste("grep -zoP '>",ensembl_protein_id,".*[\\nA-Z]*[\\Z|>]*' ",subject_file ," | tr -d '>' | tail -n +2 | tr -d '\\n' ",sep="")
        #ensembl_seq <- system(comm, wait=T, intern = T)

      }
    }

    #Performing a blastp alignment
    format = '"6 qseqid sseqid pident qcovs qlen slen qstart qend sstart send qseq sseq"' #output format from blastp

    #command line for running blastp
    comm<-paste("blastp -query ",query_file," -subject ",subject_file," -outfmt ",format,
                "| awk '$4==100 && $5==$6 && $11!~/-/ && $12!~/-/ {print $1,$2,$3,$11,$12}' | tr '|' ' ' | sort -u |awk '{printf ", "\"%s\\t%s\\t%s\\t%s\\n\"",", $2,$4,$6,$7}' | awk -F'\\t' -v OFS='\\t' '{gsub(/\\.[0-9]*/,\"\",$2)} 1'",sep="")

    colnames<-c("uniprot_blastp_id","ensembl_protein_id","uniprot_seq","ensembl_seq")


    #calling to run blast by command line
     # aux <- system(comm, wait=T, intern = T) %>%
     #   dplyr::as_tibble() %>%
     #   tidyr::separate(col = value,
     #            into = colnames,
     #            sep = "\t",
     #            convert = TRUE)
     aux<- system(comm, wait=T, intern = T)
     aux<-dplyr::as_tibble(aux)
     aux<-tidyr::separate(aux, col = value,
                           into = colnames,
                           sep = "\t",
                           convert = TRUE)


    blast_matches_seq<-rbind(blast_matches_seq,aux)

  }
  blast_matches_seq<-na.omit(blast_matches_seq)

  # If removal of *.fasta files was requested when calling matchedProteinSeqs()
  if(remove){
    #Check existence of file and remove
    if (file.exists(subject_file)) {  file.remove(subject_file) }
    if (file.exists(query_file)) { file.remove(query_file) }
  }

  return(blast_matches_seq)
}


#byGene----
#' @title byGene
#' @description Goes gene-by-gene in your list, calling the mapping of positions
#' This function will do parallel calls to 'np' processors, dividing 'bysite' in chunks
#' among them
#'
#'
#' @param genes Your list of genes
#' @param gencode_basic_uniprot Mapping table, Ensembl and UniProtKB identifiers will be taken from here
#' @param fastafolder will check if the protein sequences are already here, if not donwload
#' @param ens_version version of Ensembl
#' @param chr chromosome of interest
#' @param edb a ensembldb database connection
#' @param ccr the CCRs dataframe with sites, CCRs percentiles, etc
#' @param outfile the output file
#' @param np number of proceses to run in parallel
#' @return A dataframe with final mapping
#' @importFrom foreach %dopar%
#' @importFrom usethis use_pipe
#' @export
byGene <- function(genes,gencode_basic_uniprot,edb,ccr,outfile,chr,fastafolder,ens_version,np){
  #usethis::use_pipe()
  #Free up some memory
  gc()

  #Avoid scientific notation, otherwise genomic positions with a big number might be printed out in scientific notation
  options(scipen=999)
  #Suppressing warnings
  options(warn=-1)



  cat(paste("|-----------Going through your list of genes:","\nChromosome: ",chr,"\nGenes: ",paste(unique(genes$gene),collapse=",") , "\n\n"))

  for ( gen in unique(genes$gene))
  {
    #DEBUG:
    #gen<-"VAPA"

    cat(paste(gen," ...\n"))

    #This is just to print the name of your genes while they are being processed
    gen_pad<-stringr::str_pad(gen,width=15, side = "right") #width=15, because the longest gene name has 15 characters

    #get the list of GENCODE_basic transcripts for this gene
    selectedtx<-unique(subset(gencode_basic_uniprot, ensembl_gene_name==gen,select = c(ensembl_gene_id,ensembl_transcript_id))) #"ENST00000342175"

    # calling to perform a blast alignment to get the correspondence between ENSP and UniProt sequences
    # This will also download the protein sequences into data/mapping_tables/fastas/(chromosome number)/*.fasta
    blast_matches_seq<-matchedProteinSeqs(unique(selectedtx$ensembl_gene_id),gencode_basic_uniprot,fastafolder,remove=F,gen_pad)

    #filter the actual chromosome to only the GENCODE_basic transcripts of interest
    edbxx <- ensembldb::filter(edb, filter = TxIdFilter(selectedtx$ensembl_transcript_id))


    geneccrsdf<-NULL
    #take the ccrs for this gen
    #geneccrsdf <-subset(ccr, gene == gene)
    geneccrsdf <- ccr[ccr$gene==gen,]


    #TODO: DO NOT COLLAPSE RANGES!! or the mapping function will try to map to intronic regions which will produce empty dataframes!!!
    #geneccrsdf$start <- sapply(geneccrsdf$ranges, function(x) strsplit(strsplit(as.character(x), ",")[[1]][1], "-")[[1]][1])

    ## +1 to start as ccrs starts follow the 0-based genomic coordinates system
    geneccrsdf$start <- sapply (geneccrsdf$start, function(x) as.numeric(x)+1)

    geneccrsdf <- as.data.frame(unique(geneccrsdf))


    #creating a unique ID per ccr region
    geneccrsdf$chrom_start_end_gen<-paste(geneccrsdf$chrom,geneccrsdf$start,geneccrsdf$end,geneccrsdf$gene,sep = "-")

    #select all/some of the regions for the gene, useful for DEBUG
    abit<-as.data.frame(geneccrsdf)


    #extending positions in ranges
    bysitedf<-NULL

    abit <- abit %>% dplyr::rename (s=start)
    abit <- abit %>% dplyr::rename (e=end)

    #IMPORTANT
    #Not adding strand directly from the mapping table, strand and all the gene
    #information will be added later by using ensembldb mapping tools
    #This is because a genesymbol can have N ensembl_gene_id and even be in different strands :( (e.g. TMSB15B)


    bysitedf<-abit %>% dplyr::mutate(pos=purrr::map2(s,e,~seq(.x,.y,by=1))) %>% tidyr::unnest(pos)

    #logfile<-paste("log_chr",chr,"_",gsub(":","-",gsub(" ","_",Sys.time())),".txt",sep="")

    #comment this for running in serial mode (not parallel)
    #outfile will be where all the messages from the child processes will be dumped
    cl <- parallel::makeCluster(np, outfile="/dev/null", setup_strategy = "sequential")
    doParallel::registerDoParallel(cl)

    #writeLines(c(""), logfile)
    results<- NULL
    pos_col<-which(colnames(bysitedf)=="pos")
    cat(paste(gen_pad,"|-----------Doing the mapping in parallel using ", np," cores.... (be patient)\n",sep=""))
    system.time(
      #sending chunks of lines, only works for %dopar%.
      #Sending chunks avoids having to load and filter the EnsDb.Hsapiens.v## so many times, as it would be done by sending one line at the time
      results <- foreach::foreach(chunk=itertools::isplitRows(bysitedf,chunks=np), .packages = c('dplyr','AnnotationHub','ensembldb','R.utils','sjmisc','data.table'), .export=c('mapping'), .combine=rbind) %dopar% {

        #Creating the genomic ranges that will be mapped to amino acids
        gnm <- GenomicRanges::GRanges(chr, IRanges::IRanges(start = as.numeric(unlist(chunk[,pos_col])), width = 1))

        #Dump into a logfile the messages from functions called in parallel
        #sink(logfile, append=TRUE)

        #Call to map the genomic ranges to amino acids
        suppressWarnings(suppressPackageStartupMessages(mapping(gnm,chunk,selectedtx,gencode_basic_uniprot,blast_matches_seq,fastafolder,ens_version,chr)))
      }
    )

    #comment this if running in serial
    parallel::stopCluster(cl)

    #FINALLY dump the results of this gene to a per-chr file
    #printout gene progress
    cat(paste(gen_pad,"|-----------OK!\n",sep=""))
  } #end of going gene by gene
}

