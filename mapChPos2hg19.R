# =======================================
# Author: Olga Nikolova
# E-mail: olga.nikolova@gmail.com
# =======================================
# Map chromosome positions to genes
# Map junction coordinates to exons
# =======================================
library(biomaRt)
options(stringsAsFactors=FALSE)

# Create map
ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl")
map <- getBM(c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position"), mart = ensembl)
map$start_position <- as.numeric(map$start_position)
map$end_position <- as.numeric(map$end_position)

# List of submaps per chrom
mapchr <- lapply(as.list(c(1:22,"X","Y")),function(x){
  return(map[map$chromosome_name == x,])
})
names(mapchr) <- c(1:22,"X","Y")

# Retreive ensembl_id and hugo_symbol for 
# given chromosome and position (positive, relative)
# Output: vector of 2 strings formatted:
#     "<E1>,<E2>,...,<En>;<H1>,<H2>,...,<Hn>" "<integer>"
#   where all ensembl IDs are separated by ","
#   and all hugo symbols are separated by ","
#   and the ensembl IDs are separated from hugo symbols by ";"
# Note: ensembl IDs and hugo symbols are returned paired
# i.e. if H2 is missing we have: "E1,E2,E3,..;H1,,H3,.."
getGeneID <- function(pos, map){
  if(is.null(map)){return("NA")}
  if( !is.numeric(pos) )
    stop("pos should be numeric")
  if( !is.data.frame(map) )
    stop("map should be a data frame")
  
  g <- map[intersect(which(map$start_position <= pos), 
                   which(map$end_position >= pos)),]

 if( nrow(g) > 0){
#    paste(c("A", "C", "T"), collapse=",")
#    g_e <- paste(g$ensembl_gene_id, collapse=",")
#    g_h <- paste(g$hgnc_symbol, collapse=",")
   
    g_e <- toString(g$ensembl_gene_id, sep=",")
    g_h <- toString(g$hgnc_symbol, sep=",")
    g_ids <- paste(g_e, g_h, sep=";")
    g_l <- toString(g$end_position - g$start_position, sep=",") 

#    return(cbind(g_ids, g_l))
    return(c(g_ids, g_l))
 }else {
    return(c("NA", "NA"))
  } 
}

outputnames <- c("SAGE_sampleID","SAGE_cancerType","SAGE_sequencingID","SAGE_sequencingLane","SAGE_sequencingBarcode",
                 "SAGE_donorGeneID", "SAGE_acceptorGeneID","SAGE_donorGeneLength","SAGE_acceptorGeneLength","SAGE_donorChromosome",
                 "SAGE_acceptorChromosome","donorEnd","acceptorStart","name", "coverage","strand","itemRgb","blockCount",
                 "blockSizes","blockStarts","entropy","flank_string_case", "flank_string","min_mismatch","max_mismatch",
                 "average_mismatch","maximal_of_minimal_acceptor_site_length",
                 "maximal of minimal acceptor site length","minimal_anchor_difference","unique_read_count",
                 "multiple_reads_count","paired_reads_count","left_paired_reads_count","right_paired_reads_count",
                 "multiple_paired_reads_count","unique_paired_reads_count","single_reads_count",
                 "encompassing_read_pairs_count","donorStart","acceptorEnd","donor_isoform_structures",
                 "acceptor_isoform_structures","donor_uniformity_score","acceptor_uniformity_score",
                 "donor_uniformity_KS-test_score","acceptor_uniformity_KS-test_score","minimal_donor_isoform_length",
                 "maximal_donor_isoform_length","minimal_acceptor_isoform_length","maximal_acceptor_isoform_length",
                 "donor_site_match_to_normal_junction","acceptor_site_match_to_normal_junction","donor_sequence",
                 "acceptor_sequence","match_to_gene_strand","fusion_source","fusion_type","gene_strand",
                 "donor_annotated_gene","acceptor_annotated_gene")

clist <- as.list(c("BLCA","CESC","BRCA","HNSC","KIRP","KICH","KIRC","LGG","LIHC","LUAD","LUSC","PAAD","PRAD",
                   "SKCM","THCA"))
#clist <- as.list(c("PRAD"))
ctnames <- unlist(clist)
# input directory
path <- "/Shares/work/DAT_109__TCGA_unc_rnaseq_genefusion/Data/mapped2TCGA_final/"

d.clist <- lapply(clist, function(type){
  name <- paste(path, type,".txt", sep="")
  return(read.csv(name, header=T, sep="\t"))#, as.is=TRUE))
})

for(ct in 1:length(d.clist)){
  df <- d.clist[[ct]]
  # Get IDs of donor genes
#   df$chromosomenum <- substring(as.character(df$donorChromosome), 4)
#   
#   genesD <- lapply(as.list(unique(df$chromosomenum)), function(ch){
#     if( !(ch %in% names(mapchr)) ){
#       return(NULL)
#     }
#     thisMap <- mapchr[[ch]]
#     thisDf <- df[ df$chromosomenum == ch, ]
#     t(sapply(as.list(1:nrow(thisDf)), function(x){
#       getGeneID(thisDf$donorEnd[x], thisMap)
#     }))
#   })
#   names(genesD) <- unique(df$chromosomenum)
#   huh <- do.call(rbind, genesD)
#   names(df)
#   
  
  genes_d <- apply(df, 1, function(x){
      chromosomenum <- substring(as.character(x['donorChromosome']), 4)
    if(is.null(mapchr[[chromosomenum]])){return(c("NA", "NA"))}
    return(getGeneID(as.numeric(x['donorEnd']), mapchr[[chromosomenum]]))
  }) 
  genes_d <- t(genes_d)
  colnames(genes_d) <- c("g_ids", "g_l")
  
  #gd <- do.call('rbind', genes_d)
  # Get IDs of acceptor genes
  genes_a <- apply(df, 1, function(x){
    chromosomenum <- substring(as.character(x['acceptorChromosome']), 4)
    if(is.null(mapchr[[chromosomenum]])){return(c("NA", "NA"))}
    return(getGeneID(as.numeric(x['acceptorEnd']), mapchr[[chromosomenum]]))
  })
  genes_a <- t(genes_a)
  colnames(genes_a) <- c("g_ids", "g_l")
  
  delta <- length(outputnames) - ncol(df) - ncol(genes_d) - ncol(genes_a)
  padding <- data.frame(matrix(NA, nrow = nrow(df), ncol = delta))
  output <- cbind(df[,1:5], genes_d[,'g_ids'], genes_a[,'g_ids'], genes_d[,'g_l'], genes_a[,'g_l'], df[,6:ncol(df)], padding)
  colnames(output) <- outputnames
  filename=paste("/home/onikolov/projects/gene_fusion_unc/data_gene_length/", ctnames[ct], ".tab", sep ="")
  write.table(output, file=filename, col.names=TRUE, row.names=FALSE, sep="\t")
  
} # end for ct
