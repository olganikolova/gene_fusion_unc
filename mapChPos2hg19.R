# =======================================
# Author: Olga Nikolova
# E-mail: olga.nikolova@gmail.com
# =======================================
# Map chromosome positions to genes
# Map junction coordinates to exons
# =======================================
library(biomaRt)

# Create map
ensembl=useMart("ensembl",dataset="hsapiens_gene_ensembl")
map <- getBM(c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position"), mart = ensembl)

# List of submaps per chrom
mapchr <- lapply(as.list(c(1:22,"X","Y")),function(x){return(map[map$chromosome_name == x,])})
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
  
  g <- map[intersect(which(map$start_position <= pos), 
                   which(map$end_position >= pos)),]

 if( nrow(g) > 0){
    g_e <- toString(g$ensembl_gene_id, sep=",")
    g_h <- toString(g$hgnc_symbol, sep=",")
    g_ids <- paste(g_e, g_h, sep=";")
    g_l <- toString(g$end_position - g$start_position, sep=",") 

    return(cbind(g_ids, g_l))
  }   else {
    return("NA")
  } 
}

outputnames <- c("sampleID","cancerType","sequencingID","sequencingLane","sequencingBarcode","donorGeneID",
                 "acceptorGeneID","donorGeneLength","acceptorGeneLength","donorChromosome","acceptorChromosome","donorEnd","acceptorEnd","name",
                 "coverage","strand","itemRgb","blockCount","blockSizes","blockStarts","entropy","flank_string_case",
                 "flank_string","min_mismatch","max_mismatch","average_mismatch","maximal_of_minimal_acceptor_site_length",
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

clist <- as.list(c("BLCA","CESC","BRCA","HNSC","KIRP","KICH","KIRC","LGG","LIHC","LUAD","LUSC","PAAD","PRAD","SKCM","THCA"))
ctnames <- unlist(clist)
# input directory
path <- "/Shares/work/DAT_109__TCGA_unc_rnaseq_genefusion/Data/mapped2TCGA_final/"

d.clist <- lapply(clist, function(type){
  name <- paste(path, type,".txt", sep="")
  return(read.csv(name, header=T, sep="\t"))
})

for(ct in 1:length(d.clist)){
  df <- d.clist[[ct]]
  
  # Get IDs of donor genes
  genes_d_l <- apply(df, 1, function(x){
    chromosomenum <- substring(as.character(x['donorChromosome']), 4)
    if(is.null(mapchr[[chromosomenum]])){return("NA")}
    return(getGeneID(x['donorEnd'], mapchr[[chromosomenum]]))
  })  
  genes_d <- do.call("rbind", genes_d_l)
  
  # Get IDs of acceptor genes
  genes_a_l <- apply(df, 1, function(x){
    chromosomenum <- substring(x['acceptorChromosome'], 4)
    if(is.null(mapchr[[chromosomenum]])){return("NA")}
    return(getGeneID(x['acceptorEnd'], mapchr[[chromosomenum]]))
  })
  genes_a <- do.call("rbind", genes_a_l)
  
  delta <- length(outputnames) - ncol(df) - 4
  padding <- data.frame(matrix(NA, nrow = nrow(df), ncol = delta))
  output <- cbind(df[,1:5], genes_d[,'g_ids'], genes_a[,'g_ids'], genes_d[,'g_l'], genes_a[,'g_l'], df[,6:ncol(df)], padding)
  colnames(output) <- outputnames
  filename=paste("/home/onikolov/projects/gene_fusion_unc/data_gene_length/", ctnames[ct], ".tab", sep ="")
  write.table(output, file=filename, col.names=TRUE, row.names=FALSE, sep="\t")
  
} # end for ct
