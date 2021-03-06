#' @name GIMAP5
#' @aliases GIMAP5
#' @docType data
#' @title Example data set for LDHeatmap
#' @description SNP genotypes on HapMap founders for SNPs spanning the GIMAP5 gene.
#' @usage data(GIMAP5)
#' @format GIMAP5 is a list with three elements: snp.data, snp.support and
#'subject.support. snp.data is a \code{SnpMatrix}
#'object containing the SNP genotypes. Rows correspond to
#'subjects and columns correspond to  SNPs.
#'snp.support is a data frame with the following columns:
#'\tabular{rlll}{
#'[,1] \tab dbSNPalleles   \tab character \tab alleles at each SNP\cr
#'[,2] \tab Assignment \tab character \tab same as dbSNPalleles\cr
#'[,3] \tab Chromosome   \tab character \tab chromosome (chr7 for all)\cr
#'[,4] \tab Position    \tab numeric \tab physical position\cr
#'[,5] \tab Strand   \tab character \tab strand (all "+")\cr
#'}
#'subject.support is a one-column data frame with:
#'\tabular{rlll}{
#'  [,1] \tab pop   \tab character \tab HapMap population of each subject \cr
#'}
#'
#' @details SNP genotypes from HapMap release 27
#'for SNPs in a 10KB region spanning
#'the GIMAP5 gene. Data are on founders from each of the 11 HapMap
#'phase III populations:
#'  \tabular{ll}{
#'    ASW \tab African ancestry in Southwest USA \cr
#'    CEU \tab Utah residents with Northern and Western European ancestry from the CEPH collection \cr
#'    CHB \tab Han Chinese in Beijing, China \cr
#'    CHD \tab Chinese in Metropolitan Denver, Colorado \cr
#'    GIH \tab Gujarati Indians in Houston, Texas \cr
#'    JPT \tab Japanese in Tokyo, Japan \cr
#'    LWK \tab Luhya in Webuye, Kenya \cr
#'    MEX \tab Mexican ancestry in Los Angeles, California \cr
#'    MKK \tab Maasai in Kinyawa, Kenya \cr
#'    TSI \tab Toscani in Italia \cr
#'    YRI \tab Yoruba in Ibadan, Nigeria \cr
#'  }
#'Only those SNPs with minor allele frequency greater
#'than 5\% in all populations were retained.
#'The base positions are from NCBI build 36
#'(UCSC genome hg18).
#'
#' @seealso \code{\link{GIMAP5.CEU}}
#' @examples data(GIMAP5) 
#'#Now do a lattice plot with LDheatmaps in the panels
#'library(lattice)
#'pop<-GIMAP5$subject.support$pop
#'n<-nrow(GIMAP5$snp.data)
#'xyplot(1:n ~ 1:n | pop, type="n", scales=list(draw=FALSE), xlab="", ylab="",
#'       panel=function(x, y, subscripts,...) {
#'         LDheatmap(GIMAP5$snp.data[subscripts,],GIMAP5$snp.support$Position, 
#'                   newpage=FALSE)})
#'rm(pop,n)
#' @source International HapMap Project \url{www.hapmap.org}
#' @references The International HapMap Consortium. A haplotype map of
#'the human genome. Nature 437, 1299-1320. 2005.
#' @keywords datasets
NULL
