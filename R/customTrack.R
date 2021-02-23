#' @name customTrack
#' @aliases customTrack
#' @title Produce recombination rate plot.
#' @description Plot average rates of recombination from the deCODE genetic map for a specified genetic sequence.
#' @usage customTrack(minRange, maxRange, chromosome, genome = "hg19", vp = viewport(x = 0,
#'y = 0.99, height = 0.04, just = c("left", "top")), view = "dense")
#' @param minRange The sequence minimum range in base pairs.
#' @param maxRange The sequence maximum range in base pairs.
#' @param chromosome A character string identifying the chromosome.
#' @param genome The genome assembly to use. The default is hg19, the most recent human genome assembly on
#'the UCSC genome browser.
#' @param vp A \code{viewport}.
#' @param view Display mode. Possible values are \code{"dense"} (the default), \code{"squish"},
#'\code{"pack"} and \code{"full"}.
#' @return A \code{grob} representing custom tracks.
#' 
#' @author Virginia Fisher <vafisher@bu.edu> and more
#' @examples \dontrun{
#'grid.newpage()
#'recombRate(139000000, 140000000, "chr7", "hg18")
#'grid.newpage()
#'pushViewport(viewport(width=0.8, x=0.2, just="left"))
#'recombRate(129000000, 140000000, "chr7", "hg18", view="full")
#'popViewport()
#'}
#' @keywords hplot
#' @export

# ldheatmap - Plots measures of pairwise linkage disequilibria for SNPs
# Copyright (C) 2004  J.Shin, S. Blay, N. Lewin-Koh, J.Graham, B.McNeney

# To cite LDheatmap in publications use:
# Shin J-H, Blay S, McNeney B and Graham J (2006). LDheatmap: An R
# Function for Graphical Display of Pairwise Linkage Disequilibria
# Between Single Nucleotide Polymorphisms. J Stat Soft, 16 Code Snippet 3

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

###########################################################################



#________________Plot custom track from file_______________##
customTrack <- function(minRange, maxRange, chromosome, ctitle="Custom Track", vp=viewport(x=0, y=0.99, height=0.04, just=c("left", "top")), view="dense") {
  requireNamespace("grid")
  map_len <- convertX(vp$width, "npc", valueOnly=TRUE)
  Range <- maxRange - minRange
  vp$xscale <- c(minRange, maxRange)
  vp$name <- "customTrack"

  session <- rtracklayer::browserSession()
  GenomeInfoDb::genome(session) <- genome
  query<-rtracklayer::ucscTableQuery(session, "recombRate", rtracklayer::GRangesForUCSCGenome(genome, chromosome,IRanges::IRanges(minRange, maxRange)))
  t<-rtracklayer::getTable(query)
  if (!dim(t)[1]) {	# we got an empty table
    print ("The genetic region of the data does not correspond to any recombination rate data in the UCSC genome browser")
    return()
  }

  customTrack_title <- textGrob(ctitle,
	gp=gpar(fontsize=7, fontfamily="mono"), y=1, just=c("centre", "center"),
 				name="customTrack_title", default.units="native")

  vp_height <- convertX(vp$height, "npc", valueOnly=TRUE)
  rect_height <- 0.5
     
  customTrack <- gTree(children=gList(customTrack_title), name="customTrack", vp=vp)


  # Read recomb rate for every every row in table t
  for (i in 1:dim(t)[1]) {

       y <- 0  
       
       value <- as.integer(t[i,"decodeAvg"]*200)
       color <- as.character(cut(value, breaks=c(0,seq(length=9, from=320, to=600), 10000),
		labels=grey.colors(10)[10:1]))
       
       rect <- rectGrob(x=max(minRange, t[i, "chromStart"]), default.units="native",
       width=min(maxRange, t[i, "chromEnd"]) - max(minRange, t[i, "chromStart"]),
       y=y, height=rect_height, just=c("left", "bottom"), gp=gpar(col="transparent", fill=color),
	name=paste("rect",i, sep=""))
       customTrack <- addGrob(customTrack , rect)
       if (view == "full" || view == "pack") {
          x <- 1;        if (view == "pack") { x <- i }
          rectText <- textGrob(x=max(minRange, t[x, "chromStart"]), y=y, 
		paste(round(t[i,"decodeAvg"],1)," cM/Mb (Avg) ",sep=""), just=c("right","top"), 
        	gp=gpar(fontsize=7, fontfamily="mono"), default.units="native",
		name=paste("text",i, sep=""))
          customTrack <- addGrob(customTrack , rectText)
       }
  }
  grid.draw(customTrack)
  return(customTrack)
}

