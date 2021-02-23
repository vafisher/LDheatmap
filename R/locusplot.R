#library(LDheatmap)
#library(RColorBrewer)
#source("https://bioconductor.org/biocLite.R")
#biocLite("rtracklayer")
#library(rtracklayer)

#loci<-read.table("bmi_plot/sig.locs.txt", as.is=T)
#                 , as.is=T)$V1
#loci

#source("mod4.addScatterplot.R")
#environment(mod.LDhtm)<-environment(LDheatmap)
#environment(add2Scatterplot)<-environment(LDheatmap)


ldh4<-add2Scatterplot(ldh1,list(GWA.p.y,Wald.p.y), spc=c(1,4), ylab=c("GWAS","AnnoRE"), type="both")

ldh5<-editGrob(ldh4$LDheatmapGrob,gPath("geneMap","title"),label=NULL)
LDheatmap(ldh5)
grid.newpage()
grid.draw(ldh5)

loc<-loci[2,]


locplot<-function(loc){
  snp<-loc[1]
  gene<-loc[2]
  loc.df<-read.table(paste("data/",snp,".out.txt",sep=""), header=T, as.is=T)
  names(loc.df)
  dim(loc.df)
  head(loc.df)
  top<-loc.df[which(loc.df$Wald.p==min(loc.df$Wald.p,na.rm=T)),]
  top
  loc.df.order<-loc.df[order(loc.df$BP),]
  GWA.p.y<--log10(loc.df.order$p)
  Wald.p.y<--log10(loc.df.order$Wald.p)

  ld.long<-read.table(paste("data/",snp,".1000G.loc.hap.ld", sep=""), header=T, as.is=T)
  ld.loc<-ld.long[which(ld.long$POS1 %in% loc.df$BP & ld.long$POS2 %in% loc.df$BP),]

  LDmat<-reshape(ld.loc, direction="wide",v.names=c("R.2"),idvar=c("POS1"),timevar=c("POS2"),drop=c("CHR","D","Dprime"))
  LDmat<-rbind(cbind(NA,LDmat[,-(1:2)]),NA)
  LDmat[lower.tri(LDmat)]<-t(LDmat)[lower.tri(LDmat)]
  diag(LDmat)<-1
  LDmat<-as.matrix(LDmat)

  ldh1<-LDheatmap(LDmat, genetic.distances=loc.df.order$BP, flip=T, add.key=F, add.map=T, title=paste(gene,"locus: chr",top$CHR), color=heat.colors(20))

  tiff(paste("bmi_plot/",gene,".tiff",sep=""),width=85,height=85,units="mm",res=300)
  ldh4<-add2Scatterplot(ldh1,list(GWA.p.y,Wald.p.y), spc=c(1,4), ylab=c("GWAS","AnnoRE"), type="both")

  dev.off()
}



#sapply(loci[1:6],locplot)

#par(mfrow=c(3,2))




