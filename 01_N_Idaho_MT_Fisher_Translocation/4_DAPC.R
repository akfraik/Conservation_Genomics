# Set working directory
setwd("/Users/alexandra.fraik")
# Load R packages
library(vcfR)
library(ggplot2)
library(dplyr)
library(stringr)
library("adegenet")
library("adespatial")
library(gridExtra)
library("grid")
library("hierfstat")
library("mmod")

## Read in the metadata and subset the columns to what you are looking for
file2<-read.csv("Metadata.csv")
file2<-file2[c("ID","Haplotype","Lat","Long","pop","Year_Bracket","points","colors")]
# Read in the microsat dataset
df<-read.csv("Microsat_df_fishers.csv")
test<-cbind(df[ 1 ],mapply(paste0, df[, seq(2, 31, 2)], ":",df[, seq(3, 31, 2)]))
file4<-merge(test,file2,by="ID")
dim(file4)
genos<-file4[,c(2:16)]
genos[genos=="0:0"]<-NA
inds<-file4[,c(1)]
pop<-file4[,c(20)]
init<-df2genind(genos,sep=":",ind.names=inds,pop=pop,ploidy=2,NA.char="NA",type="codom",strata=NULL,hierarchy=NULL)
strata(init)<-file4
test_df<-setPop(init,~pop)
table(test_df@pop,test_df@strata$pop)
length(match(rownames(init@tab),file4$ID))
## Make sure this number is 397 before proceeding

## PCAs and DAPCs don't take NAs - that's why the function below uses the mean to impute
testing<-test_df$tab
f1 <- function(vec) {
  m <- mean(vec, na.rm = TRUE)
  vec[is.na(vec)] <- m
  return(vec)}
test<-apply(testing,2,f1)
genind_clusters<-find.clusters(test,max.n.clust=20,n.pca=80,n.clust=4)
genind_clusters
dapc<-dapc(test,genind_clusters$grp,n.pca=40,n.da=2)

## Make a plot of the highly contributive microsatellites
pdf(file="DAPC/K5/DAPC_Loadings.pdf",width=10,height=8)
contrib <- loadingplot(dapc$var.contr)
contrib
dev.off()

## Extract necessary information from DAPC
pca.here<-as.data.frame(dapc$pca.loadings)
pca.here<-pca.here[1:4]
indiv_loadings<-as.data.frame(dapc$loadings)
tab<-dapc$tab
tab_df<-tab[1:4]
group<-dapc$grp
coords<-dapc$ind.coord
posterior<-dapc$posterior
dapc_df<-data.frame(tab_df,group,coords,posterior)
dapc_df$ID<-rownames(dapc_df)
df<-merge(dapc_df,file4,by="ID")

## Analyse how much percent of genetic variance is explained by each axis
pca1<-dapc$pca.eig[1]
pca<-dapc$pca.eig
percent1<-pca1/sum(dapc$pca.eig)*100
percent<-pca/sum(dapc$pca.eig)*100
percent

## Let's save the DAPC plot of the variance explained
pdf(file="DAPC/K5/DAPC_Variance_Explained.pdf",width=10,height=8)
plot<-barplot(percent,ylab="Genetic variance explained by eigenvectors (%)",ylim=c(0,12),names.arg=round(percent,1))
plot
dev.off()

## We'll come back to merging metadata and such, but let's just take a quick look at this for a moment
dapc_df<-data.frame(sample.id=df$ID,Haplotype=df$Haplotype,Lat=df$Lat,Long=df$Long,pop=df$pop,
                    Year_Bracket=df$Year_Bracket,LD1=df$LD1,LD2=df$LD2,colors=df$colors,
                    points=df$points,group=df$group,stringsAsFactors=FALSE)
write.csv(dapc_df,"DAPC/K5/DAPC_K5.csv",row.names=FALSE)

## Let's once again merge the metadata file with our microsatellite genotype dadfase
neworder<-c("Pre 2000","2000-2004","2005-2009","2010-2014","2015-2022")
pca_df<-arrange(transform(dapc_df,Year_Bracket=factor(Year_Bracket,levels=neworder),Year_Bracket))
pca_df$Year_Bracket<-as.factor(pca_df$Year_Bracket)

## Let's get the colors and pop order set up
pop_order<-c("UNKBC","Chetwyn","burnslake","princegeorge","anahimlake","WilliamsLk",
             "Adamslake","Midwest","Selkirks","Cabinets","Seeley_Swan","StJoe","N_Lolo","Lolo","Clearwater",
             "EF_Granite_Creek","Bitterroot","Red_River","Salmon_Challis")
pca_df<-arrange(transform(pca_df,pop=factor(pop,levels=pop_order),pop))
colors<-as.character(pca_df$colors)
names(colors)<-as.character(pca_df$pop)

## Finally, let's get the haplotypes order and points set up
hap_order<-c("1b","4","6","7","11","9","12","1","5","10")
dapc_df<-arrange(transform(pca_df,Haplotype=factor(Haplotype,levels=hap_order),Haplotype))
points<-as.numeric(dapc_df$points)
names(points)<-as.character(dapc_df$Haplotype)
dapc_df$group<-as.factor(dapc_df$group)

################################################################################################################################################      
############################### Generate DAPC that shows group membership ##############################################################
################################################################################################################################################
## DAPC Axis 1 VS DAPC Axis 2
pdf(file="DAPC/K5/DAPC_Groups_K5_Haplotype.pdf",width=10,height=8)
plot<-ggplot(dapc_df,aes(x=group,y=Haplotype),group=group)
plot<-plot+geom_jitter(aes(color=pop,fill=pop,shape=group),size=4)
plot<-plot+labs(y="Haplotype",x="Population")
plot<-plot+scale_color_manual(name="Sampling Location",values=colors)
plot<-plot+scale_fill_manual(name="Sampling Location",values=colors)
plot<-plot+scale_shape_manual(name="DAPC Groupings",values=c(21,22,23,24,25))
plot<-plot+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plot<-plot+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-plot+theme(legend.position="right")
plot<-plot+guides(fill="none",color=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=3)),shape=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=3)),size="none")
plot
dev.off()

## DAPC Axis 1 VS DAPC Axis 2
pdf(file="DAPC/K5/DAPC_Groups_K5_Latitude.pdf",width=10,height=8)
plot<-ggplot(dapc_df,aes(x=group,y=Lat),group=group)
plot<-plot+geom_jitter(aes(color=pop,fill=pop,shape=group),size=4)
plot<-plot+labs(y="Latitude",x="Population")
plot<-plot+scale_color_manual(name="Sampling Location",values=colors)
plot<-plot+scale_fill_manual(name="Sampling Location",values=colors)
plot<-plot+scale_shape_manual(name="DAPC Groupings",values=c(21,22,23,24,25))
plot<-plot+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plot<-plot+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-plot+theme(legend.position="right")
plot<-plot+guides(fill="none",color=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=3)),shape=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=3)),size="none")
plot
dev.off()

################################################################################################################################################
######## DAPC Axis 1 VS DAPC Axis 2 ########
################################################################################################################################################
pdf(file="DAPC/K5/DAPC_K5.pdf",width=12,height=7)
plot<-ggplot(dapc_df,aes(x=LD1,y=LD2),group=group)
plot<-plot+geom_point(aes(color=pop,fill=pop,shape=group),size=5)
plot<-plot+labs(x="DAPC Axis 1",y="DAPC Axis 2")
plot<-plot+scale_color_manual(name="Populations",values=colors)
plot<-plot+scale_fill_manual(name="Populations",values=colors)
plot<-plot+scale_shape_manual(name="DAPC Groupings",values=c(21,22,23,24,25))
plot<-plot+theme(axis.text.y=element_text(color="black",size=10),axis.title.y=element_text(size=12,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=10),axis.title.x=element_text(size=12,face="bold"))
plot<-plot+theme(strip.text=element_text(face="bold",size=12),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-plot+guides(fill="none",color=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)),shape=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)),size="none")
plot<-plot+theme(legend.position="right")
plot
dev.off()

## Location Only with colors
pdf(file="DAPC/K5/Map_K5_Group.pdf",width=8,height=6)
plot<-ggplot(dapc_df,aes(x=Long,y=Lat,shape=group))
plot<-plot+geom_jitter(aes(color=pop,fill=pop,shape=group),size=4,width=0.3,height=0.3)
plot<-plot+scale_color_manual(name="Haplotypes",values=colors)
plot<-plot+scale_fill_manual(name="Haplotypes",values=colors)
plot<-plot+scale_shape_manual(name="Populations",values=c(21,22,23,24,25))
plot<-plot+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plot<-plot+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-plot+guides(fill="none",color=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)),shape=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)),size="none")
plot<-plot+theme(legend.position="right")
plot
dev.off()

## Haplotype VS Time- location Only with colors
pdf(file="DAPC/K5/Map_K5_Haplotype.pdf",width=10,height=8)
plot<-ggplot(dapc_df,aes(x=Long,y=Lat,shape=Haplotype),group=group)
plot<-plot+geom_jitter(aes(color=pop,fill=pop,shape=Haplotype),size=3,width=0.3,height=0.3)
plot<-plot+scale_color_manual(name="Haplotypes",values=colors)
plot<-plot+scale_fill_manual(name="Haplotypes",values=colors)
plot<-plot+scale_shape_manual(name="Haplotypes",values=points)
plot<-plot+facet_wrap(group,switch="x")
plot<-plot+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plot<-plot+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-plot+guides(fill="none",color=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)),shape=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)),size="none")
plot<-plot+theme(legend.position="right")
plot
dev.off()

## Haplotype VS Time
pdf(file="DAPC/K5/DAPC_K5_Hap_Time_Diff.pdf",width=10,height=8)
plot<-ggplot(dapc_df,aes(x=LD1,y=Haplotype),group=group)
plot<-plot+geom_jitter(aes(color=pop,fill=pop,shape=Haplotype),size=4,height=0.3)
plot<-plot+labs(x="DAPC Axis 1",y="haplotypes")
plot<-plot+scale_color_manual(name="Populations",values=colors)
plot<-plot+scale_fill_manual(name="Populations",values=colors)
plot<-plot+scale_shape_manual(name="Haplotypes",values=points)
plot<-plot+facet_grid(.~group,switch="x")
plot<-plot+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plot<-plot+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-plot+guides(fill="none",color=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)),shape=guide_legend(ncol=2,byrow=FALSE,override.aes=list(size=5)),size="none")
plot<-plot+theme(legend.position="right")
plot<-plot+guides(fill="none",color=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=3)),shape=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=3)),size="none")
plot
dev.off()
