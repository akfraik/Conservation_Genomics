# Set working directory
setwd("/Users/alexandra.fraik")
# Load requried packages
library(vcfR)
library(ggplot2)
library(dplyr)
library(tidyr)
library("adegenet")
library("adespatial")
library(gridExtra)
library("grid")
library("hierfstat")
library("mmod")
library(graph4lg)

# Read in the genepop file format from Kristy
file2<-read.csv("/Users/alexandra.fraik/Documents/Fisher_Stuff/Cabinets_Genetics_2023/Input_Files/K6_Pops.csv")
colnames(file2)[1]<-"ID"
file3<-file2[c("ID","Haplotype","Lat","Long","pop","Year_Bracket","colors","points")]
dim(file3)

## OK- so this is definitely what is going on, now how do we fix this...
df<-read.csv("/Users/alexandra.fraik/Documents/Fisher_Stuff/Cabinets_Genetics_2023/Input_Files/Other_other_df.csv")
test<-cbind(df[ 1 ],mapply(paste0, df[, seq(2, 31, 2)], ":",df[, seq(3, 31, 2)]))
file4<-merge(test,file3,by="ID")
dim(file4)
genos<-file4[,c(2:16)]
genos[genos=="0:0"]<-NA
inds<-file4[,c(1)]
pop<-file4[,c(20)]
init<-df2genind(genos,sep=":",ind.names=inds,pop=pop,ploidy=2,NA.char="0:0",type="codom",strata=NULL,hierarchy=NULL)
strata(init)<-file4
test_df<-setPop(init,~pop)
table(test_df@pop,test_df@strata$pop)
length(match(rownames(init@tab),file4$ID))
## Make sure this number is 397 before proceeding

## For now, let's just go ahead and do this this way
test_df<-setPop(test_df,~ID)
Ds_pop<-mat_gen_dist(test_df,dist="PCA")
pcoa(as.matrix(Ds_pop))
Ds_mod<-as.dist(Ds_pop)

## Let's try to make a PCA a different way
PCA<-dudi.pco(Ds_mod,scannf=FALSE,nf=4) 
PCA$eig
scatter(PCA,xax=1,yax=2,clab.row=1,posieig="topright")
tab<-PCA$tab
tab<-tab[c(1:5)]
tab$ID<-rownames(tab)
df<-merge(tab,file4,by="ID")

## We'll come back to merging metadata and such, but let's just take a quick look at this for a moment
dapc_df<-data.frame(sample.id=df$ID,Haplotype=df$Haplotype,Lat=df$Lat,Long=df$Long,pop=df$pop,
                    Year_Bracket=df$Year_Bracket,A1=df$A1,A2=df$A2,A3=df$A3,A4=df$A4,colors=df$colors,
                    points=df$points,stringsAsFactors=FALSE)

## Let's once again merge the metadata file with our microsatellite genotype dadfase
neworder<-c("Pre 2000","2000-2004","2005-2009","2010-2014","2015-2022")
pca_df<-arrange(transform(dapc_df,Year_Bracket=factor(Year_Bracket,levels=neworder),Year_Bracket))
pca_df$Year_Bracket<-as.factor(pca_df$Year_Bracket)

## Let's get the colors and pop order set up
pop_order<-c("Chetwyn","burnslake","princegeorge","anahimlake","WilliamsLk",
             "Adamslake","Midwest","Selkirks","Cabinets","Seeley_Swan","StJoe","N_Lolo","Lolo","Clearwater",
             "EF_Granite_Creek","Bitterroot","Red_River","Salmon_Challis")
pca_df<-arrange(transform(pca_df,pop=factor(pop,levels=pop_order),pop))
colors<-as.character(pca_df$colors)
names(colors)<-as.character(pca_df$pop)

## Finally, let's get the haplotypes order and points set up
hap_order<-c("1b","4","6","7","11","9","12","1","5","10")
df<-arrange(transform(pca_df,Haplotype=factor(Haplotype,levels=hap_order),Haplotype))
points<-as.numeric(df$points)
names(points)<-as.character(df$Haplotype)

## Location Only with colors
#pdf(file="PCA/Map_Lat_Long.pdf",width=8,height=8)
plot<-plot+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-ggplot(df,aes(x=Long,y=Lat),group=pop)
plot<-plot+geom_jitter(aes(fill=pop,color=pop,shape=Haplotype),size=5,width=0.3,height=0.3)
plot<-plot+scale_color_manual(name="Populations",values=colors)
plot<-plot+scale_fill_manual(name="Populations",values=colors)
plot<-plot+scale_shape_manual(name="Haplotypes",values=points)
plot<-plot+theme_bw()
plot<-plot+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plot<-plot+guides(fill="none",color=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)),shape=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)),size="none")
plot<-plot+theme(legend.position="right")
plot
#dev.off()

## Location Only with colors
#pdf(file="PCA/Map_Lat_Long_Less.pdf",width=8,height=8)
plot<-plot+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-ggplot(df,aes(x=Long,y=Lat),group=pop)
plot<-plot+geom_jitter(aes(fill=pop,color=pop),size=5,width=0.3,height=0.3)
plot<-plot+scale_color_manual(name="Populations",values=colors)
plot<-plot+scale_fill_manual(name="Populations",values=colors)
plot<-plot+theme_bw()
plot<-plot+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plot<-plot+guides(fill="none",color=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)),shape=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)),size="none")
plot<-plot+theme(legend.position="right")
plot
#dev.off()

## Plot these groups 1-3: Axis 2 changes depending on if you use FST or DPS
#pdf(file="PCA/PCA_12_sans_hts.pdf",width=8,height=6)
plot<-plot+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-ggplot(df,aes(x=A1,y=A2),group=pop)
plot<-plot+geom_jitter(aes(fill=pop),color="black",size=5,width=0.3,height=0.3,shape=21)
plot<-plot+scale_fill_manual(name="Populations",values=colors)
plot<-plot+labs(x="PCA Axis 1",y="PCA Axis 2")
plot<-plot+theme_bw()
plot<-plot+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plot<-plot+guides(fill=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)),shape=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)),size="none")
plot<-plot+theme(legend.position="right")
plot
#dev.off()

## Plot these groups 1-3: Axis 2 changes depending on if you use FST or DPS
#pdf(file="PCA/PCA_12.pdf",width=8,height=6)
plot<-ggplot(df,aes(x=A1,y=A2),group=pop)
plot<-plot+geom_jitter(aes(fill=pop,color=pop,shape=Haplotype),size=5)
plot<-plot+scale_color_manual(name="Populations",values=colors)
plot<-plot+scale_fill_manual(name="Populations",values=colors)
plot<-plot+scale_shape_manual(name="Haplotypes",values=points)
plot<-plot+labs(x="PCA Axis 1",y="PCA Axis 2")
plot<-plot+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plot<-plot+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-plot+theme(legend.position="right")
plot<-plot+guides(color="none",fill="none",size="none",shape=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)))
plot
#dev.off()

## Plot these groups 1-3: Axis 2 changes depending on if you use FST or DPS
pdf(file="PCA/PCA_34_sans_hts.pdf",width=8,height=6)
plot<-plot+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-ggplot(df,aes(x=A3,y=A4),group=pop)
plot<-plot+geom_jitter(aes(fill=pop),color="black",size=5,width=0.3,height=0.3,shape=21)
plot<-plot+scale_fill_manual(name="Populations",values=colors)
plot<-plot+labs(x="PCA Axis 3",y="PCA Axis 4")
plot<-plot+theme_bw()
plot<-plot+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plot<-plot+guides(fill=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)),shape=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)),size="none")
plot<-plot+theme(legend.position="right")
plot
dev.off()

## Plot these groups 1-3: Axis 2 changes depending on if you use FST or DPS
pdf(file="PCA/PCA_34.pdf",width=8,height=6)
plot<-ggplot(df,aes(x=A3,y=A4),group=pop)
plot<-plot+geom_jitter(aes(fill=pop,color=pop,shape=Haplotype),size=5)
plot<-plot+scale_color_manual(name="Populations",values=colors)
plot<-plot+scale_fill_manual(name="Populations",values=colors)
plot<-plot+scale_shape_manual(name="Haplotypes",values=points)
plot<-plot+labs(x="PCA Axis 3",y="PCA Axis 4")
plot<-plot+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plot<-plot+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-plot+theme(legend.position="right")
plot<-plot+guides(color="none",fill="none",size="none",shape=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)))
plot
dev.off()

## Analyse how much percent of genetic variance is explained by each axis
pca1<-PCA$eig[1]
PCA_eig<-PCA$eig
percent<-PCA_eig/sum(PCA$eig)*100
percent

## Let's save the PCA plot of the variance explained
pdf(file="PCA/PCA_PCA_Variance_Explained.pdf",width=10,height=8)
plot<-barplot(percent,ylab="Genetic variance explained by eigenvectors (%)",ylim=c(0,12),names.arg=round(percent,1))
plot
dev.off()
