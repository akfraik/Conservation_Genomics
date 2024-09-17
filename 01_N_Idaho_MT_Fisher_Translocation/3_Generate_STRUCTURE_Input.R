# Set working directory
setwd("/Users/alexandra.fraik")
# Load R Packages
library(ggplot2)
library(dplyr)
library("adegenet")
library("adespatial")
library(gridExtra)
library("grid")
library("hierfstat")
library("mmod")
library("poppr")
library(reshape2)
library(tidyr)
library(graph4lg)

## OK- so this is definitely what is going on, now how do we fix this...
df<-read.csv("/Users/alexandra.fraik/Documents/Fisher_Stuff/Cabinets_Genetics_2023/Input_Files/Other_other_df.csv")
test<-cbind(df[ 1 ],mapply(paste0, df[, seq(2, 31, 2)], ":",df[, seq(3, 31, 2)]))
genos<-test[,c(2:16)]
genos[genos=="0:0"]<-NA
inds<-file5[,c(1)]
pop<-file5[,c(18)]
init<-df2genind(genos,sep=":",ind.names=inds,pop=pop,ploidy=2,NA.char="NA",type="codom",strata=NULL,hierarchy=NULL)
strata(init)<-file5
test_df<-setPop(init,~group)
init<-test_df[!test_df@pop %in% c("4","5"),]
file5<-file5[c("ID","Lat","Long","Haplotype","pop","Year_Bracket","colors","points")]
dim(file5)
init

## Nope, this doesn't seem to work
genind_to_genepop(init,output="STRUCTURE/All.txt")
## Then, use this file in SPIDER.. seems to be the only way to get it to work

##########################################################################################
## Then, let's plot the output after running the GUI/running it on the command line
# I ran it on the command line, but either works
##########################################################################################
file2<-read.csv("Metadata.csv")
colnames(file2)[1]<-"ID"
file2<-file2[c("ID","pop","colors","points","Year_Bracket","Haplotype","Lat","Long")]
file2$ID2<-paste(file2$pop,file2$ID,sep="_")
  
##########################################################################################
## Read in STRUCTURE files
##########################################################################################
K2<-read.csv("STRUCTURE/Admix/K2.indfile.csv",header=FALSE)
colnames(K2)[1]<-"Order"
nosex<-read.csv("STRUCTURE/Admix/All_Plus_BC.txt",header=FALSE)
colnames(nosex)[1]<-"Order"
colnames(nosex)[2]<-"ID2"
here<-merge(K2,nosex,by="Order")
Pre<-merge(here,file2,by="ID2")
Pre<-Pre[!colnames(Pre) %in% "ID2"]
Pre<-Pre[!names(Pre) %in% "Order"]
here<- Pre %>% 
  gather(Variable,Ancestry,-ID,-pop,-colors,-points,-Year_Bracket,-Haplotype,-Lat,-Long)
here <- here %>% 
  drop_na(Variable)

## Now let's do pop order and colors
pop_order<-c("Chetwyn","burnslake","princegeorge","anahimlake","WilliamsLk",
             "Adamslake","Midwest","Selkirks","Cabinets","Seeley_Swan","StJoe","N_Lolo","Lolo","Clearwater",
             "EF_Granite_Creek","Bitterroot","Red_River","Salmon_Challis")
df<-arrange(transform(here,pop=factor(pop,levels=pop_order),pop))
colors<-c("violet","darkblue")
colors<-as.character(df$colors)
df$Variable<-as.factor(df$Variable)
names(colors)<-as.character(df$Variable)

## Finally, let's get the haplotypes order and points set up
hap_order<-c("1b","4","6","7","11","9","12","1","5","10")
man<-arrange(transform(df,Haplotype=factor(Haplotype,levels=hap_order),Haplotype))
points<-as.numeric(man$points)
names(points)<-as.character(man$Haplotype)
man$Genotype<-as.factor(man$Ancestry)
K2<-man

## K = 2
pdf("STRUCTURE/Admix/K2.pdf",height=15,width=55)
plota<-ggplot(K2,aes(x=ID,y=Ancestry),group=Variable)+theme(legend.position="none")
plota<-plota+geom_bar(aes(y=Ancestry,x=ID,fill=Variable),stat="identity")
plota<-plota+theme(axis.text.y=element_text(color="black",size=25),axis.title.y=element_text(size=40,face="bold"))
plota<-plota+theme(axis.text.x=element_blank(),axis.title.x=element_text(size=40,face="bold"))
plota<-plota+ylab(label="Ancestry")+xlab(label="Sample Locations")
plota<-plota+scale_fill_manual(values=c("darkblue","violet"))
plota<-plota+scale_color_manual(values=c("darkblue","violet"))
plota<-plota+facet_grid(~pop,switch="x",scales="free_x",space="free_x")
plota<-plota+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plota<-plota+theme(legend.text=element_text(size=rel(1.5)))
plota<-plota+theme(legend.title=element_text(face="bold",size=rel(1.5)))
plota<-plota+theme(legend.position="bottom")
plota
dev.off()
