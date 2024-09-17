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
