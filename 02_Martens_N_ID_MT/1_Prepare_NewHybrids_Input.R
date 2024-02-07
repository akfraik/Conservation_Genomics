setwd("~/MT_ID")
library("adegenet")
library("graph4lg")
library(parallelnewhybrid)
source("/Users/alexandra.fraik/Documents/Other_Small/Marten_Helen_IDMT_2024/nh_plotR.R")
library(ggplot2)
library(tidyr)
library(dplyr)

## Ok, so that doesn't seem to work well, so we're going to try to make these files ourselves
file<-read.csv("/Users/alexandra.fraik/Documents/Other_Small/Marten_Helen_IDMT_2024/STRUCTURE_12loci_k2_admix.csv")
head(file)
meta<-file[colnames(file) %in% c("IndividualName","Marten.mtDNA","Pop","k1","k2")]

## Let's make a few different versions of this file for NewHybrids, with sample and locus names
no_prior<-file[c(1,9:23)]
head(no_prior)
genos<-no_prior[,c(2:16)]
genos[genos=="0/0"]<-"000/000"
genos[genos=="/0"]<-"/000"
genos[genos=="0/"]<-"000/"
test<- genos %>% 
  separate_wider_delim(colnames(genos), "/", names_sep=".")
head(test)
inds <- nchar(test$Mer041.2) == 2
genos<-lapply(test, function(x) replace(x, nchar(x)==2, paste0('0', x[nchar(x)==2])))
here<-data.frame(ID=no_prior[c(1)],do.call(cbind, genos))
test<-cbind(here[ 1 ],mapply(paste0, here[, seq(2, 31, 2)], "",here[, seq(3, 31, 2)]))

## Take a look at this file
head(test)
dim(test)

## On top of that, let's just try to get this ready as input- let's add a column of numbers and a column of "ns"
names<-rep("n",length(test$IndividualName))
number<-1:as.numeric(length(test$IndividualName))
no_priors<-data.frame(number,names,test)
write.table(no_priors,"No_Priors.txt",quote=FALSE,row.names=FALSE)
# save the numeric and sample ID index
no_priors_index<-no_priors[c(1,3)]

################################################################################################
## NOW, let's make a file WITH priors- one with the "s"
################################################################################################
meta$Prior<-rep("UNK",length(meta$IndividualName))
meta_p <- meta %>%
  mutate(Prior=case_when(
    (k1<0.95 & k2<0.95) ~ "UNK", k1>=0.95 ~" z1s",k2>=0.95 ~ " z0s"))
table(meta_p$Prior,meta$Marten.mtDNA)
#       americana caurina
#z0s        62       2
#z1s        32     193
#UNK        16      28
meta_priors<-meta_p[c("IndividualName","Prior")]
new<-merge(no_priors,meta_priors,by="IndividualName")
#priors<-new[c(2,1,3,19,4:18)]
#head(priors)
#write.table(priors,"Priors.txt",quote=FALSE,row.names=FALSE)
alt_priors<-new[c(2,19,4:18)]
head(alt_priors)
write.table(alt_priors,"Priors.txt",quote=FALSE,row.names=FALSE)
#priors_index<-priors[c(1,2)]

################################################################################################
## NOW, let's make a file WITH priors- one without the "s"
################################################################################################
meta$Prior<-rep("UNK",length(meta$IndividualName))
meta_p <- meta %>%
  mutate(Prior=case_when(
    (k1<0.95 & k2<0.95) ~ "UNK",
    k1>=0.95 ~" z1",
    k2>=0.95~ " z0"))
meta_prior<-meta_p[c("IndividualName","Prior")]
new<-merge(no_priors,meta_prior,by="IndividualName")
#prior<-new[c(2,1,3,19,4:18)]
#head(prior)
#write.table(prior,"Prior.txt",quote=FALSE,row.names=FALSE)
alt_prior<-new[c(2,19,4:18)]
head(alt_prior)
write.table(alt_prior,"Prior.txt",quote=FALSE,row.names=FALSE)
index<-new[c(1,2)]

## Meta merge to fix the index issue
meta<-merge(meta,index,by="IndividualName")
