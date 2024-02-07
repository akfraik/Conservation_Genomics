setwd("~/MT_ID")
library("adegenet")
library("graph4lg")
source("/Users/alexandra.fraik/Documents/Other_Small/Marten_Helen_IDMT_2024/nh_plotR.R")
library(ggplot2)
library(tidyr)
library(dplyr)

#########################################################################################################################################################################
############################################### Let's start with just the individuals Helen cares about most for the SDMs ############################################### 
#########################################################################################################################################################################
file<-read.csv("~/STRUCTURE_12loci_k2_admix.csv")
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
              
##########################################################################################################################################################
## We don't have any putative individuals out of range that are "pure" in this data set, so we'll just make files with priors from STRUCTURE
## You could also use the mtDNA assignments...              
##########################################################################################################################################################
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

##########################################################################################################################################################
## Run the above again, but subset the genotype file using the code below to generate the input files (we'll call them sub form here on out)
##########################################################################################################################################################
file<-read.csv("~/STRUCTURE_12loci_k2_admix.csv")
# Helen said for the most part she is not using loci Ma3, Ma18, and Ma19 in analyses, so let's remove them
file<-file[!colnames(file) %in% c("Ma3","Ma18","Ma19","Missing")]

## Let's set up the metadata file
meta<-file[colnames(file) %in% c("IndividualName","Marten.mtDNA","Pop","k1","k2")]
meta$Prior<-rep("UNK",length(meta$IndividualName))

## Let's make a few different versions of this file for NewHybrids, with sample and locus names
no_prior<-file[c(1,8:19)]
head(no_prior)
genos<-no_prior[,c(2:13)]
genos[genos=="0/0"]<-"000/000"
genos[genos=="/0"]<-"/000"
genos[genos=="0/"]<-"000/"
test<- genos %>% 
  separate_wider_delim(colnames(genos), "/", names_sep=".")
head(test)
inds <- nchar(test$Mer041.2) == 2
genos<-lapply(test, function(x) replace(x, nchar(x)==2, paste0('0', x[nchar(x)==2])))
here<-data.frame(ID=no_prior[c(1)],do.call(cbind, genos))
test<-cbind(here[ 1 ],mapply(paste0, here[, seq(2, 25, 2)], "",here[, seq(3, 25, 2)]))
# Just make sure everything is looking good- should have 333 individuals and 12 loci
head(test)
dim(test)

## On top of that, let's just try to get this ready as input- let's add a column of numbers and a column of "ns"
names<-rep("n",length(test$IndividualName))
number<-1:as.numeric(length(test$IndividualName))
no_priors<-data.frame(number,names,test)
write.table(no_priors,"No_Priors.txt",quote=FALSE,row.names=FALSE)
# save the numeric and sample ID index
no_priors_index<-no_priors[c(1,3)]
### DO the same as above with the priors from STRUCTURE if you like (repeate lines 43-60)
              
##########################################################################################################################################################
## Now for the reference caurina and americana individuals from not nearby (outside of ID and MT), with the subset of loci that overlap
##########################################################################################################################################################
file<-read.csv("~/pure_populations_adegenet_01112024.csv")
## Let's also subset by those loci that we did in the previous file
file<-file[!colnames(file) %in% c("Ma3","Ma18","Ma19")]
head(file)

## Get your metadata file ready
meta<-file[c(1:4)]
meta$Prior<-rep("UNK",length(meta$SampleID))
table(meta$Prior)

## Now, get the body of your NewHybrids file ready
no_prior_all<-file[c(1,6:15)]
dim(no_prior_all)
genos<-no_prior_all[,c(2:11)]
genos[genos=="0/0"]<-"000/000"
genos[genos=="/0"]<-"/000"
genos[genos=="0/"]<-"000/"
test <- genos %>% 
  separate_wider_delim(colnames(genos), "/", names_sep=".")
head(test)
genos<-lapply(test, function(x) replace(x, nchar(x)==2, paste0('0', x[nchar(x)==2])))
here<-data.frame(ID=no_prior_all[c(1)],do.call(cbind, genos))
test<-cbind(here[ 1 ],mapply(paste0, here[, seq(2, 21, 2)], "",here[, seq(3, 21, 2)]))

## On top of that, let's just try to get this ready as input- let's add a column of numbers and a column of "ns"
names<-rep("n",length(test$SampleID))
number<-1:as.numeric(length(test$SampleID))
no_priors_all<-data.frame(number,names,test)
write.table(no_priors_all,"No_Priors.txt",quote=FALSE,row.names=FALSE)

################################################################################################
## Priors based on mtDNA assignment (no s since we want them inlucded in our output)
################################################################################################
meta_p <- meta %>%
  mutate(Prior=case_when(mtDNA=="caurina" ~ "z1", mtDNA=="americana" ~ "z0"))
meta_priors<-meta_p[c("SampleID","Prior")]
new<-merge(no_priors_all,meta_priors,by="SampleID")
priors_t<- new %>% 
  arrange(number)
alt_priors<-priors_t[c(2,14,4:13)]
head(alt_priors)

## Just make sure everything looks good- you should have 333 individuals that are "unknown" and 
# knowns coming from the reference data set
table(alt_priors$Prior)
write.table(alt_priors,"Priors.txt",quote=FALSE,row.names=FALSE)
priors_index<-priors_t[c(1,2)]
              
setwd("/Users/alexandra.fraik/Documents/Other_Small/Marten_Helen_IDMT_2024/All")
library("adegenet")
library("graph4lg")
library(parallelnewhybrid)
source("/Users/alexandra.fraik/Documents/Other_Small/Marten_Helen_IDMT_2024/nh_plotR.R")
library(ggplot2)
library(tidyr)
library(dplyr)

################################################################################################################
#### Finally, we'll use the combination of files and overlapping genotypes from the reference and MT/ID
#### Martens to test- both without priors and with priors             
################################################################################################################
## Ok, so that doesn't seem to work well, so we're going to try to make these files
file<-read.csv("~/all_combo_adegenet_020424.csv")
file<-file[!colnames(file) %in% c("Ma3","Ma18","Ma19")]
head(file)

## Get your metadata file ready
meta<-file[c(1:4)]
meta$Prior<-recode(meta$Prior,"n"='UNK')
table(meta$Prior)

## Now, get the body of your genotype file ready
no_prior_all<-file[c(1,5:12)]
head(no_prior_all)
genos<-no_prior_all[,c(2:9)]
genos[genos=="0/0"]<-"000/000"
genos[genos=="/0"]<-"/000"
genos[genos=="0/"]<-"000/"
test <- genos %>% 
  separate_wider_delim(colnames(genos), "/", names_sep=".")
head(test)
genos<-lapply(test, function(x) replace(x, nchar(x)==2, paste0('0', x[nchar(x)==2])))
here<-data.frame(ID=no_prior_all[c(1)],do.call(cbind, genos))
test<-cbind(here[ 1 ],mapply(paste0, here[, seq(2, 17, 2)], "",here[, seq(3, 17, 2)]))

## On top of that, let's just try to get this ready as input- let's add a column of numbers and a column of "ns"
names<-rep("n",length(test$SampleID))
number<-1:as.numeric(length(test$SampleID))
no_priors_all<-data.frame(number,names,test)
write.table(no_priors_all,"No_Priors.txt",quote=FALSE,row.names=FALSE)
# save the numeric and sample ID index
no_priors_index<-no_priors_all[c(1,3)]
head(no_priors_all)

################################################################################################
## NOW, let's make a file WITH priors- one without the "s"- just to see
################################################################################################
meta_p <- meta %>%
  filter(Prior=="y") %>%
  mutate(Prior=case_when(mtDNA=="caurina" ~ "z1", mtDNA=="americana" ~ "z0")) %>%
  rbind(meta %>% filter(Prior == "UNK"))
meta_priors<-meta_p[c("SampleID","Prior")]
new<-merge(no_priors_all,meta_priors,by="SampleID")
priors_t<- new %>% 
  arrange(number)
alt_priors<-priors_t[c(2,12,5:11)]

## Just make sure everything looks good- you should have 333 individuals that are "unknown" and 
# knowns coming from the reference data set
head(alt_priors)
write.table(alt_priors,"Priors.txt",quote=FALSE,row.names=FALSE)
priors_index<-priors_t[c(1,2)]

#######################################################################################################
## NOW, let's make a file WITH priors- one with the "s" (for the reference individuals)
## this is the output we'll probably care about comparing to the "no prior" just MT and ID individuals
#######################################################################################################
meta_p <- meta %>%
  filter(Prior=="y") %>%
  mutate(Prior=case_when(mtDNA=="caurina" ~ "z1", mtDNA=="americana" ~ "z0")) %>%
  rbind(meta %>% filter(Prior == "UNK"))
table(meta_p$Prior)
meta_priors<-meta_p[c("SampleID","Prior")]
new<-merge(no_priors_all,meta_priors,by="SampleID")
priors_t<- new %>% 
  arrange(number)
head(priors_t)
alt_priors<-priors_t[c(2,12,5:11)]

## Just make sure everything looks good- you should have 333 individuals that are "unknown" and 
# knowns coming from the reference data set
table(alt_priors$Prior)
head(alt_priors)
write.table(alt_priors,"Prior.txt",quote=FALSE,row.names=FALSE)
