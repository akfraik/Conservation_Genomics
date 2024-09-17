# Set working directory
setwd("/Users/alexandra.fraik)
# Load required R packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(graph4lg)

## Read in the metadata file and subset down to what you need
file2<-read.csv("Metadata.csv")
colnames(file2)[1]<-"ID"
file3<-file2[c("ID","Haplotype","Lat","Long","pop","Year_Bracket","colors","points","Country","group")]
dim(file3)
# Read in the micrsoat file
df<-read.csv("Microsat_df.csv")
test<-cbind(df[ 1 ],mapply(paste0, df[, seq(3, 31, 2)], ":",df[, seq(3, 31, 2)]))
file4<-merge(test,file2,by="ID")
genos<-file4[,c(2:16)]
genos[genos=="0:0"]<-NA
inds<-file4[,c(1)]
pop<-file4[,c(20)]
init<-df2genind(genos,sep=":",ind.names=inds,pop=pop,ploidy=2,NA.char="0:0",type="codom",strata=NULL,hierarchy=NULL)
# So, at this point, we have a file that has rows that match our file4 (metadata file)
strata(init)<-file4

## Now, we have a strata file appended in which every individual is the "population"
test_df<-setPop(init,~ID)
# Now, we're going to subset this to just the factors we're interested in
file5<-file4[c("ID","Haplotype","Lat","Long","pop","Year_Bracket","colors","points","Country","group","Year")]

## Now, let's subset this file to just those we know are genetically from the population 
# they were sampled in
file5<-read.csv("GSI_Genetic_Population.csv")
test_df<-test_df[rownames(test_df@tab) %in% file5$ID,]
test_df@strata<-droplevels(test_df@strata)
test_df<-setPop(test_df,~group)

## NeEstimator file for everyone for estimating LD using the genetic population assignment
#genind_to_genepop(test_df,output="Ne/All.txt")
# You should have 308 individuals now
file6<-file2[file2$ID %in% rownames(test_df@tab),]

## Now, let's recode the "year brackets" to generations- we're going to stagger by generation time
# the increments of time
file6$Year_Bracket<-recode(file6$Year_Bracket,
                            "Pre 2000"='0',
                            "2000-2004"='3',
                            "2005-2009"='4',
                            "2010-2014"='5',
                            "2015-2022"='7')

## Now, from within generation '1' we'll recode founders from the midwest as 0
# and in generation 6 we'll recode individuals from 2015-2018 as 5
test<-as.vector(cut(file6$Year,seq(1988,2000,4),include.lowest=TRUE))
file6$Year_Bracket<-ifelse(file6$Year <= 2000, test, file6$Year_Bracket)
table(file6$Year_Bracket,file6$Year)
testing <- file6 %>% 
  mutate(Year_Bracket=ifelse(Year_Bracket=="[1988,1992]","0",Year_Bracket)) %>%
  mutate(Year_Bracket=ifelse(Year_Bracket=="(1992,1996]","1",Year_Bracket)) %>%
  mutate(Year_Bracket=ifelse(Year_Bracket=="(1996,2000]","2",Year_Bracket)) %>%
  mutate(Year_Bracket=ifelse(Year=="2020" | Year=="2021" | Year=="2022","6",Year_Bracket))
table(testing$Year_Bracket)

## For now, let's just go ahead and do this this way
strata(test_df)<-testing
length(which(rownames(test_df@tab)==testing$ID)=="TRUE")
test_df<-setPop(test_df,~Year_Bracket)
test_df<-test_df[!rownames(test_df@tab) %in% c("mm20_15"),]

## Let's make some input files for NeEstimator for each population for each time frame
# here is an example for the Cabinets
cab_df<-test_df[test_df@strata$group %in% c("Cabinets"),]
cab_df<-setPop(cab_df,~Year_Bracket)
genind_to_genepop(cab_df,output="Ne/Cabinets_Time.txt")
