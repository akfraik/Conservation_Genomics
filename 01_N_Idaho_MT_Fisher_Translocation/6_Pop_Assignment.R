# Set working directory
setwd("/Users/alexandra.fraik")
# Load R packages
library(rubias)
library(tidyverse)
library(dplyr)

########################################################
## File generation for RUBIAS
########################################################
## Let's start off by reading in the massive file
ref<-read.csv("Tables/Pop_Summary_GSI_Ref.csv")
mix<-read.csv("Tables/Pop_Summary_Mixture.csv")
file<-rbind(ref,mix)
df<-read.csv("/Users/alexandra.fraik/Documents/Fisher_Stuff/Cabinets_Genetics_2023/Input_Files/Other_other_df.csv")
file4<-merge(df,file,by="ID")
file4<-file4[!colnames(file4) %in% c("group.1")]
dim(file4)
US_sub<-file4[file4$DAPC_Group %in% c("Cabinets","Clearwaters"),]
Canada<-file4[file4$DAPC_Group %in% c("Chetwynd","WilliamsLk"),]
US_sub$Country<-rep("USA",times=length(US_sub$ID))
Canada$Country<-rep("Canada",times=length(Canada$ID))
reference_df<-rbind(US_sub,Canada)

## Let's subset the US just for the early year brackets first
US_sub<-US_sub[US_sub$Year_Bracket %in% c("Pre 2000","2000-2004"),]
table(US_sub$pop,US_sub$DAPC_Group)
reference<-rbind(US_sub,Canada)

## Finally, let's drop singletons
reference<-reference[!reference$pop %in% c("Adamslake","burnslake","Selkirks","N_Lolo"),]
reference<-reference[!(reference$pop=="Chetwynd" & reference$DAPC_Group=="WilliamsLk"),]
table(reference$pop,reference$DAPC_Group)
reference<-reference[!(reference$pop=="Cabinets" & reference$DAPC_Group=="Clearwaters"),]
reference<-reference[!(reference$pop=="princegeorge" & reference$DAPC_Group=="Chetwynd"),]
reference<-reference[!reference$Type=="Mixture",]
dim(reference)

################################################################################################
############ Make a reference file
## Repunit = Country
## Collection = Sampling location
############ Include all individuals
################################################################################################
genos<-reference[,c(2:31)]
## Rename the loci
prefix<-"locus"
test<-rep(seq(1:15),each=2)
names<-c("a","b")
testa<-rep(names,times=15)
loci.names<-paste(prefix,test,names,sep="_")
colnames(genos)<-loci.names
genos[genos=="0"]<-NA

### Now, create new meatadata file
dim(reference)
sample_type<-as.character(rep("reference",times=147))
df_ref<-data.frame(sample_type=sample_type,repunit=reference$DAPC_Group,collection=reference$pop,indiv=reference$ID,genos)
df_ref$sample_type<-as.character(df_ref$sample_type)
df_ref$repunit<-as.character(df_ref$repunit)
df_ref$collection<-as.character(df_ref$collection)
df_ref$indiv<-as.character(df_ref$indiv)

## Assessing the self-assignments from the reference
## Log likelihood is the log probability of the fish's genotype given it is from the collection using the leave-out approach
## The posterior probability is that of assigning the fish to the collection given an equal prior on every collection in the reference
## Z-scores can be used to assess the distribution of the z-score statistic for fish from known reference populations
sa<-self_assign(reference=df_ref,gen_start_col=5)
head(sa,n=100)
max_log_likelihood <- sa %>%
  group_by(collection, repunit) %>%
  filter(log_likelihood == max(log_likelihood)) 
test<-as.data.frame(max_log_likelihood)
test<-test[c(1:8)]

### Population orders
newerorder<-c("Chetwynd","princegeorge","WilliamsLk","Midwest","Cabinets","Seeley_Swan",
              "Lolo","Clearwater","EF_Granite_Creek","Bitterroot")
length(newerorder)

## Now let's summarize the LOO approach by the individual mixture inferred collection 
## The scaled_likelihood is the posterior prob of assigning the fish to the inferred_collection given an equal prior on 
# every collection in the reference
sa_to_indiv<-sa %>%
  select(indiv,collection,inferred_collection,repunit,inferred_repunit,z_score,log_likelihood,scaled_likelihood)
write.csv(sa_to_indiv,"Assessment_GSI_References_Individuals.csv",row.names=FALSE)

## Now let's summarize the leave-one-out approach by the individual mixture inferred collection 
## The scaled likelihood is the likehood of assigning a fish back to the inferred collection
## The repu scaled likelihood is the summation of the likelihoods of assigning a fish back to the inferred repunit
sa_to_collec <- sa_to_indiv %>% 
  group_by(collection, inferred_collection) %>% 
  summarise(collec_mean_scaled_like = mean(scaled_likelihood)) %>% 
  spread(inferred_collection,collec_mean_scaled_like) %>% 
  arrange(factor(collection,levels=newerorder),collection) %>% 
  select(as.vector(newerorder))
write.csv(sa_to_collec,"Assessment_GSI_References_Collections.csv",row.names=FALSE)

## Now let's summarize the leave-one-out approach by the individual mixture inferred collection 
# The scaled likelihood is the likehood of assigning a fish back to the inferred collection
# The repu scaled likelihood is the summation of the likelihoods of assigning a fish back to the inferred repunit
sa_to_repu <- sa %>%
  group_by(indiv, collection, repunit, inferred_repunit) %>%
  summarise(repu_scaled_like = sum(scaled_likelihood))
write.csv(sa_to_repu,"Assessment_GSI_References_Repunit.csv",row.names=FALSE)

## By collection- takes the average of the scaled likelihoods (likelihood of membership in inferred collection) 
# per inferred rep unit
sa_to_collec <- sa %>%
  group_by(indiv, collection, repunit, inferred_repunit) %>%
  summarise(repu_scaled_like = sum(scaled_likelihood)) %>%
  group_by(collection, repunit, inferred_repunit) %>%
  summarise(repu_mean_scaled_like = mean(repu_scaled_like)) %>% 
  spread(inferred_repunit,repu_mean_scaled_like) %>% 
  arrange(factor(collection,levels=newerorder),collection)
write.csv(sa_to_collec,"Assessment_GSI_References_Inferred_Repunits.csv",row.names=FALSE)

## Let's also use the LOO approach to simulate different mixtures, now, call those repunits and then sum up over reporting units
# resampling_unit = "gene_copies"
omy_sims<-assess_reference_loo(reference=df_ref,gen_start_col=5,reps=500,mixsize=50,return_indiv_posteriors=TRUE)
mix_props<-omy_sims$mixing_proportions
tmp <- mix_props %>%
  group_by(iter, repunit) %>%
  summarise(true_repprop = sum(true_pi), 
            reprop_posterior_mean = sum(post_mean_pi),
            repu_n = sum(n)) %>%
  mutate(repu_n_prop = repu_n / sum(repu_n))

## So, really I need to determine how many SDs from the mean a fisher has to be (for its PofZ) 
# I am going to compare the distributions for the number of fish that have a PofZ for their max posterior membership
sa_to_repu_z_scores <- sa %>%
  group_by(indiv, collection, inferred_collection, repunit, inferred_repunit, z_score) %>%
  summarise(repu_scaled_like = sum(scaled_likelihood))
write.csv(sa_to_repu_z_scores,"Individual_Posterior.csv",row.names=FALSE)

################################################################################################################
### Plot the mixutre and simulated proportions per repunit
################################################################################################################
pdf(file="Simulations_Genetic_Population_Reassign.pdf",width=8,height=8)
plota<-ggplot(tmp,aes(x=true_repprop,y=reprop_posterior_mean,colour=repunit))+geom_point() 
plota<-plota+scale_color_manual(name="Reference Reporting Units",values=c("aquamarine2","darkmagenta"))
plota<-plota+geom_abline(intercept=0,slope=1)
plota<-plota+facet_wrap(~repunit)
plota<-plota+labs(x="True Posterior Mixture Proprtion",y="Simulated Posterior Mixture Proprtion")
plota<-plota+scale_x_continuous(limits=c(0,1))
plota<-plota+scale_y_continuous(limits=c(0,1))
plota<-plota+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plota<-plota+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plota<-plota+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plota<-plota+theme(legend.text=element_text(size=rel(1)))
plota<-plota+theme(legend.title=element_text(face="bold",size=rel(1.5)))
plota<-plota+theme(legend.position="none")
plota
dev.off()

### Simulated_repunit is the reporting unit the individual was simulated from 
### Simulated_collection is the collection the simulated genotype came from 
### PofZ is the mean posterior membership over the MCMC of the posterior probability that the individual 
# originated from the COLLECTION, not repunit
indiv_props<-omy_sims$indiv_posteriors

## We can see what the distribution of posteriors to the correct reporting unit is for fisher from the different 
# simulated collections
## We’ll do that with a boxplot, coloring by repunit, first aggregating over reporting units
repu_pofzs <- indiv_props %>%
  filter(repunit == simulated_repunit) %>%
  group_by(iter, indiv, simulated_collection, repunit) %>%  
  summarise(repu_PofZ = sum(PofZ)) %>%
  ungroup() %>%
  arrange(repunit, simulated_collection) %>%
  mutate(simulated_collection = factor(simulated_collection, levels = unique(simulated_collection)))
df<-arrange(transform(repu_pofzs,simulated_collection=factor(simulated_collection,levels=newerorder),simulated_collection))
repu_pofzs<- df %>% 
  drop_na(repu_PofZ)

# also get the number of simulated individuals from each collection
num_simmed <- indiv_props %>%
  group_by(iter, indiv) %>%
  slice(1) %>%
  ungroup() %>%
  count(simulated_collection)

## Set the colors for the collections
colors_sites<-c("aquamarine4","deepskyblue4","darkturquoise","goldenrod","darkblue","cornflowerblue","deeppink4",
                             "violetred","brown1","pink")
names<-c("Chetwynd","princegeorge","WilliamsLk","Midwest","Cabinets","Seeley_Swan","Lolo",
         "Clearwater","EF_Granite_Creek","Bitterroot")
names(colors_sites)<-names

################################################################################################################
### Plot the mixutre and simulated proportions per individual per collection
# note, the last few steps make simulated collection a factor so that collections within
# the same repunit are grouped together in the plot.
################################################################################################################
pdf(file="Simulations_Reference_Collections.pdf",width=8,height=8)
plota<-ggplot(repu_pofzs, aes(x = simulated_collection, y = repu_PofZ))
plota<-plota+geom_boxplot(aes(colour = repunit))
plota<-plota+labs(x="Simulated Sampling Location",y="Posterior Probability of Reassignment")
plota<-plota+geom_text(data = num_simmed, mapping = aes(y = 1.025, label = n), angle = 90, hjust = 0, vjust = 0.5, size = 3)
plota<-plota+theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9, vjust = 0.5))
plota<-plota+ylim(c(0, 1.05))
plota<-plota+scale_color_manual(name="Reference Reporting Units",values=c("violet","darkblue","darkturquoise","aquamarine4"))
plota<-plota+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plota<-plota+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plota<-plota+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plota<-plota+theme(legend.text=element_text(size=rel(1)))
plota<-plota+theme(legend.title=element_text(face="bold",size=rel(1.5)))
plota<-plota+theme(legend.position="none")
plota
dev.off()

########################################################
## Mixture file
########################################################
## Make a reference file
mixture<-reference_df[!reference_df$ID %in% reference$ID,]
dim(mixture)
meta<-mixture[c("ID","pop","Year_Bracket","Haplotype","Country")]
genos<-mixture[,c(2:31)]
genos[genos=="0:0"]<-NA

## Rename the loci
prefix<-"locus"
test<-rep(seq(1:15),each=2)
names<-c("a","b")
testa<-rep(names,times=15)
loci.names<-paste(prefix,test,names,sep="_")
colnames(genos)<-loci.names

### Now, create a mixture file
repunit<-as.character(rep("here",times=length(mixture$ID)))
sample_type<-as.character(rep("mixture",times=length(mixture$ID)))
df_mix<-data.frame(sample_type=sample_type,repunit=repunit,collection=file5$Haplotype,
                   indiv=file5$ID,genos)
df_mix$sample_type<-as.character(df_mix$sample_type)
df_mix$repunit<-as.character(df_mix$repunit)
df_mix$collection<-as.character(df_mix$collection)
df_mix$indiv<-as.character(df_mix$indiv)
df_mix[df_mix=="here"]<-NA

## For the first round, let's try this
## Method "PB" runs the full Bayesian model and looks for 
mix_est<-infer_mixture(reference=df_ref,mixture=df_mix,gen_start_col=5,
                       burn_in=100,pb_iter=100,method="PB")
lapply(mix_est, head)
location_post<-as.data.frame(mix_est$indiv_posteriors)
location_post<-location_post[c(1:7)]
write.csv(location_post,"Individual_Posterior_Location.csv",row.names=FALSE)

## Aggregating collections into reporting units for mixing proportions 
# adding mixing proportions over collections in the repunit
rep_mix_ests <- mix_est$mixing_proportions %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = sum(pi))  
write.csv(rep_mix_ests,"Rep_Mix_Ests_Location.csv",row.names=FALSE)

## For individuals posteriors
rep_indiv_ests <- mix_est$indiv_posteriors %>%
  group_by(mixture_collection, indiv, repunit) %>%
  summarise(rep_pofz = sum(PofZ))
colnames(rep_indiv_ests)[2]<-"sample.id"
write.csv(rep_indiv_ests,"Rep_Indiv_Ests_Location.csv",row.names=FALSE)

####################################################################################################
#### Let's plot the output 
## Let's start by reading in the DAPC and the STRUCTURE files
####################################################################################################
mix<-read.csv("Tables/Pop_Summary_Mix.csv")
file<-rbind(ref,mix)
df<-read.csv("/Users/alexandra.fraik/Documents/Fisher_Stuff/Cabinets_Genetics_2023/Input_Files/Other_other_df.csv")
file4<-merge(df,file,by="ID")
file4<-file4[!colnames(file4) %in% c("group.1")]
US_sub<-file4[file4$DAPC_Group %in% c("Cabinets","Clearwaters"),]
Canada<-file4[file4$DAPC_Group %in% c("Chetwynd","WilliamsLk"),]
US_sub$Country<-rep("USA",times=length(US_sub$ID))
Canada$Country<-rep("Canada",times=length(Canada$ID))
reference_df<-rbind(US_sub,Canada)

######## So, let's start by just plotting the weird ones (that didn't have concordant DAPC and STRUCTURE results)
mixed<-reference_df[reference_df$GSI=="Mixture",]
dim(mixed)
sub1_df<-sub1[sub1$sample.id %in% mixed$ID,]
dim(sub1_df)
merged_1<-merge(sub1_df,mixed,by.x="sample.id",by.y="ID",all.x=TRUE)
dim(merged_1)

####### STOP- we need to come back and add some lines of code about posterior probabilities to
## So NOW, let's try to plot
pdf(file="RUBIAS/Sub1/Mixture_Sims_Country.pdf",height=8,width=6)
plota<-ggplot(merged_1,aes(fill=repunit,y=rep_pofz,x=sample.id),group=pop)
plota<-plota+geom_bar(position="fill",stat="identity")
plota<-plota+scale_fill_manual(values=c("aquamarine2","darkmagenta"))
plota<-plota+ylab(label="Inferred Genetic Ancestries")
plota<-plota+xlab(label="Individuals")
plota<-plota+theme_bw()
plota<-plota+facet_grid(.~pop,switch="x",scales="free_x",space="free_x")
plota<-plota+theme(strip.text.x=element_text(color="black",size=12,angle=90),axis.text.x=element_blank())
plota<-plota+theme(axis.text.y=element_text(color="black",size=14,face="bold"),axis.title.y=element_text(color="black",size=16,face="bold"))
plota
dev.off()

## Subset to just the American ones (that we care about)
merged_11<-merged_1[merged_1$Country=="USA",]

## Let's also just plot the unknown American ones
## So NOW, let's try to plot
pdf(file="RUBIAS/Sub1/Mixture_Sims_USA_Country.pdf",height=8,width=6)
plota<-ggplot(merged_11,aes(fill=repunit,y=rep_pofz,x=sample.id),group=pop)
plota<-plota+geom_bar(position="fill",stat="identity")
plota<-plota+scale_fill_manual(values=c("aquamarine2","darkmagenta"))
plota<-plota+ylab(label="Inferred Genetic Ancestries")
plota<-plota+xlab(label="Individuals")
plota<-plota+theme_bw()
plota<-plota+facet_grid(.~pop,switch="x",scales="free_x",space="free_x")
plota<-plota+theme(strip.text.x=element_text(color="black",size=12,angle=90),axis.text.x=element_blank())
plota<-plota+theme(axis.text.y=element_text(color="black",size=14,face="bold"),axis.title.y=element_text(color="black",size=16,face="bold"))
plota
dev.off()

######## So, let's start by just plotting the weird ones (that didn't have concordant DAPC and STRUCTURE results)
mixed<-reference_df[reference_df$GSI=="Mixture",]
dim(mixed)
sub2_df<-sub2[sub2$sample.id %in% mixed$ID,]
dim(sub2_df)
merged_2<-merge(sub2_df,mixed,by.x="sample.id",by.y="ID",all.x=TRUE)
dim(merged_2)

## So NOW, let's try to plot
pdf(file="RUBIAS/Sub2/Mixture_Sims_GenPop.pdf",height=8,width=10)
plota<-ggplot(merged_2,aes(fill=repunit,y=rep_pofz,x=sample.id),group=pop)
plota<-plota+geom_bar(position="fill",stat="identity")
plota<-plota+scale_fill_manual(values=c("darkblue","aquamarine4","violet","darkturquoise"))
plota<-plota+ylab(label="Inferred Genetic Ancestries")
plota<-plota+xlab(label="Individuals")
plota<-plota+theme_bw()
plota<-plota+facet_grid(.~pop,switch="x",scales="free_x",space="free_x")
plota<-plota+theme(strip.text.x=element_text(color="black",size=12,angle=90),axis.text.x=element_blank())
plota<-plota+theme(axis.text.y=element_text(color="black",size=14,face="bold"),axis.title.y=element_text(color="black",size=16,face="bold"))
plota
dev.off()
write.csv(merged_2,"RUBIAS/Sub2/Summary_GSI.csv")

## Subset to just the American ones (that we care about)
merged_22<-merged_2[merged_2$Country=="USA",]

## Let's also just plot the unknown American ones
## So NOW, let's try to plot
pdf(file="RUBIAS/Sub2/Mixture_Sims_USA_Genetic_Population.pdf",height=8,width=6)
plota<-ggplot(merged_22,aes(fill=repunit,y=rep_pofz,x=sample.id),group=pop)
plota<-plota+geom_bar(position="fill",stat="identity")
plota<-plota+scale_fill_manual(values=c("darkblue","aquamarine4","violet","darkturquoise"))
plota<-plota+ylab(label="Inferred Genetic Ancestries")
plota<-plota+xlab(label="Individuals")
plota<-plota+theme_bw()
plota<-plota+facet_grid(.~pop,switch="x",scales="free_x",space="free_x")
plota<-plota+theme(strip.text.x=element_text(color="black",size=12,angle=90),axis.text.x=element_blank())
plota<-plota+theme(axis.text.y=element_text(color="black",size=14,face="bold"),axis.title.y=element_text(color="black",size=16,face="bold"))
plota
dev.off()

######## So, let's start by just plotting the weird ones (that didn't have concordant DAPC and STRUCTURE results)
mixed<-reference_df[reference_df$GSI=="Mixture",]
dim(mixed)
Sub3_df<-sub3[sub3$sample.id %in% mixed$ID,]
dim(Sub3_df)
merged_2<-merge(Sub3_df,mixed,by.x="sample.id",by.y="ID",all.x=TRUE)
dim(merged_2)

## So NOW, let's try to plot
pdf(file="RUBIAS/Sub3/Mixture_Sims_GenPop.pdf",height=8,width=10)
plota<-ggplot(merged_2,aes(fill=repunit,y=rep_pofz,x=sample.id),group=pop)
plota<-plota+geom_bar(position="fill",stat="identity")
plota<-plota+scale_fill_manual(values=c("darkblue","aquamarine4","violet","darkturquoise"))
plota<-plota+ylab(label="Inferred Genetic Ancestries")
plota<-plota+xlab(label="Individuals")
plota<-plota+theme_bw()
plota<-plota+facet_grid(.~pop,switch="x",scales="free_x",space="free_x")
plota<-plota+theme(strip.text.x=element_text(color="black",size=12,angle=90),axis.text.x=element_blank())
plota<-plota+theme(axis.text.y=element_text(color="black",size=14,face="bold"),axis.title.y=element_text(color="black",size=16,face="bold"))
plota
dev.off()
write.csv(merged_2,"RUBIAS/Sub3/Summary_GSI.csv")

## Subset to just the American ones (that we care about)
merged_22<-merged_2[merged_2$Country=="USA",]

## Let's also just plot the unknown American ones
## So NOW, let's try to plot
pdf(file="RUBIAS/Sub3/Mixture_Sims_USA_Genetic_Population.pdf",height=8,width=6)
plota<-ggplot(merged_22,aes(fill=repunit,y=rep_pofz,x=sample.id),group=pop)
plota<-plota+geom_bar(position="fill",stat="identity")
plota<-plota+scale_fill_manual(values=c("darkblue","aquamarine4","violet","darkturquoise"))
plota<-plota+ylab(label="Inferred Genetic Ancestries")
plota<-plota+xlab(label="Individuals")
plota<-plota+theme_bw()
plota<-plota+facet_grid(.~pop,switch="x",scales="free_x",space="free_x")
plota<-plota+theme(strip.text.x=element_text(color="black",size=12,angle=90),axis.text.x=element_blank())
plota<-plota+theme(axis.text.y=element_text(color="black",size=14,face="bold"),axis.title.y=element_text(color="black",size=16,face="bold"))
plota
dev.off()
