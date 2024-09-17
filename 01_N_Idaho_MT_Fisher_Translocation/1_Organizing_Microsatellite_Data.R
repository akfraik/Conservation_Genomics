setwd("/Users/alexandra.fraik/Documents/Fisher_Stuff/Cabinets_Genetics_2023/All_Plus_BC")
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
library(pegas)
library(graph4lg)

## We'll start with the output from the PCoA (where we recoded the number of pop/gen pops)
file3<-read.csv("/Users/alexandra.fraik/Documents/Fisher_Stuff/Cabinets_Genetics_2023/Input_Files/K4_Pops.csv")
colnames(file3)[1]<-"ID"
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
length(match(rownames(init@tab),file4$ID))
file5<-file4[c("ID","Haplotype","Lat","Long","pop","Year_Bracket","colors","points","Country","group")]

## For now, let's just go ahead and do this this way
strata(init)<-file5
test_df<-setPop(init,~ID)
table(test_df@pop,test_df@strata$pop)
length(match(rownames(test_df@tab),file4$ID))
# make sure this is 397
names(test_df@strata)
test_df<-setPop(test_df,~pop)

################################################################################################
## Look at missing data proportions
################################################################################################
## All of these loci are genotyped at > 98% of the loci
locmiss_test<-propTyped(test_df,by="loc")
pdf(file="Hierfstat/Missing_Data.pdf",width=7,height=7)
barplot(locmiss_test,ylim = c(0,1), ylab = "Complete genotypes (proportion)",xlab = "Locus", las = 2, cex.names = 0.7)
dev.off()
indmiss_test<-propTyped(test_df,by="ind")

## Now let's see if anything could be a duplicate
dups_test<-mlg.id(test_df)
# for each element in the list object
for (i in dups_test){
# if the length is greater than 1
  if (length(dups_test[i]) > 1){
    # print individuals that are duplicates
    print(i) } }
table(test_df$loc.fac)

## Print the number of private alleles per population
private_alleles<-private_alleles(test_df)
private_alleles<-colSums(t(private_alleles) != 0)
## Now, let's make a function to count them
barplot(private_alleles,ylab="Private Alleles",xlab = "Locus", las = 2, cex.names = 0.7)

## Let's make a figure that just shows these three pops in space
init.hfstat<-genind2hierfstat(test_df)
init.hfstat$ID<-row.names(init.hfstat)
test<-merge(init.hfstat,file5,by="ID")
names(test)
hfstat.man<-test[!names(test) %in% c("pop.x")]
names(hfstat.man)
colnames(hfstat.man)[20]<-"pop"

## Let's get the colors and pop order set up
pop_order<-c("Chetwyn","WilliamsLk","Midwest","Selkirks","Cabinets","Seeley_Swan",
             "StJoe","N_Lolo","Lolo","Clearwater","EF_Granite_Creek","Bitterroot",
             "Red_River","Salmon_Challis")
pca_df<-arrange(transform(hfstat.man,pop=factor(pop,levels=pop_order),pop))
colors<-as.character(pca_df$colors)
names(colors)<-as.character(pca_df$pop)

## Finally, let's get the haplotypes order and points set up
hap_order<-c("1b","4","6","7","11","9","12","1","5","10")
man<-arrange(transform(pca_df,Haplotype=factor(Haplotype,levels=hap_order),Haplotype))
points<-as.numeric(man$points)
names(points)<-as.character(man$Haplotype)
man<-droplevels(man)

## Plot these groups 1-3 
#pdf(file="Hierfstat/Map_Pops.pdf",width=8,height=8)
plot<-ggplot(man,aes(x=Long,y=Lat))
plot<-plot+geom_jitter(aes(color=pop,fill=pop),size=5,width=0.6,height=0.6)
plot<-plot+labs(x="Longitude",y="Latitude")
plot<-plot+scale_color_manual(name="Populations",values=colors)
plot<-plot+scale_fill_manual(name="Populations",values=colors)
plot<-plot+theme_bw()
plot<-plot+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plot<-plot+guides(fill="none",shape=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)),size="none")
plot<-plot+theme(legend.position="right")
plot
#dev.off()

## Convert to hierfstat file from a genind file
test_df<-setPop(test_df,~group)
hw.test(test_df, B = 10000)
init.hfstat<-genind2hierfstat(test_df)
table(init.hfstat$pop)
bs.fish<-basic.stats(init.hfstat,diploid=TRUE,digits=2)
Ho<-apply(bs.fish$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
boxplot(bs.fish$Ho)
He<-apply(bs.fish$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
Fis<-apply(bs.fish$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
betas<-pairwise.betas(init.hfstat,diploid=TRUE)
betas
## Seemes like it's just another estimator for pairwise FST...

## Count the number of alleles for each pop
allele_summ<-nb.alleles(init.hfstat)
colnames(allele_summ)<-levels(init.hfstat$pop)

## Plot
boxplot(bs.fish$Hs)
boxplot(bs.fish$Fis)
names(bs.fish)
bartlett.test(list(bs.fish$Hs,bs.fish$Ho))
barplot(betas(init.hfstat)$betaiovl)
# over loci of pairwise Fst
boot.ppfst(init.hfstat,nboot=1000)
boot.ppbetas(init.hfstat,nboot=1000)
boot.ppfis(init.hfstat,nboot=1000)

## Convert hierfstat object to dataframe for looking at pops, not genetic clusters
test_df<-setPop(test_df,~pop)
init.hfstat<-genind2hierfstat(test_df)
table(init.hfstat$pop)
bs.fish<-basic.stats(init.hfstat,diploid=TRUE,digits=2)

## Convert this into an hierfstat object from a genind object
Ho<-apply(bs.fish$Ho, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
boxplot(bs.fish$Ho)
He<-apply(bs.fish$Hs, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
Fis<-apply(bs.fish$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 2)
boxplot(bs.fish$Hs)
boxplot(bs.fish$Fis)
names(bs.fish)
bartlett.test(list(bs.fish$Hs,bs.fish$Ho))
barplot(betas(init.hfstat)$betaiovl)
betas(init.hfstat)
betas<-pairwise.betas(init.hfstat,diploid=TRUE)
betas(init.hfstat)$betaiovl
boot.ppbetas(init.hfstat,nboot=1000)
#We observed the greatest population-specific moments estimate of genetic differentiation 
 #(FST) between fishers sampled geographically from Adams Lake and all other pairs of sampling locations (Average βWTi= 0.40; Supplementary Table 2). 

# over loci of pairwise Fst
boot.ppfst(init.hfstat,nboot=1000)
# pairwise betas FST estimates

boot.ppfis(init.hfstat,nboot=1000) 

## Create a data.frame of site names, Ho and He and then convert to long format
Het_df<-data.frame(Genetic_Cluster=names(Ho),Ho=Ho,He=He) %>% 
  melt(id.vars="Genetic_Cluster")

## Let's try this allelic richness measure again
#pdf(file="Hierfstat/Heterozygosity_Genetic_Cluster.pdf",width=15,height=5)
plot<-ggplot(data=Het_df,aes(x=Genetic_Cluster,y=value,fill=variable))
plot<-plot+geom_bar(stat="identity",position=position_dodge(width = 0.6),colour="black")
plot<-plot+scale_fill_manual(values=c("royalblue","grey"))
plot<-plot+labs(x="Genetic Cluster",y="Heterozygosity")
plot<-plot+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plot<-plot+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-plot+theme(legend.position="bottom")
plot
#dev.off()

## Next, let's look at pairwise FST
fst<-genet.dist(init.hfstat,method="WC84") %>% 
  round(digits = 3)
# Set desired order of labels
lab_order<-pop_order
# Change order of rows and cols
fst.mat<-as.matrix(fst)
fst.mat1<-fst.mat[lab_order,]
fst.mat2<-fst.mat1[,lab_order]
# Convert the matrix into a data.frame
ind<-which(upper.tri(fst.mat2),arr.ind=TRUE)
fst.df<-data.frame(Cluster1=dimnames(fst.mat2)[[2]][ind[,2]],Cluster2=dimnames(fst.mat2)[[1]][ind[,1]],Fst=fst.mat2[ind])
fst.df$Cluster1<-factor(fst.df$Cluster1, levels = unique(fst.df$Cluster1))
fst.df$Cluster2<-factor(fst.df$Cluster2, levels = unique(fst.df$Cluster2))
fst.df$Fst[fst.df$Fst<0]=0

## Print data.frame summary and extract middle Fst value
mid<-max(fst.df$Fst) / 2
max<-max(fst.df$Fst)
# Plot FST heatmap
pdf(file="Hierfstat/FST_Cluster.pdf",width=10,height=10)
plot<-ggplot(data=fst.df,aes(x=Cluster1,y=Cluster2,fill=Fst))
plot<-plot+geom_tile(colour="black")
plot<-plot+geom_text(aes(label=Fst),color="black",size=5)
plot<-plot+scale_fill_gradient2(low="blue",mid="pink",high="red",midpoint=mid,limits=c(0,max),breaks=c(0,max))
plot<-plot+scale_x_discrete(expand=c(0,0))
plot<-plot+scale_y_discrete(expand=c(0,0),position="right")
plot<-plot+theme(axis.text.x=element_text(colour="black",size=15,angle=90,face="bold"),axis.title=element_blank())
plot<-plot+theme(axis.text=element_text(colour="black",size=15,face="bold"),axis.title=element_blank())
plot<-plot+theme(panel.grid = element_blank(),panel.background=element_blank())
plot
dev.off()

## But remember, there was some weirdness in the private alleles calculation, so maybe let's check
boxplot(t(allelic.richness(init.hfstat)$Ar))
boxplot(allelic.richness(init.hfstat)$Ar)
AR_df<-allelic.richness(init.hfstat)$Ar
summ<-colMeans(AR_df)
# Now, let's just make a dataframe for this
AR_df$Locus<-rownames(AR_df)
df <- AR_df %>% 
  gather(group,AR,-Locus)
head(df)
df$group<-as.factor(df$group)
df$Locus<-as.factor(df$Locus)

## Let's output a file for NeEstimator for genetic pops
genind_to_genepop(test_df,output="Ne/All.txt")

######## Mike also wants a PCoA that includes each individual (not by group) to see how it differs from PCA #######
test_df<-setPop(init,~ID)
Ds_pop<-mat_gen_dist(test_df,dist="PCA")
pcoa(as.matrix(Ds_pop))
Ds_mod<-as.dist(Ds_pop)

## Let's try to make a PCoA a different way
pcoa<-dudi.pco(Ds_mod,scannf=FALSE,nf=4) 
pcoa$eig
scatter(pcoa,xax=1,yax=2,clab.row=1,posieig="topright")
tab<-pcoa$tab
tab$ID<-rownames(tab)
merged<-merge(tab,file4,by="ID")

## Let's get the colors and pop order set up
colors<-as.character(merged$colors)
names(colors)<-as.character(merged$pop)
pop_order<-c("Chetwyn","WilliamsLk","Midwest","Selkirks","Cabinets","Seeley_Swan",
             "StJoe","N_Lolo","Lolo","Clearwater","EF_Granite_Creek","Bitterroot",
             "Red_River","Salmon_Challis")
pca_df<-arrange(transform(merged,pop=factor(pop,levels=pop_order),pop))

## Now, let's fix the groups
pca_df$group<-as.factor(pca_df$group)
orders<-c("Clearwaters","Cabinets","Chetwyn","WilliamsLk")
merged<-arrange(transform(pca_df,group=factor(group,levels=orders),group))
neworder<-c("Pre 2000","2000-2004","2005-2009","2010-2014","2015-2022")
merged<-arrange(transform(merged,Year_Bracket=factor(Year_Bracket,levels=neworder),Year_Bracket))

## Plot these groups 1-3: Axis 2 changes depending on if you use FST or DPS
pdf(file="Hierfstat/PCOA_PCA_Inds_12.pdf",width=8,height=6)
plot<-ggplot(merged,aes(x=A1,y=A2),group=pop)
plot<-plot+geom_jitter(aes(fill=pop,color=pop,shape=Haplotype),size=5)
plot<-plot+scale_color_manual(name="Populations",values=colors)
plot<-plot+scale_fill_manual(name="Populations",values=colors)
plot<-plot+scale_shape_manual(name="Haplotypes",values=points)
plot<-plot+labs(x="PCOA Axis 1",y="PCOA Axis 2")
plot<-plot+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plot<-plot+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-plot+theme(legend.position="right")
plot<-plot+guides(color="none",fill="none",size="none",shape=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)))
plot
dev.off()

## Plot these groups 1-3: Axis 2 changes depending on if you use FST or DPS
pdf(file="Hierfstat/PCOA_PCA_Inds_12_sans_hts.pdf",width=8,height=6)
plot<-ggplot(merged,aes(x=A1,y=A2),group=pop)
plot<-plot+geom_point(aes(fill=pop),color="black",size=5,width=0.3,height=0.3,shape=21)
plot<-plot+scale_fill_manual(name="Populations",values=colors)
plot<-plot+labs(x="PCOA Axis 1",y="PCOA Axis 2")
plot<-plot+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plot<-plot+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-plot+theme(legend.position="right")
plot
dev.off()

## Plot these groups 1-3: Axis 2 changes depending on if you use FST or DPS
pdf(file="Hierfstat/PCOA_PCA_Inds_34.pdf",width=6,height=6)
plot<-ggplot(merged,aes(x=A3,y=A4),group=pop)
plot<-plot+geom_jitter(aes(fill=pop,color=pop,shape=Haplotype),size=5)
plot<-plot+scale_color_manual(name="Populationss",values=colors)
plot<-plot+scale_fill_manual(name="Populations",values=colors)
plot<-plot+scale_shape_manual(name="Haplotypes",values=points)
plot<-plot+labs(x="PCOA Axis 3",y="PCOA Axis 4")
plot<-plot+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plot<-plot+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-plot+theme(legend.position="right")
plot<-plot+guides(color="none",fill="none",size="none",shape=guide_legend(ncol=1,byrow=FALSE,override.aes=list(size=5)))
plot
dev.off()

## Plot these axes 3 and 4 changes depending on if you use FST or DPS
pdf(file="Hierfstat/PCOA_PCA_Inds_34_sans_hts.pdf",width=6,height=6)
plot<-ggplot(merged,aes(x=A3,y=A4),group=pop)
plot<-plot+geom_point(aes(fill=pop),color="black",size=5,width=0.3,height=0.3,shape=21)
plot<-plot+scale_fill_manual(name="Populations",values=colors)
plot<-plot+labs(x="PCOA Axis 3",y="PCOA Axis 4")
plot<-plot+theme(axis.text.y=element_text(color="black",size=12),axis.title.y=element_text(size=15,face="bold"))
plot<-plot+theme(axis.text.x=element_text(color="black",size=12),axis.title.x=element_text(size=15,face="bold"))
plot<-plot+theme(strip.text=element_text(face="bold",size=15),strip.background=element_rect(color="black"),panel.border=element_rect(color="black",fill="NA"))
plot<-plot+theme(legend.position="right")
plot
dev.off()

## Analyse how much percent of genetic variance is explained by each axis
pca1<-pcoa$eig[1]
pcoa_eig<-pcoa$eig
percent<-pcoa_eig/sum(pcoa$eig)*100
percent

## Let's save the PCA plot of the variance explained
pdf(file="Hierfstat/PCoA_PCA_Variance_Explained.pdf",width=10,height=8)
plot<-barplot(percent,ylab="Genetic variance explained by eigenvectors (%)",ylim=c(0,12),names.arg=round(percent,1))
plot
dev.off()

## Let's output a file for NeEstimator for genetic pops
test_df<-setPop(init,~Country)
test_df<-test_df[test_df@strata$Country %in% c("USA"),]
genind_to_genepop(test_df,output="Ne/All_Country.txt")
