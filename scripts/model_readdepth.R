#Linear model for TE data
#November 2024

#viral read depth
#use bowtie2 data

rm(list=ls())

library(tidyverse)
library(ggpubr)
library(lme4)
library(boot)
library(car)

####read counts/normalised read counts####

#metadata

metadata<-read.csv("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/metadata/sampleIDs_TESpikeIn.csv",header=TRUE)

metadata2 <- metadata %>%
  select(Sample.ID,Number.of.read.pairs..quality.adaptor.trimmed.) %>%
  rename(Sample_id = Sample.ID) %>%
  rename(QC_reads = Number.of.read.pairs..quality.adaptor.trimmed.)

#import and combine read count files
#deduplicated and non-deduplicated

counts_dedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readcount_per_sample_and_virus_bowtie2_dedup_atcc_ref.tsv", sep = "\t", header = TRUE)

counts_nodedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readcount_per_sample_and_virus_bowtie2_nodedup_atcc_ref.tsv", sep = "\t", header = TRUE)

counts_bt_all<-rbind(counts_dedup_bt,counts_nodedup_bt)

counts_reads<-left_join(counts_bt_all,metadata2,by="Sample_id")

####read depths####

#import and combine read depth files

depth_dedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_sample_and_virus_bowtie2_dedup_atcc_ref.tsv", sep = "\t", header = TRUE)

depth_nodedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_sample_and_virus_bowtie2_nodedup_atcc_ref.tsv", sep = "\t", header = TRUE)

depths_bt_all<-rbind(depth_dedup_bt,depth_nodedup_bt)

depths_reads<-left_join(depths_bt_all,metadata2,by="Sample_id")%>%
  mutate(genome_structure = case_when((virus == "Human_adenovirus_40"| virus == "Human_betaherpesvirus") ~ "DNA",
                                      ,.default = "RNA"))


####model for spike in viruses - mean read depth####
#read count and read depth are more similar to count data (but read depth not an integer so gives warnings)

depths_reads_sub<-depths_reads %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(type == "dedup_TE") %>%
  mutate(log_depth = log(mean_depth)) %>%
  mutate(log_viral_load = log10(Viral.load))

#inspect data
hist(depths_reads_sub$mean_depth)
hist(depths_reads_sub$Viral.load)
hist(depths_reads_sub$log_depth)
hist(depths_reads_sub$log_viral_load)

#lm1<-glm(mean_depth~Viral.load*virus,data=depths_reads_sub,family=poisson)
#similar to count data but gives warning because there are non-integers

lm1<-glm(mean_depth~Viral.load+virus,data=depths_reads_sub,family=gaussian(link="log"))

lm2<-glm(log_depth~log_viral_load*virus,data=depths_reads_sub,family=gaussian)

summary(lm1)

summary(lm2)

glm.diag.plots(lm1)

glm.diag.plots(lm2)

pred<-predict(lm1,type="response",se.fit=TRUE)

pdf<-data.frame(depths_reads_sub$Viral.load)
pdf$fit<-(pred$fit)^10
pdf$se<-(pred$se.fit)^10

ggplot()+
  geom_point(data=depths_reads_sub,aes(x=Viral.load,y=mean_depth,color=virus))+
  geom_line(data=pdf,aes(y=fit,x=depths_reads_sub.Viral.load))+
  geom_ribbon(data=pdf,aes(y=fit,ymin=fit-1.96*se,ymax=fit+1.96*se,x=depths_reads_sub.Viral.load),fill="blue",alpha=0.3)+
  scale_x_log10()+
  scale_y_log10()

#residuals plot shows funnel shape - homogeneity of variance not met. need to check how to resolve that

lm2<-lmer(mean_depth~Viral.load+ (1|virus),data=depths_reads_sub)

summary(lm3)
Anova(lm3)

#could also try with read count which is more straightforward to model potentially

counts_reads_sub<-counts_reads %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(type == "dedup_TE")

hist(counts_reads_sub$matched)

lm4<-glm(matched ~ Viral.load+virus,data=counts_reads_sub,family=poisson)

summary(lm4)

glm.diag.plots(lm4)

lm5<-glm(matched ~ Viral.load+virus,data=counts_reads_sub)

summary(lm5)

glm.diag.plots(lm5)
