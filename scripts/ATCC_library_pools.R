##plots of library capture pool comparisons
#ATCC genomes updated version
#October 2024

#use bowtie2 data

rm(list=ls())

library(tidyverse)
library(ggpubr)

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

#normalise by both raw read count and genome length - should be the same as mean read depth

counts_reads_norm<-counts_reads %>%
  mutate(norm_counts1 = matched/QC_reads) %>%
  mutate(norm_counts2 = matched/QC_reads/length) %>%
  mutate(norm_counts3 = matched/length) %>%
  mutate(genome_structure = case_when((virus == "Human_adenovirus_40"| virus == "Human_betaherpesvirus") ~ "DNA",
                                      ,.default = "RNA"))

####compare library capture pool####

cols3<-c("#228833","#AA3377")

facet_names<-c("dedup_TE" = "Deduplicated",
               "nodedup_TE"= "Non-Deduplicated")

metadata3 <- metadata %>%
  select(Sample.ID,Pool.for.sequencing) %>%
  rename(Sample_id = Sample.ID) %>%
  rename(Pool = Pool.for.sequencing)

counts_pool<-left_join(counts_reads_norm,metadata3,by="Sample_id")

pool_readcounts<-counts_pool %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=virus,y=log(matched)))+
  geom_boxplot()+
  facet_grid(Pool~type)+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_text(size=10))+
  ylab("Log(Viral Reads)")+
  scale_color_manual(values=cols3,labels=c("Pool 1","Pool 2"))

pool_counts_sum<-counts_pool %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=as.character(Viral.load),y=log(matched),colour=Pool))+
  geom_boxplot()+
  facet_grid(~type,labeller=as_labeller(facet_names))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  ylab("Log(Viral Reads)")+
  scale_color_manual(values=cols3,labels=c("Pool 1","Pool 2"))

pool_readcounts_norm<-counts_pool %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=virus,y=log(norm_counts1)))+
  geom_boxplot()+
  facet_grid(Pool~type)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_text(size=10))+
  ylab("Log(viral reads/cleaned reads)")

pool_norm1_sum<-counts_pool %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=as.character(Viral.load),y=log(matched),colour=Pool))+
  geom_boxplot()+
  facet_grid(~type,labeller=as_labeller(facet_names))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  ylab("Log(viral reads/cleaned reads)")+
  scale_color_manual(values=cols3,labels=c("Pool 1","Pool 2"))

pool_readcounts_norm2<-counts_pool %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=virus,y=log(norm_counts2)))+
  geom_boxplot()+
  facet_grid(Pool~type)+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  ylab("Log(viral reads/cleaned reads/genome length)")

pool_norm2_sum<-counts_pool %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=as.character(Viral.load),y=log(matched),colour=Pool))+
  geom_boxplot()+
  facet_grid(~type,labeller=as_labeller(facet_names))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  ylab("Log(viral reads/cleaned reads/genome length)")+
  scale_color_manual(values=cols3,labels=c("Pool 1","Pool 2"))

#import and combine read depth files

depth_dedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_sample_and_virus_bowtie2_dedup_atcc_ref.tsv", sep = "\t", header = TRUE)

depth_nodedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_sample_and_virus_bowtie2_nodedup_atcc_ref.tsv", sep = "\t", header = TRUE)

depths_bt_all<-rbind(depth_dedup_bt,depth_nodedup_bt)

depths_reads<-left_join(depths_bt_all,metadata2,by="Sample_id")%>%
  mutate(genome_structure = case_when((virus == "Human_adenovirus_40"| virus == "Human_betaherpesvirus") ~ "DNA",
                                      ,.default = "RNA"))

depths_pool<-left_join(depths_reads,metadata3,by="Sample_id")

pool_readdepths<-depths_pool %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=virus,y=log(mean_depth)))+
  geom_boxplot()+
  facet_grid(Pool~type)+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  ylab("Log(mean read depth)")

pool_depths_sum<-depths_pool %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=as.character(Viral.load),y=log(mean_depth),colour=Pool))+
  geom_boxplot()+
  facet_grid(~type,labeller=as_labeller(facet_names))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  ylab("Log(mean read depth)")+
  scale_color_manual(values=cols3,labels=c("Pool 1","Pool 2"))

ggarrange(pool_readcounts,pool_readcounts_norm,pool_readcounts_norm2,pool_readdepths,nrow=2,ncol=2,heights=c(1,2),common.legend=TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/pool_compare_viruses.pdf")

ggarrange(pool_counts_sum,pool_norm1_sum,pool_norm2_sum,pool_depths_sum,nrow=2,ncol=2,common.legend = TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/pool_compare_viruses_sum.png",width=10,height=7)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/FigureS2.pdf",width=10,height=7)
