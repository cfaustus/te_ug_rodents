##plots of read counts and viral read counts
#November 2024

rm(list=ls())

library(tidyverse)
library(ggpubr)

#metadata

metadata<-read.csv("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/metadata/sampleIDs_TESpikeIn.csv",header=TRUE)

metadata2 <- metadata %>%
  select(Sample.ID,Number.of.read.pairs..quality.adaptor.trimmed.) %>%
  rename(Sample_id = Sample.ID) %>%
  rename(QC_reads = Number.of.read.pairs..quality.adaptor.trimmed.)

#import read count files

counts_dedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readcount_per_sample_and_virus_bowtie2_dedup_atcc_ref.tsv", sep = "\t", header = TRUE)

counts_nodedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readcount_per_sample_and_virus_bowtie2_nodedup_atcc_ref.tsv", sep = "\t", header = TRUE)

#import viral reads mapped (to calculate proportions)

viral_reads_dedup<-read.csv("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/total_virus_mapped_reads_per_sample_dedup_atcc_ref_20241108.csv",header=TRUE)

viral_reads_nodedup<-read.csv("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/total_virus_mapped_reads_per_sample_nodedup_atcc_ref_20241108.csv",header=TRUE)

cols2<-c("#BB5566","#004488")

facet_names<-c("dedup_TE" = "Deduplicated",
               "nodedup_TE"= "Non-Deduplicated")

reads_metadata_dedup<-left_join(counts_dedup_bt,metadata2,by="Sample_id")

reads_viral_dedup<-left_join(reads_metadata_dedup,viral_reads_dedup,by="Sample_id") 

reads_metadata_nodedup<-left_join(counts_nodedup_bt,metadata2,by="Sample_id")

reads_viral_nodedup<-left_join(reads_metadata_nodedup,viral_reads_nodedup,by="Sample_id") 

reads_viral_all<-rbind(reads_viral_dedup,reads_viral_nodedup)

reads_plot<-reads_viral_all %>%
  group_by(Background,Sample_id,Viral.load,type) %>%
  summarise(total_reads = (QC_reads*2),
            viral_reads = total_virus_reads,
            ATCC_reads = sum(matched),
            prop_ATCC = ATCC_reads/total_reads,
            prop_viral = total_virus_reads/total_reads,
            diff = viral_reads - ATCC_reads) %>%
  unique()

total_reads<-reads_plot %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=log(total_reads),colour=Background))+
  geom_point()+
  geom_smooth(aes(group=Background),se=FALSE)+
  facet_grid(~type,labeller=as_labeller(facet_names))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ggtitle("Total Reads")+
  ylab("Log(Total Reads)")+
  scale_color_manual(values=cols2,labels=c("ME_P1","ME_P2"))

viral_reads<-reads_plot %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=log(viral_reads),colour=Background))+
  geom_point()+
  geom_smooth(aes(group=Background),se=FALSE)+
  theme_bw()+
  facet_grid(~type,labeller=as_labeller(facet_names))+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ggtitle("Total Viral Reads")+
  ylab("Log(Total Viral Reads)")+
  scale_color_manual(values=cols2,labels=c("ME_P1","ME_P2"))

ATCC_reads<-reads_plot %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=log(ATCC_reads),colour=Background))+
  geom_point()+
  geom_smooth(aes(group=Background),se=FALSE)+
  theme_bw()+
  facet_grid(~type,labeller=as_labeller(facet_names))+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ggtitle("Spike In Viral Reads")+
  ylab("Log(Spike In Viral Reads)")+
  scale_color_manual(values=cols2,labels=c("ME_P1","ME_P2"))

prop_ATCC<-reads_plot %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=prop_ATCC,colour=Background))+
  geom_point()+
  geom_smooth(aes(group=Background),se=FALSE)+
  facet_grid(~type,labeller=as_labeller(facet_names))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ggtitle("Proportion Spike In Viral")+
  ylab("Spike In Viral Reads/Total Reads")+
  scale_color_manual(values=cols2,labels=c("ME_P1","ME_P2"))

prop_viral<-reads_plot %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=prop_viral,colour=Background))+
  geom_point()+
  geom_smooth(aes(group=Background),se=FALSE)+
  theme_bw()+
  facet_grid(~type,labeller=as_labeller(facet_names))+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ggtitle("Proportion Viral")+
  ylab("Viral Reads/Total Reads")+
  scale_color_manual(values=cols2,labels=c("ME_P1","ME_P2"))

#ggarrange(total_reads,viral_reads,ATCC_reads,prop_viral,prop_ATCC,nrow=2,ncol=3,common.legend = TRUE)

ggarrange(total_reads,viral_reads,prop_viral,nrow=3,common.legend = TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/backgrounds_reads.png")

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/FigureS3.pdf")

####compare proportion viral reads per pool with shotgun####

metadata3 <- metadata %>%
  select(Sample.ID,Number.of.read.pairs..quality.adaptor.trimmed.,Pool.for.sequencing) %>%
  rename(Sample_id = Sample.ID) %>%
  rename(QC_reads = Number.of.read.pairs..quality.adaptor.trimmed.)

reads_viral_dedup3<-full_join(metadata3,viral_reads_dedup,by="Sample_id") 

reads_viral_nodedup3<-full_join(metadata3,viral_reads_nodedup,by="Sample_id") 

reads_viral_dedup3 %>%
  group_by(Pool.for.sequencing) %>%
  summarise(viral_reads = sum(total_virus_reads),
            total_reads = sum(QC_reads*2),
            prop_viral = viral_reads/total_reads)

polyomics<-read.csv("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_polyomics/total_virus_mapped_reads_per_sample_dedup.csv")

#read data from polyomics_indexes_stefano document
total_reads<-data.frame(total_reads = c(5942760,7714275))

keeps<-c("RNA-Msp-p2","RNA-Msp-p8")

polyomics_read_prop<-polyomics %>%
  filter(Sample_id %in% keeps)

polyomics_read_prop2<-cbind(polyomics_read_prop,total_reads)

polyomics_read_prop2 %>%
summarise(prop_viral = total_virus_reads/total_reads)
