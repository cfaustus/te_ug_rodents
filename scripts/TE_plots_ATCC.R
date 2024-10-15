#plots of TE experiment data - ATCC genomes updated version
#October 2024

#viral load vs read count normalised using different methods, comparing deduplicated and non-deduplicated
#use bowtie2 data

rm(list=ls())

library(tidyverse)
library(ggpubr)

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

#normalise by both raw read count and genome length - should be the same as mean read depth

counts_reads_norm<-counts_reads %>%
  mutate(norm_counts1 = matched/QC_reads) %>%
  mutate(norm_counts2 = matched/QC_reads/length) %>%
  mutate(norm_counts3 = matched/length)
 
#plot viral read counts (non normalised - log and not log scale)

read_count_raw<-counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=matched,colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus),se=FALSE)+
  facet_grid(~type)+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ylab("Viral Reads")

read_count<- counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
ggplot(aes(x=Viral.load,y=log(matched),colour = virus))+
  geom_point()+
   geom_smooth(aes(group=virus),se=FALSE)+
   facet_grid(~type)+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ylab("Log(Viral Reads)")

#normalise by raw read count

read_count_norm1<- counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=log(norm_counts1),colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus),se=FALSE)+
  facet_grid(~type)+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ylab("Log(viral reads/cleaned reads)")

raw_read_count_norm1<- counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=norm_counts1,colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus),se=FALSE)+
  facet_grid(~type)+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ylab("Viral reads/cleaned reads")

#normalise by raw read count and genome length

read_count_norm2<- counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=log(norm_counts2),colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus),se=FALSE)+
  facet_grid(~type)+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ylab("Log(viral reads/cleaned reads/genome length)")

raw_read_count_norm2<- counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=norm_counts2,colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus),se=FALSE)+
  facet_grid(~type)+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ylab("Viral reads/cleaned reads/genome length")

read_count_norm3<- counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=log(norm_counts3),colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus),se=FALSE)+
  facet_grid(~type)+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ylab("Log(viral reads/genome length)")

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/norm_genome_length.pdf")

####read depths####

#import and combine read depth files

depth_dedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_sample_and_virus_bowtie2_dedup_atcc_ref.tsv", sep = "\t", header = TRUE)

depth_nodedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_sample_and_virus_bowtie2_nodedup_atcc_ref.tsv", sep = "\t", header = TRUE)

depths_bt_all<-rbind(depth_dedup_bt,depth_nodedup_bt)

depths_reads<-left_join(depths_bt_all,metadata2,by="Sample_id")

read_depths<- depths_reads %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=log(mean_depth),colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus),se=FALSE)+
  facet_grid(~type)+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ylab("Log(mean read depth)")

raw_read_depths<- depths_reads %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=mean_depth,colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus),se=FALSE)+
  facet_grid(~type)+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ylab("Mean read depth")

ggarrange(read_count,read_count_norm1,read_count_norm2,read_depths,nrow=2,ncol=2,common.legend=TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/log_read_count_depth.pdf")

ggarrange(read_count_raw,raw_read_count_norm1,raw_read_count_norm2,raw_read_depths,nrow=2,ncol=2,common.legend=TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/read_count_depth.pdf")

#apparently normalised read count using method 2 should not actually be the same

####read counts/depths split by background####

#break down by Background x virus

background_counts_reads<-counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=log(matched),colour=virus))+
  geom_point()+
  geom_smooth(aes(group=virus),se=FALSE)+
  facet_grid(Background~type)+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ylab("Log(Viral Reads)")

backgrounds_counts_reads_norm<-counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=log(norm_counts1),colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus),se=FALSE)+
  facet_grid(Background~type)+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ylab("Log(viral reads/cleaned reads)")

backgrounds_counts_reads_norm2<-counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=log(norm_counts2),colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus),se=FALSE)+
  facet_grid(Background~type)+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ylab("Log(viral reads/cleaned reads/genome length)")

background_read_depths<- depths_reads %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=log(mean_depth),colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus),se=FALSE)+
  facet_grid(Background~type)+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ylab("Log(mean read depth)")

ggarrange(background_counts_reads,backgrounds_counts_reads_norm,backgrounds_counts_reads_norm2,background_read_depths,nrow=2,ncol=2,common.legend=TRUE)

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/backgrounds_read_count_depth.pdf")

#there is a difference between backgrounds in where plateauing occurs - look at whether total or viral read counts can explain this discrepancy?

total_reads_backgrounds<-counts_reads_norm %>%
  ggplot(aes(x=Viral.load,y=log(QC_reads),colour=Background))+
  geom_point()+
  geom_smooth(aes(group=Background),se=FALSE)+
  facet_grid(~type)+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ggtitle("Total Reads")+
  ylab("Log(Total Reads)")

viral_reads_backgrounds<-counts_reads_norm %>%
  group_by(Background,Viral.load,type) %>%
  summarise(viral_reads = sum(matched)) %>%
  ggplot(aes(x=Viral.load,y=log(viral_reads),colour=Background))+
  geom_point()+
  geom_smooth(aes(group=Background),se=FALSE)+
  facet_grid(~type)+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ggtitle("Viral Reads")+
  ylab("Log(Viral Reads)")

ggarrange(total_reads_backgrounds,viral_reads_backgrounds,nrow=2,common.legend = TRUE)
  
ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/backgrounds_reads.pdf")

####compare library capture pool####

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
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_text(size=10))+
  ylab("Log(Viral Reads)")

pool_counts_sum<-counts_pool %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Pool,y=log(matched)))+
  geom_boxplot()+
  facet_grid(~type)+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  ylab("Log(Viral Reads)")

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
  ggplot(aes(x=Pool,y=log(norm_counts1)))+
  geom_boxplot()+
  facet_grid(~type)+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  ylab("Log(viral reads/cleaned reads)")

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
  ggplot(aes(x=Pool,y=log(norm_counts2)))+
  geom_boxplot()+
  facet_grid(~type)+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  ylab("Log(viral reads/cleaned reads/genome length)")

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
  ggplot(aes(x=Pool,y=log(mean_depth)))+
  geom_boxplot()+
  facet_grid(~type)+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  ylab("Log(mean read depth)")

ggarrange(pool_readcounts,pool_readcounts_norm,pool_readcounts_norm2,pool_readdepths,nrow=2,ncol=2,heights=c(1,2),common.legend=TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/pool_compare_viruses.pdf")

ggarrange(pool_counts_sum,pool_norm1_sum,pool_norm2_sum,pool_depths_sum,nrow=2,ncol=2,common.legend = TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/pool_compare_viruses_sum.pdf")

####split by viruses####

counts_split<-counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
ggplot(aes(x=virus,y=log(matched),colour=type))+
  geom_boxplot()+
  facet_grid(~as.character(Viral.load))+
  ylab("Log(Viral Reads)")+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_text(size=10))

counts_norm_split<-counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=virus,y=log(norm_counts1),colour=type))+
  geom_boxplot()+
  facet_grid(~as.character(Viral.load))+
  ylab("Log(viral reads/cleaned reads)")+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_text(size=10))

counts_norm2_split<-counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=virus,y=log(norm_counts2),colour=type))+
  geom_boxplot()+
  facet_grid(~as.character(Viral.load))+
  ylab("Log(viral reads/cleaned reads/genome length)")+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))

depth_split<-depths_pool %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=virus,y=log(mean_depth),colour=type))+
  geom_boxplot()+
  facet_grid(~as.character(Viral.load))+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  ylab("Log(mean read depth)")+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))

ggarrange(counts_split,counts_norm_split,counts_norm2_split,depth_split,nrow=2,ncol=2,heights=c(1,2),common.legend = TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/compare_viruses_split.pdf")

####genome coverage####

dedup_per_site<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_site_sample_and_virus_bowtie2_dedup_atcc_ref.tsv", sep = "\t", header = TRUE)

nodedup_per_site<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_site_sample_and_virus_bowtie2_nodedup_atcc_ref.tsv", sep = "\t", header = TRUE)

#add length column from read depth file to per site file

lengths<-counts_reads_norm %>%
  select(virus,length) %>%
  distinct()

#calculate genome coverage

persite_coverage_dedup<-dedup_per_site %>%
  group_by(Sample_id,virus,Viral.load,Background,type) %>%
  filter(coverage > 0) %>%
  summarise(genome_coverage = n()) %>%
  left_join(lengths) %>%
  mutate(percent_coverage = genome_coverage/length)

persite_coverage_nodedup<-nodedup_per_site %>%
  group_by(Sample_id,virus,Viral.load,Background,type) %>%
  filter(coverage > 0) %>%
  summarise(genome_coverage = n()) %>%
  left_join(lengths) %>%
  mutate(percent_coverage = genome_coverage/length)

persite_coverage_both<-rbind(persite_coverage_dedup,persite_coverage_nodedup)

ggplot(persite_coverage,aes(x=virus,y=percent_coverage,colour = Background))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  facet_grid(type~as.character(Viral.load))+
  theme(axis.title.y=element_text(size=10))+
  ylab("Genome coverage")+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/genome_cov_split.pdf")

##are the deduplicated and non deduplicated datasets identical??

ggplot(persite_coverage,aes(x=as.character(Viral.load),y=percent_coverage))+
  geom_point()+
  geom_boxplot(outlier.shape=NA)+
  facet_wrap(~type)+
  ylab("Genome coverage")+
  xlab("Viral load")

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/genome_cov_combined.pdf")

persite_both<-rbind(dedup_per_site,nodedup_per_site)

persite_reads<-left_join(persite_both,metadata2,by="Sample_id")

persite_norm<-persite_reads %>%
    mutate(norm_cov = coverage/QC_reads)

RSV<-persite_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(virus == "Human_respiratory_syncytial_virus") %>%
  ggplot(aes(x=site,y=norm_cov,fill=as.character(Viral.load)))+
  geom_col(position="dodge")+
  facet_wrap(type~Background)+
  ylab("Normalised coverage")+
  ggtitle("Human RSV")+
  guides(fill=guide_legend(title="Viral load"))

ZIKV<-persite_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(virus == "Zika_virus") %>%
  ggplot(aes(x=site,y=norm_cov,fill=as.character(Viral.load)))+
  geom_col(position="dodge")+
  facet_wrap(type~Background)+
  ylab("Normalised coverage")+
  ggtitle("Zika virus")+
  guides(fill=guide_legend(title="Viral load"))

HADV<-persite_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(virus == "Human_adenovirus_40") %>%
  ggplot(aes(x=site,y=norm_cov,fill=as.character(Viral.load)))+
  geom_col(position="dodge")+
  facet_grid(type~Background)+
  ylab("Normalised coverage")+
  ggtitle("Human adenovirus")+
  guides(fill=guide_legend(title="Viral load"))

HBHV<-persite_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(virus == "Human_betaherpesvirus") %>%
  ggplot(aes(x=site,y=norm_cov,fill=as.character(Viral.load)))+
  geom_col(position="dodge")+
  facet_grid(type~Background)+
  ylab("Normalised coverage")+
  ggtitle("Human betaherpesvirus")+
  guides(fill=guide_legend(title="Viral load"))

#have to do these differently due to segmentation
FLUB_dedup<-persite_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(virus == "Influenza_B_virus") %>%
  filter(type == "dedup_TE") %>%
  ggplot(aes(x=site,y=norm_cov,fill=as.character(Viral.load)))+
  geom_col(position="dodge")+
  facet_grid(seg~Background)+
  ylab("Normalised coverage")+
  ggtitle("Influenza B deduplicated")+
  guides(fill=guide_legend(title="Viral load"))

FLUB_nodedup<-persite_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(virus == "Influenza_B_virus") %>%
  filter(type == "nodedup_TE") %>%
  ggplot(aes(x=site,y=norm_cov,fill=as.character(Viral.load)))+
  geom_col(position="dodge")+
  facet_grid(seg~Background)+
  ylab("Normalised coverage")+
  ggtitle("Influenza B non deduplicated")+
  guides(fill=guide_legend(title="Viral load"))

ggarrange(FLUB_dedup,FLUB_nodedup,ncol=2,common.legend=TRUE)

REOV_dedup<-persite_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(virus == "Mammalian_orthoreovirus3") %>%
  filter(type == "dedup_TE") %>%
  ggplot(aes(x=site,y=norm_cov,fill=as.character(Viral.load)))+
  geom_col(position="dodge")+
  facet_grid(seg~Background)+
  ylab("Normalised coverage")+
  ggtitle("Orthoreovirus deduplicated")+
  guides(fill=guide_legend(title="Viral load"))

REOV_nodedup<-persite_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(virus == "Mammalian_orthoreovirus3") %>%
  filter(type == "nodedup_TE") %>%
  ggplot(aes(x=site,y=norm_cov,fill=as.character(Viral.load)))+
  geom_col(position="dodge")+
  facet_grid(seg~Background)+
  ylab("Normalised coverage")+
  ggtitle("Orthoreovirus non deduplicated")+
  guides(fill=guide_legend(title="Viral load"))

ggarrange(REOV_dedup,REOV_nodedup,ncol=2,common.legend=TRUE)
