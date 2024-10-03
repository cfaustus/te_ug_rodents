#plots of TE experiment data - ATCC genomes updated version
#July 2024

library(tidyverse)
library(ggpubr)

#metadata

metadata<-read.csv("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/metadata/sampleIDs_TESpikeIn.csv",header=TRUE)

metadata2 <- metadata %>%
  select(Sample.ID,Number.of.read.pairs..quality.adaptor.trimmed.) %>%
  rename(Sample_id = Sample.ID) %>%
  rename(QC_reads = Number.of.read.pairs..quality.adaptor.trimmed.)

#import and combine read count files
#deduplicated

counts_dedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readcount_per_sample_and_virus_bowtie2_dedup_atcc_ref.tsv", sep = "\t", header = TRUE)

counts_dedup_ngm<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readcount_per_sample_and_virus_ngm_dedup_atcc_ref.tsv", sep = "\t", header = TRUE)

counts_dedup_all<-rbind(counts_dedup_bt,counts_dedup_ngm)

counts_reads<-left_join(counts_dedup_all,metadata2,by="Sample_id")

counts_reads_norm<-counts_reads %>%
  mutate(norm_counts = matched/QC_reads)

p1<-counts_reads_norm %>%
  filter(mapper == "ngm") %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
ggplot(aes(x=virus,y=log(norm_counts),colour = Background))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  facet_grid(~as.character(Viral.load))+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  ylab("Normalised Read Count(log(count))")

p2<-counts_reads_norm %>%
  filter(mapper == "ngm") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=as.character(Viral.load),y=log(norm_counts)))+
  geom_point()+
  geom_boxplot(outlier.shape=NA)+
  ylab("Normalised Read Count(log(count))")+
  xlab("Viral load")

#import and combine read depth files - deduplicated

depth_dedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_sample_and_virus_bowtie2_dedup_atcc_ref.tsv", sep = "\t", header = TRUE)

depth_dedup_ngm<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_sample_and_virus_ngm_dedup_atcc_ref.tsv", sep = "\t", header = TRUE)

depths_dedup_all<-rbind(depth_dedup_bt,depth_dedup_ngm)

depths_reads<-left_join(depths_dedup_all,metadata2,by="Sample_id")

depths_norm<-depths_reads %>%
  mutate(mean_depth_norm = mean_depth/QC_reads) %>%
  mutate(median_depth_norm = median_depth/QC_reads)

p3<-depths_norm %>%
  filter(mapper == "ngm") %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%  
ggplot(aes(x=virus,y=log(mean_depth_norm),colour = Background))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  facet_grid(~as.character(Viral.load))+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  ylab("Normalised Depth(log(mean depth))")

p4<-depths_norm %>%
  filter(mapper == "ngm") %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%  
  ggplot(aes(x=virus,y=log(median_depth_norm),colour = Background))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  facet_grid(~as.character(Viral.load))+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  ylab("Normalised Depth(log(median depth))")

p5<-depths_norm %>%
  filter(mapper == "ngm") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=as.character(Viral.load),y=log(mean_depth_norm)))+
  geom_point()+
  geom_boxplot(outlier.shape=NA)+
  ylab("Normalised Depth(log(mean depth))")+
  xlab("Viral load")

p6<-depths_norm %>%
  filter(mapper == "ngm") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=as.character(Viral.load),y=log(median_depth_norm)))+
  geom_point()+
  geom_boxplot(outlier.shape=NA)+
  ylab("Normalised Depth(log(median depth))")+
  xlab("Viral load")

ggarrange(p2,p6)

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/norm_read_count_depth.pdf")

#import and combine genome coverage files - deduplicated

cov_dedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_site_sample_and_virus_bowtie2_dedup_atcc_ref.tsv", sep = "\t", header = TRUE)

cov_dedup_ngm<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_site_sample_and_virus_ngm_dedup_atcc_ref.tsv", sep = "\t", header = TRUE)


##compare deduplicated and non-deduplicated

#import and combine read count files
#non-deduplicated

counts_nondedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readcount_per_sample_and_virus_bowtie2_nodedup_atcc_ref.tsv", sep = "\t", header = TRUE)

counts_nondedup_ngm<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readcount_per_sample_and_virus_ngm_nodedup_atcc_ref.tsv", sep = "\t", header = TRUE)

counts_nondedup_all<-rbind(counts_nondedup_bt,counts_nondedup_ngm)

counts_all<-rbind(counts_dedup_all,counts_nondedup_all)

p7<-counts_all %>%
  filter(mapper == "ngm") %>%
  filter(Background != "p6") %>%
ggplot(aes(x=as.character(Viral.load),y=log(matched),color=type))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  ylab("Read Count(log(count))")+
  xlab("Viral load")

#import and combine read depth files
#non-deduplicated

depth_nondedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_sample_and_virus_bowtie2_nodedup_atcc_ref.tsv", sep = "\t", header = TRUE)

depth_nondedup_ngm<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_sample_and_virus_ngm_nodedup_atcc_ref.tsv", sep = "\t", header = TRUE)

depth_nondedup_all<-rbind(depth_nondedup_bt,depth_nondedup_ngm)

depth_all<-rbind(depths_dedup_all,depth_nondedup_all)

p8<-depth_all %>%
  filter(mapper == "ngm") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=as.character(Viral.load),y=log(mean_depth),color=type))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  ylab("Read Depth(log(mean depth))")+
  xlab("Viral load")

ggarrange(p7,p8,common.legend = TRUE)

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/dedup_compare.pdf")

#compare library capture pool

metadata3 <- metadata %>%
  select(Sample.ID,Pool.for.sequencing) %>%
  rename(Sample_id = Sample.ID) %>%
  rename(Pool = Pool.for.sequencing)

counts_pool<-left_join(counts_dedup_all,metadata3,by="Sample_id")

p9<-counts_pool %>%
  filter(mapper == "ngm") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=virus,y=log(matched),colour = Pool))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Read Count (log(count)")

p10<-counts_pool %>%
  filter(mapper == "ngm") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=Pool,y=log(matched)))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  ylab("Read Count(log(count))")+
  xlab("Hybridization pool")

depths_pool<-left_join(depths_dedup_all,metadata3,by="Sample_id")

p11<-depths_pool %>%
  filter(mapper == "ngm") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=virus,y=log(median_depth),colour = Pool))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Read Depth (log(median depth)")

p12<-depths_pool %>%
  filter(mapper == "ngm") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=Pool,y=log(median_depth)))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  ylab("Read Depth (log(median depth)")
  xlab("Hybridization pool")

ggarrange(p10,p12)  

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/pool_compare.pdf")
