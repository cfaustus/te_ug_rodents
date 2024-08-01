#preliminary plots of TE experiment data
#July 2024

library(tidyverse)
library(ggpubr)

#deduplicated plot - coverage, depth, reads

#genome coverage

cov_dedup<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_genomeCoverage_spike-in_virus_dedup.tsv", sep = "\t", header = TRUE)

p1<-ggplot(cov_dedup,aes(x=virus,y=genome_coverage,colour = Background))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  facet_grid(~as.character(Viral.load))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_text(size=10))+
  ylab("Genome coverage")

#read depths

depth_dedup<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_sample_and_virus_dedup.tsv", sep = "\t", header = TRUE)

p2<-ggplot(depth_dedup,aes(x=virus,y=log(mean_depth),colour = Background))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  facet_grid(~as.character(Viral.load))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_text(size=10))+
  ylab("Depth(log(mean depth))")


#read counts

counts_dedup<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readcount_per_sample_and_virus_dedup.tsv", sep = "\t", header = TRUE)

p3<-ggplot(counts_dedup,aes(x=virus,y=log(matched),colour = Background))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  facet_grid(~as.character(Viral.load))+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  ylab("Read Count(log(count))")

ggarrange(p1,p2,p3,nrow=3,common.legend=TRUE,legend="right",heights=c(1,1,1.5))

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/dedup_viruses_sep.pdf")

#summarised across viruses

p1.1<-ggplot(cov_dedup,aes(x=as.character(Viral.load),y=genome_coverage))+
  geom_point()+
  geom_boxplot(outlier.shape=NA)+
  ylab("Genome coverage")+
  xlab("Viral load")

p2.1<-ggplot(depth_dedup,aes(x=as.character(Viral.load),y=log(mean_depth)))+
  geom_point()+
  geom_boxplot(outlier.shape=NA)+
  ylab("Depth(log(mean depth))")+
  xlab("Viral load")


p3.1<-ggplot(counts_dedup,aes(x=as.character(Viral.load),y=log(matched)))+
  geom_point()+
  geom_boxplot(outlier.shape=NA)+
  ylab("Read Count(log(count))")+
  xlab("Viral load")


ggarrange(p1.1,p2.1,p3.1,ncol=3)

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/dedup_viruses_comb.pdf")

#not deduplicated plot - coverage, depth, reads

#genome coverage

cov_no_dedup<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_genomeCoverage_spike-in_virus_nodedup.tsv", sep = "\t", header = TRUE)

p4<-ggplot(cov_no_dedup,aes(x=virus,y=genome_coverage,colour = Background))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  facet_grid(~as.character(Viral.load))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_text(size=10))+
  ylab("Genome coverage")

#read depths

depth_no_dedup<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_sample_and_virus_nodedup.tsv", sep = "\t", header = TRUE)

p5<-ggplot(depth_no_dedup,aes(x=virus,y=log(mean_depth),colour = Background))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  facet_grid(~as.character(Viral.load))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_text(size=10))+
  ylab(" Depth(log(mean depth))")


#read counts

counts_no_dedup<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readcount_per_sample_and_virus_nodedup.tsv", sep = "\t", header = TRUE)

p6<-ggplot(counts_no_dedup,aes(x=virus,y=log(matched),colour = Background))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  facet_grid(~as.character(Viral.load))+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  ylab("Read Count(log(count))")

ggarrange(p4,p5,p6,nrow=3,common.legend=TRUE,legend="right",heights=c(1,1,1.5))

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/non_dedup_viruses_sep.pdf")

#summarised across viruses

p4.1<-ggplot(cov_no_dedup,aes(x=as.character(Viral.load),y=genome_coverage))+
  geom_point()+
  geom_boxplot()+
  ylab("Genome coverage")+
  xlab("Viral load")

p5.1<-ggplot(depth_no_dedup,aes(x=as.character(Viral.load),y=log(mean_depth)))+
  geom_point()+
  geom_boxplot()+
  ylab("Depth(log(mean depth))")+
  xlab("Viral load")

p6.1<-ggplot(counts_no_dedup,aes(x=as.character(Viral.load),y=log(matched)))+
  geom_point()+
  geom_boxplot()+
  ylab("Read Count(log(count))")+
  xlab("Viral load")

ggarrange(p4.1,p5.1,p6.1,ncol=3)

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/non_dedup_viruses_comb.pdf")

#compare with or without deduplication

cov_both<-rbind(cov_dedup,cov_no_dedup)

p7<-ggplot(cov_both,aes(x=as.character(Viral.load),y=genome_coverage,color=dataset))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  ylab("Genome coverage")+
  xlab("Viral load")

depth_both<-rbind(depth_dedup,depth_no_dedup)

p8<-ggplot(depth_both,aes(x=as.character(Viral.load),y=log(mean_depth),color=dataset))+geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  ylab("Depth(log(mean depth))")+
  xlab("Viral load")

counts_both<-rbind(counts_dedup,counts_no_dedup)

p9<-ggplot(counts_both,aes(x=as.character(Viral.load),y=log(matched),color=dataset))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  ylab("Read Count(log(count))")+
  xlab("Viral load")

ggarrange(p7,p8,p9,ncol=3,common.legend=TRUE)

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_dedup.pdf")

#lower coverage of Zika virus across read depths, especially in deduplicated dataset

ggplot(cov_both,aes(x=as.character(Viral.load),y=genome_coverage,colour=dataset))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  facet_grid(~virus)+
  ylab("Genome coverage")+
  xlab("Viral load")

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/coverage_per_virus.pdf")

#library capture pool

metadata<-read.csv("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/metadata/sampleIDs_TESpikeIn.csv",header=TRUE)

metadata_sub<- metadata %>%
  select("Sample.ID","Pool.for.sequencing") %>%
  rename(Sample_id = Sample.ID)

cov_pools<-left_join(cov_dedup,metadata_sub,by="Sample_id")

p10<-ggplot(cov_pools,aes(x=virus,y=genome_coverage,colour = Pool.for.sequencing))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Genome coverage")

depth_pools<-left_join(depth_dedup,metadata_sub,by="Sample_id")

p11<-ggplot(depth_pools,aes(x=virus,y=log(mean_depth),colour = Pool.for.sequencing))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Depth(log(mean depth))")

counts_pools<-left_join(counts_dedup,metadata_sub,by="Sample_id")

p12<-ggplot(counts_pools,aes(x=virus,y=log(matched),colour = Pool.for.sequencing))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1))+
  ylab("Read Count(log(count))")

ggarrange(p10,p11,p12,ncol=3,common.legend=TRUE)

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_pools_hyb.pdf")

#combined across viruses

p10.1<-ggplot(cov_pools,aes(x=Pool.for.sequencing,y=genome_coverage))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  ylab("Genome coverage")+
  xlab("Hybridization pool")

p11.1<-ggplot(depth_pools,aes(x=Pool.for.sequencing,y=log(mean_depth)))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  ylab("Depth(log(mean depth))")+
  xlab("Hybridization pool")

p12.1<-ggplot(counts_pools,aes(x=Pool.for.sequencing,y=log(matched)))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  ylab("Read Count(log(count))")+
  xlab("Hybridization pool")

ggarrange(p10.1,p11.1,p12.1,ncol=3)

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_pools_hyb_combined.pdf")
