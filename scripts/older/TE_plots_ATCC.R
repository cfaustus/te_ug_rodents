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
  mutate(norm_counts3 = matched/length) %>%
  mutate(genome_structure = case_when((virus == "Human_adenovirus_40"| virus == "Human_betaherpesvirus") ~ "DNA",
                                     ,.default = "RNA"))

#change labels in facet plots

facet_names<-c("dedup_TE" = "Deduplicated",
               "nodedup_TE"= "Non-Deduplicated")

cols<-c("#4477AA","#66CCEE","#228833","#CCBB44","#EE6677","#AA3377")
 
#plot viral read counts (non normalised - log and not log scale)

read_count_raw<-counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=matched,colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus),se=FALSE)+
  facet_grid(~type,labeller=as_labeller(facet_names))+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  ylab("Viral Reads")

read_count<- counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
ggplot(aes(x=Viral.load,y=log(matched),colour = virus))+
  geom_point()+
   geom_smooth(aes(group=virus,linetype=genome_structure),se=FALSE)+
   facet_grid(~type,labeller=as_labeller(facet_names))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10),legend.title=element_blank(),legend.position="top")+
  scale_x_log10()+
  scale_color_manual(values=cols,labels=c("Human adenovirus 40","Human betaherpesvirus","Human respiratory syncytial virus","Influenza B virus","Mammalian orthoreovirus 3","Zika virus"))+
  ylab("Log(Viral Reads)")

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/read_counts.pdf",width=8,height=6)

#normalise by raw read count

read_count_norm1<- counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=log(norm_counts1),colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus,linetype=genome_structure),se=FALSE)+
  facet_grid(~type,labeller=as_labeller(facet_names))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10),legend.title=element_blank())+
  scale_x_log10()+
  scale_color_manual(values=cols,labels=c("Human adenovirus 40","Human betaherpesvirus","Human respiratory syncytial virus","Influenza B virus","Mammalian orthoreovirus 3","Zika virus"))+
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
  geom_smooth(aes(group=virus,linetype=genome_structure),se=FALSE)+
  facet_grid(~type,labeller=as_labeller(facet_names))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10),legend.title=element_blank())+
  scale_x_log10()+
  scale_color_manual(values=cols,labels=c("Human adenovirus 40","Human betaherpesvirus","Human respiratory syncytial virus","Influenza B virus","Mammalian orthoreovirus 3","Zika virus"))+
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
  facet_grid(~type,labeller=as_labeller(facet_names))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10))+
  scale_x_log10()+
  scale_color_hue(labels=c("Human adenovirus 40","Human betaherpesvirus","Human respiratory syncytial virus","Influenza B virus","Mammalian orthoreovirus 3","Zika virus"))+
  ylab("Log(viral reads/genome length)")

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/norm_genome_length.pdf")

#plot just the two lefthand panels from the original figure

ggarrange(read_count,read_count_norm2,nrow=2,common.legend = TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/log_read_count_depth_2panels.png",width=10,height=7)

####read depths####

#import and combine read depth files

depth_dedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_sample_and_virus_bowtie2_dedup_atcc_ref.tsv", sep = "\t", header = TRUE)

depth_nodedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_sample_and_virus_bowtie2_nodedup_atcc_ref.tsv", sep = "\t", header = TRUE)

depths_bt_all<-rbind(depth_dedup_bt,depth_nodedup_bt)

depths_reads<-left_join(depths_bt_all,metadata2,by="Sample_id")%>%
  mutate(genome_structure = case_when((virus == "Human_adenovirus_40"| virus == "Human_betaherpesvirus") ~ "DNA",
                                      ,.default = "RNA"))

read_depths<- depths_reads %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=log(mean_depth),colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus,linetype=genome_structure),se=FALSE)+
  facet_grid(~type,labeller=as_labeller(facet_names))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10),legend.title=element_blank(),legend.position="top")+
  scale_color_manual(values=cols,labels=c("Human adenovirus 40","Human betaherpesvirus","Human respiratory syncytial virus","Influenza B virus","Mammalian orthoreovirus 3","Zika virus"))+
  scale_x_log10()+
  ylab("Log(mean read depth)")

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/mean_depth.pdf",width=8,height=6)

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

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/log_read_count_depth.png",width=10,height=7)

ggarrange(read_count_raw,raw_read_count_norm1,raw_read_count_norm2,raw_read_depths,nrow=2,ncol=2,common.legend=TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/read_count_depth.pdf")

#apparently normalised read count using method 2 should not actually be the same

####read counts/depths split by background####

#change labels in facet plots

facet_names_bg<-c("dedup_TE" = "Deduplicated",
               "nodedup_TE"= "Non-Deduplicated",
               "p2" = "ME_P1",
               "p8" = "ME_P2")

#break down by Background x virus

background_counts_reads<-counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=log(matched),colour=virus))+
  geom_point()+
  geom_smooth(aes(group=virus,linetype=genome_structure),se=FALSE)+
  facet_grid(Background~type,labeller=as_labeller(facet_names_bg))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10),legend.title=element_blank())+
  scale_color_manual(values=cols,labels=c("Human adenovirus 40","Human betaherpesvirus","Human respiratory syncytial virus","Influenza B virus","Mammalian orthoreovirus 3","Zika virus"))+
  scale_x_log10()+
  ylab("Log(Viral Reads)")

backgrounds_counts_reads_norm<-counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=log(norm_counts1),colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus,linetype=genome_structure),se=FALSE)+
  facet_grid(Background~type,labeller=as_labeller(facet_names_bg))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10),legend.title=element_blank())+
  scale_color_manual(values=cols,labels=c("Human adenovirus 40","Human betaherpesvirus","Human respiratory syncytial virus","Influenza B virus","Mammalian orthoreovirus 3","Zika virus"))+
  scale_x_log10()+
  ylab("Log(viral reads/cleaned reads)")

backgrounds_counts_reads_norm2<-counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=log(norm_counts2),colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus,linetype=genome_structure),se=FALSE)+
  facet_grid(Background~type,labeller=as_labeller(facet_names_bg))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10),legend.title=element_blank())+
  scale_color_manual(values=cols,labels=c("Human adenovirus 40","Human betaherpesvirus","Human respiratory syncytial virus","Influenza B virus","Mammalian orthoreovirus 3","Zika virus"))+
  scale_x_log10()+
  ylab("Log(viral reads/cleaned reads/genome length)")

background_read_depths<- depths_reads %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=Viral.load,y=log(mean_depth),colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus,linetype=genome_structure),se=FALSE)+
  facet_grid(Background~type,labeller=as_labeller(facet_names_bg))+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10),legend.title=element_blank())+
  scale_color_manual(values=cols,labels=c("Human adenovirus 40","Human betaherpesvirus","Human respiratory syncytial virus","Influenza B virus","Mammalian orthoreovirus 3","Zika virus"))+
  scale_x_log10()+
  ylab("Log(mean read depth)")

ggarrange(background_counts_reads,backgrounds_counts_reads_norm,backgrounds_counts_reads_norm2,background_read_depths,nrow=2,ncol=2,common.legend=TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/backgrounds_read_count_depth.png",width=10,height=7)

####model for spike in viruses - mean read depth####
#what kind of data is read depth, is it more like a count or a proportion?

library(lme4)
library(boot)
library(car)

depths_reads_sub<-depths_reads %>%
filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(type == "dedup_TE")

lm1<-glm(mean_depth~Viral.load+virus,data=depths_reads_sub)

summary(lm1)

glm.diag.plots(lm1)

predict.glm(lm1)

hist(depths_reads_sub$mean_depth)

#residuals plot shows funnel shape - homogeneity of variance not met. need to check how to resolve that

lm2<-lmer(mean_depth~Viral.load+ (1|virus),data=depths_reads_sub)

summary(lm3)
Anova(lm3)

#could also try with read count which is more straightforward to model potentially

counts_reads_sub<-counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(type == "dedup_TE")

hist(counts_reads_sub$matched)

lm4<-glm(matched ~ Viral.load+virus,data=counts_reads_sub)

summary(lm4)

lm5<-glm(matched ~ Viral.load+virus,data=counts_reads_sub,family="Gamma")

summary(lm5)

####compare library capture pool####

cols3<-c("#228833","#AA3377")

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

coverage_labels<-c("0" = "Control",
                   "100" = "1e+02",
                   "1000" = "1e+03",
                   "100000" = "1e+05",
    "dedup_TE" = "Deduplicated",
                   "nodedup_TE"= "Non-Deduplicated")

coverage_labels2<-c("0" = "Control",
                   "100" = "1e+02",
                   "1000" = "1e+03",
                   "100000" = "1e+05")

persite_coverage_both %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
ggplot(aes(x=virus,y=percent_coverage,colour = Background))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  facet_grid(type~as.character(Viral.load),labeller=as_labeller(coverage_labels))+
  theme_bw()+
  theme(axis.title.y=element_text(size=10))+
  ylab("Genome coverage")+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=15),axis.text=element_text(size=12),strip.text=element_text(size=12))+
  scale_color_manual(values=cols2,labels=c("ME_P1","ME_P2"))+
  scale_x_discrete(labels=c("Human adenovirus 40","Human betaherpesvirus","Human respiratory syncytial virus","Influenza B virus","Mammalian orthoreovirus 3","Zika virus"))+
  scale_y_continuous(limits=c(0,1))

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/genome_cov_split.png",width=15,height=6)

##are the deduplicated and non deduplicated datasets identical??
##if they are the same can maybe just show one
#show deduplicated

persite_coverage_both %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(type == "dedup_TE") %>%
  ggplot(aes(x=virus,y=percent_coverage,colour = Background))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  facet_grid(~as.character(Viral.load),labeller=as_labeller(coverage_labels2))+
  theme_bw()+
  ylab("Genome coverage")+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=15),axis.text=element_text(size=12),strip.text=element_text(size=12),legend.text = element_text(size=10))+
  scale_color_manual(values=cols2,labels=c("ME_P1","ME_P2"))+
  scale_x_discrete(labels=c("Human adenovirus 40","Human betaherpesvirus","Human respiratory syncytial virus","Influenza B virus","Mammalian orthoreovirus 3","Zika virus"))+
  scale_y_continuous(limits=c(0,1))

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/genome_cov_deduponly.png",width=15,height=6)

ggplot(persite_coverage_both,aes(x=as.character(Viral.load),y=percent_coverage))+
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

coverage_labels3<-c("100" = "1e+02",
                   "1000" = "1e+03",
                   "100000" = "1e+05",
                   "p2" = "ME_P1",
                   "p8"= "ME_P2")

cols4<-c("#DDAA33","#BB5566","#004488")

persite_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(type == "dedup_TE") %>%
  filter(virus == "Human_respiratory_syncytial_virus") %>%
  ggplot(aes(x=site,y=coverage,fill=Background))+
  geom_col()+
  facet_wrap(Viral.load~Sample_id)+
  ylab("Log(Coverage)")+
  ggtitle("Human RSV (deduplicated)")+
  guides(fill=guide_legend(title="Background"))+
  theme_bw()+
  scale_y_log10()+
scale_fill_manual(values=cols2,labels=c("ME_P1","ME_P2"))

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/RSV_coverage_persample.pdf",width=10,height=8)


persite_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(type == "dedup_TE") %>%
  filter(virus == "Zika_virus") %>%
  ggplot(aes(x=site,y=coverage,fill=Background))+
  geom_col()+
  facet_wrap(Viral.load~Sample_id)+
  ylab("Log(Coverage)")+
  ggtitle("Zika virus (deduplicated)")+
  guides(fill=guide_legend(title="Background"))+
  theme_bw()+
  scale_y_log10()+
  scale_fill_manual(values=cols2,labels=c("ME_P1","ME_P2"))

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/ZIKV_coverage_per_sample.pdf",width=10,height=8)

persite_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(type == "dedup_TE") %>%
  filter(virus == "Human_adenovirus_40") %>%
  ggplot(aes(x=site,y=coverage,fill=Background))+
  geom_col()+
  facet_wrap(Viral.load~Sample_id)+
  ylab("Log(Coverage)")+
  ggtitle("Human Adenovirus (deduplicated)")+
  guides(fill=guide_legend(title="Background"))+
  theme_bw()+
  scale_y_log10()+
  scale_fill_manual(values=cols2,labels=c("ME_P1","ME_P2"))

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/HADV_coverage_per_sample.pdf",width=20,height=18)

persite_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(type == "dedup_TE") %>%
  filter(virus == "Human_betaherpesvirus") %>%
  ggplot(aes(x=site,y=coverage,fill=Background))+
  geom_col()+
  facet_wrap(Viral.load~Sample_id)+
  ylab("Log(Coverage)")+
  ggtitle("Human Betaherpesvirus (deduplicated)")+
  guides(fill=guide_legend(title="Background"))+
  theme_bw()+
  scale_y_log10()+
  scale_fill_manual(values=cols2,labels=c("ME_P1","ME_P2"))

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/HBHV_coverage_per_sample_100.pdf",width=20,height=18)

#have to do these differently due to segmentation

FLU_ME1<-persite_norm %>%
  filter(virus == "Influenza_B_virus") %>%
  filter(type == "dedup_TE") %>%
  filter(Background == "p2") %>%
  ggplot(aes(x=site,y=coverage,fill=as.character(Viral.load)))+
  geom_col()+
  facet_grid(seg~Sample_id)+
  ylab("Log(Coverage)")+
  ggtitle("Influenza B virus (Background ME_1 deduplicated)")+
  guides(fill=guide_legend(title="Viral load"))+
  theme_bw()+
  scale_y_log10()+
  scale_fill_manual(values=cols4,labels=c("1e+02","1e+03","1e+05"))

FLU_ME2<-persite_norm %>%
  filter(virus == "Influenza_B_virus") %>%
  filter(type == "dedup_TE") %>%
  filter(Background == "p8") %>%
  ggplot(aes(x=site,y=coverage,fill=as.character(Viral.load)))+
  geom_col()+
  facet_grid(seg~Sample_id)+
  ylab("Log(Coverage)")+
  ggtitle("Influenza B virus (Background ME_2 deduplicated)")+
  guides(fill=guide_legend(title="Viral load"))+
  theme_bw()+
  scale_y_log10()+
  scale_fill_manual(values=cols4,labels=c("1e+02","1e+03","1e+05"))

ggarrange(FLU_ME1,FLU_ME2,ncol=2,common.legend=TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/FLUB_cov_per_sample.pdf",width=20,height=18)

REO_ME1<-persite_norm %>%
  filter(virus == "Mammalian_orthoreovirus3") %>%
  filter(type == "dedup_TE") %>%
  filter(Background == "p2") %>%
  ggplot(aes(x=site,y=coverage,fill=as.character(Viral.load)))+
  geom_col()+
  facet_grid(seg~Sample_id)+
  ylab("Log(Coverage)")+
  ggtitle("Mammalian orthoreovirus 3 (Background ME_1 deduplicated)")+
  guides(fill=guide_legend(title="Viral load"))+
  theme_bw()+
  scale_y_log10()+
  scale_fill_manual(values=cols4,labels=c("1e+02","1e+03","1e+05"))

REO_ME2<-persite_norm %>%
  filter(virus == "Mammalian_orthoreovirus3") %>%
  filter(type == "dedup_TE") %>%
  filter(Background == "p8") %>%
  ggplot(aes(x=site,y=coverage,fill=as.character(Viral.load)))+
  geom_col()+
  facet_grid(seg~Sample_id)+
  ylab("Log(Coverage)")+
  ggtitle("Mammalian orthoreovirus 3 (Background ME_2 deduplicated)")+
  guides(fill=guide_legend(title="Viral load"))+
  theme_bw()+
  scale_y_log10()+
  scale_fill_manual(values=cols4,labels=c("1e+02","1e+03","1e+05"))

ggarrange(REO_ME1,REO_ME2,ncol=2,common.legend=TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/REOV_cov_per_sample.pdf",width=20,height=18)
