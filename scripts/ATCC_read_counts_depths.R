#plots of TE experiment data - read counts and read depths
#ATCC genomes updated version
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

#plot viral read counts (non normalised)

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

#plot just the two lefthand panels from the original figure

ggarrange(read_count,read_count_norm2,nrow=2,common.legend = TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/log_read_count_depth_2panels.png",width=10,height=7)

##deduplicated only

#plot viral read counts (non normalised)

read_count_dedup<- counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(type == "dedup_TE") %>%
  ggplot(aes(x=Viral.load,y=log(matched),colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus,linetype=genome_structure),se=FALSE)+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10),legend.title=element_blank(),legend.position="top")+
  scale_x_log10()+
  scale_color_manual(values=cols,labels=c("Human adenovirus 40","Human betaherpesvirus","Human respiratory syncytial virus","Influenza B virus","Mammalian orthoreovirus 3","Zika virus"))+
  ylab("Log(Viral Reads)")

#normalise by raw read count and genome length

read_count_norm2_dedup<- counts_reads_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(type == "dedup_TE") %>%
  ggplot(aes(x=Viral.load,y=log(norm_counts2),colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus,linetype=genome_structure),se=FALSE)+
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10),legend.title=element_blank())+
  scale_x_log10()+
  scale_color_manual(values=cols,labels=c("Human adenovirus 40","Human betaherpesvirus","Human respiratory syncytial virus","Influenza B virus","Mammalian orthoreovirus 3","Zika virus"))+
  ylab("Log(viral reads/cleaned reads/genome length)")

#plot just the two lefthand panels from the original figure

ggarrange(read_count_dedup,read_count_norm2_dedup,nrow=2,common.legend = TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/log_read_count_depth_2panels_dedup.png",width=10,height=7)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/Figure1.pdf",width=10,height=7)

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

ggarrange(read_count,read_count_norm1,read_count_norm2,read_depths,nrow=2,ncol=2,common.legend=TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/log_read_count_depth.png",width=10,height=7)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/FigureS1.pdf",width=10,height=7)

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

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/FigureS4.pdf",width=10,height=7)

####extras####

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

##just genome length

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

ggarrange(read_count_raw,raw_read_count_norm1,raw_read_count_norm2,raw_read_depths,nrow=2,ncol=2,common.legend=TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/read_count_depth.pdf")

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
