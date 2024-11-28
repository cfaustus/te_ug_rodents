#plots of Genome Medicine experiment data
#November 2024

#viral load vs read count normalised using different methods
#use bowtie2 data

rm(list=ls())

library(tidyverse)
library(ggpubr)

####read counts/normalised read counts####

#import and combine read count files

gm_readcount_dedup<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_genome_medicine/genome_medicine_TE_sequencing_experiment_readcount_per_sample_and_virus_bowtie2_dedup_atcc_ref.tsv", sep = "\t", header = TRUE)

gm_readcount_nodedup<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_genome_medicine/genome_medicine_TE_sequencing_experiment_readcount_per_sample_and_virus_bowtie2_nodedup_atcc_ref.tsv", sep = "\t", header = TRUE)

gm_counts_all<-rbind(gm_readcount_dedup,gm_readcount_nodedup)

#normalise by both raw read count and genome length - should be the same as mean read depth

facet_names<-c("dedup_GM_TE" = "Deduplicated",
               "nodedup_GM_TE"= "Non-Deduplicated")

cols<-c("#4477AA","#66CCEE","#228833","#CCBB44","#EE6677","#AA3377")

gm_counts_norm<-gm_counts_all %>%
  mutate(norm_counts3 = matched/length)

read_counts<-gm_counts_norm %>%
  filter(Viral.load != 0) %>%
  ggplot(aes(x=as.character(Viral.load),y=log(matched),colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus),se=FALSE)+
  facet_grid(~type,labeller=as_labeller(facet_names))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10),legend.title=element_blank(),legend.position="top")+
  scale_color_manual(values=cols,labels=c("Human adenovirus 40","Human betaherpesvirus","Human respiratory syncytial virus","Influenza B virus","Mammalian orthoreovirus 3","Zika virus"))+
  ylab("Log(Viral Reads)")+
  xlab("Spike In Viral load")

#import and combine read depth files

depth_dedup_gm<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_genome_medicine/genome_medicine_TE_sequencing_experiment_readdepth_per_sample_and_virus_bowtie2_dedup_atcc_ref.tsv", sep = "\t", header = TRUE)

depth_nodedup_gm<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_genome_medicine/genome_medicine_TE_sequencing_experiment_readdepth_per_sample_and_virus_bowtie2_nodedup_atcc_ref.tsv", sep = "\t", header = TRUE)

depths_gm_all<-rbind(depth_dedup_gm,depth_nodedup_gm)

read_depths<-depths_gm_all %>%
  filter(Viral.load != 0) %>%
  ggplot(aes(x=as.character(Viral.load),y=log(mean_depth),colour = virus))+
  geom_point()+
  geom_smooth(aes(group=virus),se=FALSE)+
  facet_grid(~type,labeller=as_labeller(facet_names))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10),legend.title=element_blank(),legend.position="top")+
  scale_color_manual(values=cols,labels=c("Human adenovirus 40","Human betaherpesvirus","Human respiratory syncytial virus","Influenza B virus","Mammalian orthoreovirus 3","Zika virus"))+
  ylab("Log(mean read depth)")+
  xlab("Spike In Viral load")

ggarrange(read_counts,read_depths,nrow=2,common.legend = TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/genome_medicine/read_counts_depths.png")

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/FigureS11.pdf")

####genome coverage####

dedup_per_site<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_genome_medicine/genome_medicine_TE_sequencing_experiment_readdepth_per_site_sample_and_virus_bowtie2_dedup_atcc_ref.tsv", sep = "\t", header = TRUE)

#View(dedup_per_site)

counts_dedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readcount_per_sample_and_virus_bowtie2_dedup_atcc_ref.tsv", sep = "\t", header = TRUE)

#add length column from read depth file to per site file

lengths<-counts_dedup_bt %>%
  select(virus,length) %>%
  distinct()

#calculate genome coverage

coverage<-dedup_per_site %>%
  group_by(sample_id,virus,Viral.load,Replicate) %>%
  filter(coverage > 0) %>%
  summarise(genome_coverage = n()) %>%
  left_join(lengths) %>%
  mutate(percent_coverage = genome_coverage/length)

coverage %>%
  filter(Viral.load != 0) %>%
  ggplot(aes(x=virus,y=percent_coverage))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position = position_dodge(width = .75))+
  facet_grid(~as.character(Viral.load))+
  theme_bw()+
  theme(axis.title.y=element_text(size=10))+
  ylab("Genome coverage")+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=15),axis.text=element_text(size=12),strip.text=element_text(size=12))+
  scale_y_continuous(limits=c(0,1))

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/genome_medicine/all_genome_coverage.png",width=14,height=6)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/FigureS12.pdf",width=20,height=6)

####per site plot AdV only####

dedup_per_site %>%
  filter(virus == "Human_adenovirus_40") %>%
  filter(Viral.load != 0) %>%
  ggplot(aes(x=site,y=coverage))+
  geom_col()+
  facet_grid(Viral.load~Replicate)+
  ylab("Log(Coverage per site)")+
  ggtitle("Human Adenovirus (deduplicated)")+
  theme_bw()+
  scale_y_log10()+
  xlab("")

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/genome_medicine/HADV_genome_coverage.pdf",width=20,height=16)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/FigureS14.pdf",width=20,height=16)
