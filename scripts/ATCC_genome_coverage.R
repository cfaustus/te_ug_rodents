#plots of TE experiment data - genome coverage & per site coverage
#also separated by viruses
#ATCC genomes updated version
#October 2024

#use bowtie2 data

rm(list=ls())

library(tidyverse)
library(ggpubr)

####genome coverage####

dedup_per_site<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_site_sample_and_virus_bowtie2_dedup_atcc_ref.tsv", sep = "\t", header = TRUE)

nodedup_per_site<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_site_sample_and_virus_bowtie2_nodedup_atcc_ref.tsv", sep = "\t", header = TRUE)

#add length column from read depth file to per site file

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

cols2<-c("#BB5566","#004488")

cols4<-c("#DDAA33","#BB5566","#004488")


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

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/Figure2.pdf",width=15,height=6)

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

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/FigureS5.pdf",width=10,height=8)

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

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/FigureS6.pdf",width=10,height=8)

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

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/FigureS7.pdf",width=20,height=18)

persite_norm %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(type == "dedup_TE") %>%
  filter(virus == "Human_betaherpesvirus") %>%
  ggplot(aes(x=site,y=coverage,fill=Background))+
  geom_polygon()+
  facet_wrap(Viral.load~Sample_id)+
  ylab("Log(Coverage)")+
  ggtitle("Human Betaherpesvirus (deduplicated)")+
  guides(fill=guide_legend(title="Background"))+
  theme_bw()+
  scale_y_log10()+
  scale_fill_manual(values=cols2,labels=c("ME_P1","ME_P2"))

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/HBHV_coverage_test.pdf",width=20,height=18)

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

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/FigureS9.pdf",width=20,height=18)

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

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/FigureS10.pdf",width=20,height=18)
