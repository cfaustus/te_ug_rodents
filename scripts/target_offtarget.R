##all viruses - target and off target

rm(list=ls())

library(tidyverse)

#metadata

metadata<-read.csv("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/metadata/sampleIDs_TESpikeIn.csv",header=TRUE)

metadata2<-metadata %>%
  select(Sample.ID,Background.sample,Viral.load,Number.of.read.pairs..quality.adaptor.trimmed.) %>%
  rename(Sample_id = Sample.ID) %>%
  rename(QC_reads = Number.of.read.pairs..quality.adaptor.trimmed.)

#import viral reads mapped (to calculate proportions)

viral_reads_dedup<-read.csv("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/total_virus_mapped_reads_per_sample_dedup_atcc_ref_20241108.csv",header=TRUE)

viral_reads_nodedup<-read.csv("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/total_virus_mapped_reads_per_sample_nodedup_atcc_ref_20241108.csv",header=TRUE)

#import contingency tables

contingency_dedup<-read.csv("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/contingency_table_mapped_virus_reads_per_family_genus_sample_dedup_atcc_ref_20241106.csv", header = TRUE)

contingency_nodedup<-read.csv("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/contingency_table_mapped_virus_reads_per_family_genus_sample_nodedup_atcc_ref_20241106.csv", header = TRUE)

contaminants<-c("Betacoronavirus","Alphainfluenzavirus","Gammaretrovirus")

dedup_long<-contingency_dedup %>%
  filter(!genus %in% contaminants) %>%
  filter(genus != "NA") %>%
  pivot_longer(cols = A:P,
               names_to = "Sample_id",
               values_to = "count") %>%
  mutate(target = case_when(
    genus == "Cyclovirus" ~ "off_target",
    genus == "Cardiovirus" ~ "off_target",
    genus == "Kobuvirus" ~ "off_target",
  .default = "on_target")) %>%
  mutate(genome_structure = case_when((genus == "Mastadenovirus"| genus == "Cytomegalovirus") ~ "DNA",
                                      ,.default = "RNA"))

no_dedup_long<-contingency_nodedup %>%
  filter(!genus %in% contaminants) %>%
  filter(genus != "NA") %>%
  pivot_longer(cols = A:P,
               names_to = "Sample_id",
               values_to = "count") %>%
  mutate(target = case_when(
    genus == "Cyclovirus" ~ "off_target",
    genus == "Cardiovirus" ~ "off_target",
    genus == "Kobuvirus" ~ "off_target",
    .default = "on_target")) %>%
  mutate(genome_structure = case_when((genus == "Mastadenovirus"| genus == "Cytomegalovirus") ~ "DNA",
                                      ,.default = "RNA"))


metadata_dedup<-left_join(dedup_long,metadata2,by="Sample_id")

dedup_viral_reads<-left_join(metadata_dedup,viral_reads_dedup,by="Sample_id") %>%
  mutate(total_reads = QC_reads*2) %>%
  mutate(prop_total = count/total_reads) %>%
  mutate(prop_viral = count/total_virus_reads)

metadata_no_dedup<-left_join(no_dedup_long,metadata2,by="Sample_id")

nodedup_viral_reads<-left_join(metadata_no_dedup,viral_reads_nodedup,by="Sample_id") %>%
  mutate(total_reads = QC_reads*2) %>%
  mutate(prop_total = count/total_reads) %>%
  mutate(prop_viral = count/total_virus_reads)

facet_names<-c("Msp_p2" = "ME_P1",
               "Msp_p8" = "ME_P2",
               "off_target" = "Off target viruses",
               "on_target" = "On target viruses")

cols<-c("#CCBB44","#332288","#EE7733","#66CCEE","#882255","#4477AA","#AA3377","#228833","#EE6677")

dedup_reads<-dedup_viral_reads %>%
  filter(Background.sample != "Neg_control") %>%
  filter(Background.sample != "Msp_p6") %>%
  ggplot(aes(x=Viral.load,y=log(count),colour = genus))+
  geom_point()+
  geom_smooth(aes(group=genus,linetype=genome_structure),se=FALSE)+
  facet_grid(target~Background.sample,labeller=as_labeller(facet_names)) +
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10),legend.title=element_blank(),legend.position="top")+
  scale_x_log10()+
  ylab("Log(Viral Reads)")+
  scale_color_manual(values=cols)

dedup_prop_viral<-dedup_viral_reads %>%
  filter(Background.sample != "Neg_control") %>%
  filter(Background.sample != "Msp_p6") %>%
  ggplot(aes(x=Viral.load,y=prop_viral,colour = genus))+
  geom_point()+
  geom_smooth(aes(group=genus,linetype=genome_structure),se=FALSE)+
  facet_grid(target~Background.sample,labeller=as_labeller(facet_names)) +
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10),legend.title=element_blank(),legend.position="top")+
  scale_x_log10()+
  scale_y_log10()+
  ylab("Proportion of viral reads")+
  scale_color_manual(values=cols)

dedup_prop_all<-dedup_viral_reads %>%
  filter(Background.sample != "Neg_control") %>%
  filter(Background.sample != "Msp_p6") %>%
  ggplot(aes(x=Viral.load,y=prop_total,colour = genus))+
  geom_point()+
  geom_smooth(aes(group=genus,linetype=genome_structure),se=FALSE)+
  facet_grid(target~Background.sample,labeller=as_labeller(facet_names)) +
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10),legend.title=element_blank(),legend.position="top")+
  scale_x_log10()+
  scale_y_log10()+
  ylab("Proportion of total reads")+
  scale_color_manual(values=cols)

ggarrange(dedup_reads,dedup_prop_all,dedup_prop_viral,nrow=3,common.legend = TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/target_offtarget_dedup.png",width=8,height=10)

##just plot viral reads and prop viral

ggarrange(dedup_reads,dedup_prop_viral,nrow=2,common.legend = TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/target_offtarget_dedup.png",width=10,height=8)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/Figure3.pdf",width=10,height=8)

##extras - no dedup

nodedup_reads<- nodedup_viral_reads %>%
  filter(Background.sample != "Neg_control") %>%
  filter(Background.sample != "Msp_p6") %>%
  ggplot(aes(x=Viral.load,y=log(count),colour = genus))+
  geom_point()+
  geom_smooth(aes(group=genus,linetype=genome_structure),se=FALSE)+
  facet_grid(target~Background.sample,labeller=as_labeller(facet_names)) +
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10),legend.title=element_blank(),legend.position="top")+
  scale_x_log10()+
  ylab("Log(Viral Reads)")+
  scale_color_manual(values=cols)

nodedup_prop_viral<-nodedup_viral_reads %>%
  filter(Background.sample != "Neg_control") %>%
  filter(Background.sample != "Msp_p6") %>%
  ggplot(aes(x=Viral.load,y=prop_viral,colour = genus))+
  geom_point()+
  geom_smooth(aes(group=genus,linetype=genome_structure),se=FALSE)+
  facet_grid(target~Background.sample,labeller=as_labeller(facet_names)) +
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10),legend.title=element_blank(),legend.position="top")+
  scale_x_log10()+
  scale_y_log10()+
  ylab("Proportion of viral reads")+
  scale_color_manual(values=cols)

nodedup_prop_all<-nodedup_viral_reads %>%
  filter(Background.sample != "Neg_control") %>%
  filter(Background.sample != "Msp_p6") %>%
  ggplot(aes(x=Viral.load,y=prop_total,colour = genus))+
  geom_point()+
  geom_smooth(aes(group=genus,linetype=genome_structure),se=FALSE)+
  facet_grid(target~Background.sample,labeller=as_labeller(facet_names)) +
  theme_bw()+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10),legend.title=element_blank(),legend.position="top")+
  scale_x_log10()+
  scale_y_log10()+
  ylab("Proportion of total reads")+
  scale_color_manual(values=cols)

ggarrange(nodedup_reads,nodedup_prop_all,nodedup_prop_viral,nrow=3,common.legend = TRUE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/target_offtarget_nodedup.pdf",width=8,height=10)



