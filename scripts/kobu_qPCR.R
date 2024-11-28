#Kobuvirus qPCR data
#November 2024

rm(list=ls())

library(tidyverse)

kobu_qpcr<-read.csv("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_kobuvirus/qPCR_data.csv",header=TRUE)

kobu_avg<-kobu_qpcr %>%
  group_by(Sample_name) %>%
  mutate(avg_Ct = mean(c(Ct1,Ct2,Ct3)),
         sd_Ct = sd(c(Ct1,Ct2,Ct3)))

kobu_long <- kobu_avg %>%
  group_by(Sample_name) %>%
  pivot_longer(cols = c(Ct1,Ct2,Ct3),names_to="replicate",values_to="Ct_value")

expt<-c("Msp-p2","Msp-p8")

kobu_test<-kobu_long %>%
  select(Sample_name,replicate,Ct_value) %>%
  filter(Sample_name %in% expt)

t.test(Ct_value~Sample_name,kobu_test)

#check with kobu depth

kobu_depth<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_kobuvirus/kobuvirus_TE_polyomics_readdepth_20241007.tsv", sep = "\t", header = TRUE) %>%
  select(!Number.of.paired.end.reads..QT.) %>%
  filter(mapper == "bowtie2") %>%
  filter(Background != "p6")

kobu_depth$Viral.load <- kobu_depth$Viral.load %>%
  replace_na(., 0)

expt2<-c("p2","p8")

depth_summary<-kobu_depth %>%
  filter(type == "dedup") %>%
  group_by(Background,Viral.load) %>%
  summarise(mean_depth = mean(mean_depth ))

depth_test <- kobu_depth %>%
  filter(Viral.load == 0)

