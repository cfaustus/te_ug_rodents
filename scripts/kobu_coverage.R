#kobuvirus genome coverage
#August 2024

library(tidyverse)

kobu_per_site<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_kobuvirus/kobuvirus_TE_polyomics_readdepth_persite.tsv", sep = "\t", header = TRUE)

#TE data only - compare coverage between backgrounds/viral loads

kobu_per_site %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=site,y=coverage,fill=Background))+
  geom_col(position="dodge")+
  facet_wrap(~Sample_id)

kobu_per_site %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=site,y=coverage,fill=as.character(Viral.load)))+
  geom_col(position="dodge")+
  facet_wrap(~Sample_id)

TE_coverage<-kobu_per_site %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=site,y=coverage))+
  geom_col(position="dodge")+
  facet_grid(as.character(Viral.load)~Background)

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_TE_coverage.pdf")

#TE vs shotgun

cov_TE_shotgun<-kobu_per_site %>%
  filter(mapper == "bowtie2") %>%
  filter(type == "dedup") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=site,y=coverage))+
  geom_col(position="dodge")+
  facet_grid(expt~Background)

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_TE_shotgun_coverage.pdf")

#TE vs shotgun with p6

cov_TE_shotgun_p6<-kobu_per_site %>%
  filter(mapper == "bowtie2") %>%
  filter(type == "dedup") %>%
  ggplot(aes(x=site,y=coverage))+
  geom_col(position="dodge")+
  facet_grid(expt~Background)

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_TE_shotgun_coverage_p6.pdf")
