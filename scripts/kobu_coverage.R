#kobuvirus genome coverage
#August 2024

library(tidyverse)

kobu_per_site<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_kobuvirus/kobuvirus_TE_polyomics_readdepth_persite.tsv", sep = "\t", header = TRUE)

#TE data only - compare coverage between backgrounds/viral loads

#calculate the normalised coverage per site

norm_cov_per_site<-kobu_per_site %>%
  mutate(norm_cov = coverage/Number.of.paired.end.reads..QT.)

norm_cov_per_site %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=site,y=norm_cov,fill=Background))+
  geom_col(position="dodge")+
  facet_wrap(~Sample_id)

TE_coverage_sep<-norm_cov_per_site %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=site,y=norm_cov,fill=as.character(Viral.load)))+
  geom_col(position="dodge")+
  facet_wrap(~Sample_id)+
  ylab("Normalised coverage")+
  guides(fill=guide_legend(title="Viral load"))

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_TE_coverage_separate.pdf")

TE_coverage<-norm_cov_per_site %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=site,y=norm_cov))+
  geom_col(position="dodge")+
  facet_grid(as.character(Viral.load)~Background)+
  ylab("Normalised coverage")

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_TE_coverage.pdf")

#TE vs shotgun

cov_TE_shotgun<-norm_cov_per_site %>%
  filter(mapper == "bowtie2") %>%
  filter(type == "dedup") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=site,y=norm_cov))+
  geom_col(position="dodge")+
  facet_grid(expt~Background)+
  ylab("Normalised coverage")

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_TE_shotgun_coverage.pdf")

#TE vs shotgun with p6

cov_TE_shotgun_p6<-norm_cov_per_site %>%
  filter(mapper == "bowtie2") %>%
  filter(type == "dedup") %>%
  ggplot(aes(x=site,y=norm_cov))+
  geom_col(position="dodge")+
  facet_grid(expt~Background)+
  ylab("Normalised coverage")

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_TE_shotgun_coverage_p6.pdf")
