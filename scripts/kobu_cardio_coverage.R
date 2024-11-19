#kobuvirus and cardiovirus genome coverage
#November 2024

rm(list=ls())

library(tidyverse)

kobu_per_site<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_kobuvirus/kobuvirus_TE_polyomics_readdepth_persite_20241007.tsv", sep = "\t", header = TRUE) %>%
  select(virus,seg,site,coverage,n_sites,Sample_id,Background,Viral.load,mapper,type)


kobu_per_site$Viral.load<-kobu_per_site$Viral.load %>%
  replace_na(., 0)
  
kobu_per_site$Sample_id <- kobu_per_site$Sample_id %>%
  case_match("RNA-Msp-p2" ~ "ME_P1",
             "RNA-Msp-p8" ~ "ME_P2",
             .default = kobu_per_site$Sample_id)

cardio_per_site<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_cardiovirus/cardiovirus_denovo_bowtie2_read_depth_per_site_and_sample_dedup.tsv", sep = ",", header = TRUE)

persite_both<-rbind(kobu_per_site,cardio_per_site)

persite_both$type <- persite_both$type %>%
  case_match("dedup_TE" ~ "dedup",
             .default = persite_both$type)

#TE data only - compare coverage between backgrounds/viral loads

cols2<-c("#BB5566","#004488")

coverage_labels<-c("100" = "1e+02",
                    "1000" = "1e+03",
                    "100000" = "1e+05",
                    "0" = "Shotgun",
                   "Kobuvirus" = "Kobuvirus",
                   "Cardiovirus" = "Cardiovirus")

coverage_labels2<-c("100" = "1e+02",
                   "1000" = "1e+03",
                   "100000" = "1e+05",
                   "0" = "Shotgun",
                   "A" = "A",
                   "B" = "B",
                   "C" = "C",
                   "D" = "D",
                   "F" = "F",
                   "G" = "G",
                   "H" = "H",
                   "I" = "I",
                   "J" = "J",
                   "K" = "K",
                   "L" = "L",
                   "M" = "M",
                   "N" = "N",
                   "O" = "O",
                   "P" = "P",
                   "ME_P1" = "ME_P1",
                   "ME_P2" = "ME_P2")

kobu_coverage<-kobu_per_site %>%
  filter(mapper == "bowtie2") %>%
  filter(type == "dedup") %>%
  filter(Background != "p6") %>%
  group_by(Sample_id,Viral.load,Background) %>%
  filter(coverage > 0) %>%
  summarise(genome_coverage = n()) %>%
  mutate(percent_coverage = genome_coverage/8467) %>%
  mutate(virus = "Kobuvirus")

cardio_coverage<-cardio_per_site %>%
  filter(Background != "p6") %>%
  group_by(Sample_id,Viral.load,Background) %>%
  filter(coverage > 0) %>%
  summarise(genome_coverage = n()) %>%
  mutate(percent_coverage = genome_coverage/6945) %>%
  mutate(virus = "Cardiovirus")

coverage_both<-rbind(kobu_coverage,cardio_coverage)

#plot together

coverage_both %>%
  filter(Viral.load != 0) %>%
  ggplot(aes(x=Viral.load,y=percent_coverage,color=Background))+
  geom_point()+
  ylab("Genome coverage")+
  xlab("Spike In Viral Load")+
  facet_grid(~virus)+
  theme_bw()+
  scale_color_manual(values=cols2,labels=c("ME_P1","ME_P2"))+
  scale_x_log10()+
  scale_y_continuous(limits=c(0,1))

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_coverage_new//kobuvirus_cardiovirus_coverage.png")

#TE and shotgun

coverage_both %>%
  ggplot(aes(x=Sample_id,y=percent_coverage,color=Background))+
  geom_point()+
  ylab("Genome coverage")+
  xlab("Sample ID")+
  facet_grid(virus~Viral.load,scales="free_x",labeller=as_labeller(coverage_labels))+
  theme_bw()+
  scale_color_manual(values=cols2,labels=c("ME_P1","ME_P2"))+
  scale_y_continuous(limits=c(0,1))

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_coverage_new//kobuvirus_TE_shotgun_coverage.png")

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/FigureS14.pdf")

#per site genome coverage - Kobu

persite_both %>%
  filter(virus == "k97_19971") %>%
  filter(mapper == "bowtie2") %>%
  filter(type == "dedup") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=site,y=log(coverage),fill=Background))+
  geom_col()+
  theme_bw()+
  facet_wrap(Viral.load~Sample_id,labeller=as_labeller(coverage_labels2))+
  scale_fill_manual(values=cols2,labels=c("ME_P1","ME_P2"))+
  ylab("Log(Kobuvirus Coverage)")+
  xlab("") +
  ggtitle("Kobuvirus")

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_coverage_new//kobuvirus_TE_shotgun_coverage_per_site.pdf",width=12,height=7)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/FigureS15.pdf",width=12,height=7)

#per site genome coverage - Cardio

persite_both %>%
  filter(virus == "cardiovirus") %>%
  filter(mapper == "bowtie2") %>%
  filter(type == "dedup") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=site,y=log(coverage),fill=Background))+
  geom_col()+
  facet_wrap(Viral.load~Sample_id,labeller=as_labeller(coverage_labels2))+
  scale_fill_manual(values=cols2,labels=c("ME_P1","ME_P2"))+
  theme_bw()+
  ylab("Log(Cardiovirus Coverage)")+
  xlab("") +
  ggtitle("Cardiovirus")

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_coverage_new/cardio_TE_coverage_per_site.pdf",width=12,height=7)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/FigureS17.pdf",width=12,height=7)
