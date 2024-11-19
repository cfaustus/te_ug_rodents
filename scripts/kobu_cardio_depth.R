#kobuvirus cardiovirus read depth
#November 2024

rm(list=ls())

library(tidyverse)
library(ggpubr)
library(scales)

kobu_depth<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_kobuvirus/kobuvirus_TE_polyomics_readdepth_20241007.tsv", sep = "\t", header = TRUE) %>%
  select(!expt) %>%
  select(!Number.of.paired.end.reads..QT.) %>%

cardio_depth<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_cardiovirus/cardiovirus_denovo_bowtie2_read_depth_per_sample_dedup.tsv",sep =",",header=TRUE)

depth_both<-rbind(kobu_depth,cardio_depth)

depth_both$type <- depth_both$type %>%
  case_match("dedup_TE" ~ "dedup",
             .default = depth_both$type)

cols2<-c("#BB5566","#004488")

coverage_labels<-c("k97_19971" = "Kobuvirus",
                   "cardiovirus" = "Cardiovirus",
                   "p2" = "ME_P1",
                   "p8" = "ME_P2")

depth_both %>%
  filter(Background != "p6") %>%
  filter(type == "dedup") %>%
  filter(mapper == "bowtie2") %>%
  filter(Viral.load != "NA") %>%
  ggplot(aes(x=Viral.load,y=log(mean_depth),colour = Background))+
  geom_point()+
  theme_bw()+
  facet_grid(Background~virus,labeller=as_labeller(coverage_labels))+
  theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_text(size=10),legend.position="none")+
  scale_color_manual(values=cols2)+
  scale_x_log10()+
  ylab("Log(mean read depth)")

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_coverage_new/kobu_cardio_depth.png")

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/FigureS16.pdf",width=12,height=7)
