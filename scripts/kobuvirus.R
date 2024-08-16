#kobuvirus data - plot read counts and read depths per sample
#August 2019

library(tidyverse)
library(ggpubr)
library(scales)

####read counts####

kobu_counts<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_kobuvirus/kobuvirus_TE_polyomics_readcount.tsv", sep = "\t", header=TRUE)

#normalise kobuvirus read counts by total counts per library

kobu_norm<-kobu_counts %>%
  mutate(norm_counts = matched/Number.of.paired.end.reads..QT.)

#number of viral reads vs total reads
#some relationship but not always the case that more reads = more depth (esp for non dedup)

kobu_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  ggplot(aes(x=Number.of.paired.end.reads..QT.,y=matched,colour=type))+
  geom_point()

kobu_norm %>%
  filter(mapper == "ngm") %>%
  filter(expt == "TE") %>%
  ggplot(aes(x=Number.of.paired.end.reads..QT.,y=matched,colour=type))+
  geom_point()

#for each non p6 sample, plot read count splitting by viral load & background

readcount_norm_split<-kobu_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=reorder(Sample_id,-norm_counts),y=norm_counts,fill = Background))+
  geom_col()+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  facet_grid(~as.character(Viral.load),scales="free")+
  ylab("Read count (normalised)")

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_readcount_norm_split.pdf")

#with p6

readcount_norm_p6<-kobu_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=reorder(Sample_id,-norm_counts),y=norm_counts,fill = Background))+
  geom_col()+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  facet_grid(~as.character(Viral.load),scales="free")+
  ylab("Read count (normalised)")

#combined across backgrounds

readcounts_norm_comb<-kobu_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  filter(Background != "control") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=as.character(Viral.load),y=norm_counts,colour = Background))+
  geom_boxplot()+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  ylab("Read count (normalised)")

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_readcounts_combined.pdf")

#read counts separated by background including p6
#p6 behaving more in an expected way here compared to virome mix data

readcounts_background<-kobu_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=reorder(Sample_id,-norm_counts),y=norm_counts,fill = as.character(Viral.load)))+
  geom_col()+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  facet_grid(~Background,scales="free")+
  ylab("Read count (normalised)")

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_readcounts_background_p6.pdf")

#compare bowtie2 vs ngm
#higher counts with ngm for all samples

kobu_norm %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=as.character(Viral.load),y=norm_counts,fill = mapper))+
  geom_col(position="dodge")+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  facet_grid(~Background,scales="free")

#compare dedup vs nondedup
#higher counts in nondedup samples

kobu_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=as.character(Viral.load),y=norm_counts,fill = type))+
  geom_col(position="dodge")+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  facet_grid(~Background,scales="free")

####read depth####
#first check that the normalised mean/median is the same as the mean/median of the normalised per site depth

kobu_depth<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_kobuvirus/kobuvirus_TE_polyomics_readdepth.tsv", sep = "\t", header = TRUE)

kobu_per_site<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_kobuvirus/kobuvirus_TE_polyomics_readdepth_persite.tsv", sep = "\t", header = TRUE)

kobu_depth_norm<-kobu_depth %>%
  mutate(norm_mean = mean_depth/Number.of.paired.end.reads..QT.) %>%
  mutate(norm_med = median_depth/Number.of.paired.end.reads..QT.)

#calculate the normalised coverage per site

norm_cov_per_site<-kobu_per_site %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  mutate(norm_cov = coverage/Number.of.paired.end.reads..QT.) %>%
  group_by(Sample_id,type) %>%
  summarise(mean_cov_site = mean(norm_cov), med_cov_site = median(norm_cov))

#check for just one sample

kobu_depth_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(Sample_id == "A")  %>%
  filter(type == "dedup")

kobu_per_site %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(Sample_id == "A") %>%
  filter(type == "dedup") %>%
  mutate(norm_cov = coverage/Number.of.paired.end.reads..QT.) %>%
  summarise(mean_cov_site = mean(norm_cov), med_cov_site = median(norm_cov))

#for a comparison across all samples, merge the datasets

compare_cov<-left_join(kobu_depth_norm,norm_cov_per_site,by=c("Sample_id","type"))

compare_cov_subset<- compare_cov %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE")

#after rounding the mean values, both are identical

identical(round(compare_cov_subset[['norm_mean']],digits=5),round(compare_cov_subset[['mean_cov_site']],digits=5))

identical(compare_cov_subset[['norm_med']],compare_cov_subset[['med_cov_site']])

compare_cov_subset %>%
  ggplot(aes(x=norm_mean,y=mean_cov_site))+
  geom_point()

compare_cov_subset %>%
  ggplot(aes(x=norm_med,y=med_cov_site))+
  geom_point()

#number of reads vs depth

kobu_depth_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
ggplot(aes(x=Number.of.paired.end.reads..QT.,y=mean_depth,colour=type))+
  geom_point()

kobu_depth_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  ggplot(aes(x=Number.of.paired.end.reads..QT.,y=median_depth,colour=type))+
  geom_point()

kobu_depth_norm %>%
  filter(mapper == "ngm") %>%
  filter(expt == "TE") %>%
ggplot(aes(x=Number.of.paired.end.reads..QT.,y=mean_depth,colour=type))+
  geom_point()

kobu_depth_norm %>%
  filter(mapper == "ngm") %>%
  filter(expt == "TE") %>%
  ggplot(aes(x=Number.of.paired.end.reads..QT.,y=median_depth,colour = type))+
  geom_point()

#for each non p6 sample, plot read depth splitting by viral load & background

mean_depth_norm<-kobu_depth_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=reorder(Sample_id,-norm_mean),y=norm_mean,fill = Background))+
  geom_col()+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  facet_grid(~as.character(Viral.load),scales="free")+
  ylab("Mean read depth (normalised)")

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_readdepth_mean_norm_split.pdf")

med_depth_norm<-kobu_depth_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=reorder(Sample_id,-norm_med),y=norm_med,fill = Background))+
  geom_col()+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  facet_grid(~as.character(Viral.load),scales="free")+
  ylab("Median read depth (normalised)")

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_readdepth_med_norm_split.pdf")

#with p6

depth_mean_p6<-kobu_depth_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  ggplot(aes(x=reorder(Sample_id,-norm_mean),y=norm_mean,fill = Background))+
  geom_col()+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  facet_grid(~as.character(Viral.load),scales="free")+
  ylab("Mean read depth (normalised)")

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_readdepth_mean_p6.pdf")

depth_med_p6<-kobu_depth_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  ggplot(aes(x=reorder(Sample_id,-norm_med),y=norm_med,fill = Background))+
  geom_col()+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  facet_grid(~as.character(Viral.load),scales="free")+
  ylab("Median read depth (normalised)")

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_readdepth_med_p6.pdf")

#summarised across samples

mean_depth_comb<-kobu_depth_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  filter(Background != "p6") %>%
ggplot(aes(x=as.character(Viral.load),y=norm_mean,colour = Background))+
  geom_boxplot()+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  ylab("Mean read depth (normalised)")

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_readdepth_mean_combined.pdf")

med_depth_comb<-kobu_depth_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  filter(Background != "p6") %>%
  ggplot(aes(x=as.character(Viral.load),y=norm_med,colour = Background))+
  geom_boxplot()+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  ylab("Median read depth (normalised)")

ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/kobuvirus_readdepth_med_combined.pdf")

#plot all samples including p6 split by background

kobu_depth_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  ggplot(aes(x=reorder(Sample_id,-norm_mean),y=norm_mean,fill = as.character(Viral.load)))+
  geom_col()+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  facet_grid(~Background,scales="free")+
  scale_y_continuous(labels = scales::label_comma())

kobu_depth_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  ggplot(aes(x=reorder(Sample_id,-norm_med),y=norm_med,fill = as.character(Viral.load)))+
  geom_col()+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  facet_grid(~Background,scales="free")+
  scale_y_continuous(labels = scales::label_comma())

#compare bowtie2 vs ngm
#higher depth with ngm for all samples

kobu_depth_norm %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  ggplot(aes(x=as.character(Viral.load),y=norm_med,fill = mapper))+
  geom_col(position="dodge")+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  facet_grid(~Background,scales="free")+
  scale_y_continuous(labels = scales::label_comma())

kobu_depth_norm %>%
  filter(expt == "TE") %>%
  filter(type == "dedup") %>%
  ggplot(aes(x=as.character(Viral.load),y=norm_mean,fill = mapper))+
  geom_col(position="dodge")+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  facet_grid(~Background,scales="free")+
  scale_y_continuous(labels = scales::label_comma())

#compare dedup vs nondedup
#more in nondedup samples

kobu_depth_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  ggplot(aes(x=as.character(Viral.load),y=norm_mean,fill = type))+
  geom_col(position="dodge")+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  facet_grid(~Background,scales="free")+
  scale_y_continuous(labels = scales::label_comma())

kobu_depth_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(expt == "TE") %>%
  ggplot(aes(x=as.character(Viral.load),y=norm_med,fill = type))+
  geom_col(position="dodge")+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  facet_grid(~Background,scales="free")+
  scale_y_continuous(labels = scales::label_comma())

####TE vs shotgun####

#for each sample, plot read count splitting by viral load & background

kobu_norm %>%
  filter(mapper == "bowtie2") %>%
  filter(type == "dedup") %>%
  filter(Background != "control") %>%
  ggplot(aes(x=reorder(Sample_id,-norm_counts),y=norm_counts,fill = Background))+
  geom_col()+
  theme(axis.title.x=element_blank(),axis.title.y=element_text(size=10))+
  facet_grid(~as.character(Viral.load),scales="free")+
  ylab("Read count (normalised)")

kobu_norm %>%
  filter(expt == "polyomics")
