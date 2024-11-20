#Linear model for TE data - separated by virus
#November 2024

#viral read depth
#use bowtie2 data

rm(list=ls())

library(tidyverse)
library(broom)

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

####read depths####

#import and combine read depth files

depth_dedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_sample_and_virus_bowtie2_dedup_atcc_ref.tsv", sep = "\t", header = TRUE)

depth_nodedup_bt<-read.table("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/data_TE/TE_sequencing_experiment_readdepth_per_sample_and_virus_bowtie2_nodedup_atcc_ref.tsv", sep = "\t", header = TRUE)

depths_bt_all<-rbind(depth_dedup_bt,depth_nodedup_bt)

depths_reads<-left_join(depths_bt_all,metadata2,by="Sample_id")%>%
  mutate(genome_structure = case_when((virus == "Human_adenovirus_40"| virus == "Human_betaherpesvirus") ~ "DNA",
                                      ,.default = "RNA"))

####linear model for spike in viruses - mean read depth####

depths_reads_sub<-depths_reads %>%
  filter(Background != "p6") %>%
  filter(Background != "control") %>%
  filter(type == "dedup_TE") %>%
  mutate(log_depth = log(mean_depth)) %>%
  mutate(log_viral_load = log10(Viral.load))

#inspect data
hist(depths_reads_sub$mean_depth)
hist(depths_reads_sub$Viral.load)
hist(depths_reads_sub$log_depth)
hist(depths_reads_sub$log_viral_load)

summary(depths_reads_sub$log_depth)
summary(depths_reads_sub$log_viral_load)

#all viruses together

lm1<-glm(log_depth~log_viral_load*virus,data=depths_reads_sub,family="gaussian")

summary(lm1)

glm.diag.plots(lm1)

pred<-predict(lm1,type="response")

rsq <-function (x,y) cor(x,y)^2

rsq(depths_reads_sub$log_depth,pred)

##split by viruses

cols<-c("#4477AA","#66CCEE","#228833","#CCBB44","#EE6677","#AA3377")

facet_names<-c("Human_adenovirus_40" = "Human adenovirus 40",
               "Human_betaherpesvirus"= "Human betaherpesvirus",
               "Human_respiratory_syncytial_virus"= "Human respiratory syncytial virus",
               "Influenza_B_virus"= "Influenza B virus",
               "Mammalian_orthoreovirus3"= "Mammalian orthoreovirus 3",
               "Zika_virus"= "Zika virus")

virus<-(unique(depths_reads_sub$virus)) %>%
  rep(.,each=2) %>%
  data.frame()

list_models<-depths_reads_sub %>%
  group_split(virus) %>%
  map(~lm(log_depth~log_viral_load,data=.))

lm_tidy<-map(list_models,broom::tidy) %>%
  do.call(rbind.data.frame,.) %>%
  cbind(virus,.) %>%
  rename(virus = ".") %>%
  select(virus,term,estimate) %>%
  pivot_wider(names_from = "term",values_from = "estimate")

lm_summary<-depths_reads_sub %>%
  group_by(virus) %>%
  summarise(Intercept = lm(log_depth~log_viral_load)$coefficients[1],
            Coeff_x1 = lm(log_depth~log_viral_load)$coefficients[2],
            R2 = summary(lm(log_depth~log_viral_load))$r.squared,
            pvalue = summary(lm(log_depth~log_viral_load))$coefficients["log_viral_load",4])

lm_combined<-left_join(lm_tidy,lm_summary,by="virus")

ggplot() +
  geom_point(data=depths_reads_sub,mapping=aes(log_viral_load,log_depth,color=virus))+
  geom_abline(data = lm_combined,aes(intercept = `(Intercept)`,slope = log_viral_load)) +
  geom_text(data=lm_combined,colour="black",aes(x=4.7,y=3,label=round(Coeff_x1,digits=2)) )+
  geom_text(data=lm_combined,colour="black",aes(x=4,y=3,label="Slope ="))+
  geom_text(data=lm_combined,colour="black",aes(x=4.7,y=2,label=round(R2,digits=2)) )+
  geom_text(data=lm_combined,colour="black",aes(x=4,y=2,label="R2 ="))+
  facet_wrap(~virus,labeller=as_labeller(facet_names)) +
  theme_bw() +
  xlab("Log(Viral Load)") +
  ylab("Log(Mean Read Depth)")+
  scale_color_manual(values=cols)+
  guides(color = FALSE)

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/compare_spike_ins_atcc/model_readdepth_dedup.png")

#ggsave("/Users/laura/Dropbox/glasgow/github/te_ug_rodents/figures/manuscript_figures_pdf/FigureS5.pdf")
