rm(list=ls())
################################################################
library(tidyverse)
library(magrittr)
library(rstatix)
library(ggplot2)
library(ggpubr)

# Code for RDQ
## ITT data set
rdq_itt <- readxl::read_xlsx("RDQ score.xlsx",sheet = "ITT")

#Convert the average score to the total score
rdq_itt <- rdq_itt %>% mutate_if(is.numeric,~.*8)

# statistic test
rdq_itt_stat <- rdq_itt %>% select(-ID) %>%
  group_by(Group) %>%
  get_summary_stats(`0W`,`4W`,`8W`,`12W`,type = "mean_se") %>%
  rename(Time=variable)


rdq_itt_long <- rdq_itt %>% select(-ID) %>%
  pivot_longer(cols = -Group, names_to = "Time", values_to = "Score")
rdq_itt_long$Time <- factor(rdq_itt_long$Time,levels = c("0W","4W","8W","12W"))

stat.test <- rdq_itt_long %>% group_by(Time) %>%
  wilcox_test(Score~Group) %>% add_x_position(x="Time")


#create line graph
ggplot(rdq_itt_stat,aes(x=Time,y=mean,fill=Group,colour=Group,group=Group,shape=Group))+
  geom_point(size=2) + geom_line()+
  scale_shape_manual(values = c(21,24))+scale_fill_manual(values=c("#2B6A99","#F16C23"))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.1)+
  scale_y_continuous(name='RDQ Score') +
  scale_x_discrete(name = "",)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0,
                                   vjust = 0.5
        )
  )+scale_color_manual(values = c("#2B6A99","#F16C23"))+
  stat_pvalue_manual(data=stat.test,x="x",y.position=11.8,label = "italic(P)~'='~{p}",
                     hide.ns=T,inherit.aes = FALSE,parse = T,remove.bracket = T)

ggsave(filename = "RDQ_ITT_plot.pdf",width = 6,height = 5)

## PPS data set
rdq <- readxl::read_xlsx("RDQ score.xlsx",sheet = "PPS")

# statistic test
rdq_stat <- rdq %>% select(-ID) %>%
  group_by(Group) %>%
  get_summary_stats(`0W`,`4W`,`8W`,`12W`,type = "mean_se") %>%
  rename(Time=variable)

rdq_long <- rdq %>% select(-ID) %>%
  pivot_longer(-Group,names_to = "Time",values_to = "Score") 
rdq_long$Time <- factor(rdq_long$Time,levels = c("0W","4W","8W","12W"))

stat.test <- rdq_long %>% group_by(Time) %>%
  wilcox_test(Score~Group) %>% add_x_position(x="Time") %>%
  mutate(p=round(p,digits = 3))

#create line graph
ggplot(rdq_stat,aes(x=Time,y=mean,fill=Group,colour=Group,group=Group,shape=Group))+
  geom_point(size=2) + geom_line()+
  scale_shape_manual(values = c(21,24),labels=c("pla","pro"))+
  scale_fill_manual(values=c("#2B6A99","#F16C23"),labels=c("pla","pro"))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.1)+
  scale_y_continuous(name='RDQ Score') +
  scale_x_discrete(name = "",)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0,
                                   vjust = 0.5
        )
  )+scale_color_manual(values = c("#2B6A99","#F16C23"),labels=c("pla","pro"))+
  stat_pvalue_manual(data=stat.test,x="x",y.position=11.8,label = "italic(P)~'='~{p}",
                     hide.ns=T,inherit.aes = FALSE,parse = T,remove.bracket = T)

ggsave(filename = "RDQ_PPS_plot.pdf",width = 6,height = 5)


# Code for GSRS ITT set data
## create violin plot
gsrs_violin_plot <- function(data="GSRS score.xlsx",sheet,delcol="sample",ylab,filename){
  gsrs <- readxl::read_xlsx(data,sheet = sheet)
  gsrs_long <- gsrs %>% select(-!!sym(delcol)) %>%
    set_colnames(c("group","0W","4W","8W","12W")) %>%
    pivot_longer(cols = -group, names_to = "Time", values_to = "Score")
  gsrs_long$Time <- factor(gsrs_long$Time,levels = c("0W","4W","8W","12W"))
  
  stat.test <- gsrs_long %>% group_by(Time) %>%
    wilcox_test(Score~group) %>%
    add_xy_position(x = "Time", dodge = 0.9)
  
  if(min(stat.test$p)<0.05){
    groupmax <- gsrs %>% select(which(stat.test$p < 0.05)+2) %>%
      summarise_all(max)
  }else{
    groupmax <- 0
  }
  
  if((max(gsrs_long$Score)-max(groupmax))<0.2){
    ggplot(gsrs_long,aes(x=Time,y=Score)) +
      geom_violin(aes(fill=group),
                  position = position_dodge(width = 0.9),
                  alpha = 0.7,trim = F,scale = "area") +
      geom_boxplot(aes(color=group),width=0.1,
                   outlier.shape = NA,
                   position = position_dodge(0.9))+
      scale_color_manual(values = rep("black",2),name="group",
                         labels=c("pla","pro")) +
      scale_y_continuous(name=ylab,expand = expansion(mult = c(0.05, 0.15))) +
      scale_x_discrete(name = "") +
      theme_classic() +
      theme(axis.text.x = element_text(vjust = 0.5)) +
      scale_fill_manual(values = c("#2B6A99","#F16C23"),
                        labels=c("pla","pro")) +
      stat_pvalue_manual(stat.test,label = "p",hide.ns = T,
                         tip.length = 0.02,bracket.nudge.y = 0.85)
  }else{
    ggplot(gsrs_long,aes(x=Time,y=Score)) +
      geom_violin(aes(fill=group),
                  position = position_dodge(width = 0.9),
                  alpha = 0.7,trim = F,scale = "area") +
      geom_boxplot(aes(color=group),width=0.1,
                   outlier.shape = NA,
                   position = position_dodge(0.9))+
      scale_color_manual(values = rep("black",2),name="group",
                         labels=c("pla","pro")) +
      scale_y_continuous(name=ylab) +
      scale_x_discrete(name = "") +
      theme_classic() +
      theme(axis.text.x = element_text(vjust = 0.5)) +
      scale_fill_manual(values = c("#2B6A99","#F16C23"),
                        labels=c("pla","pro")) +
      stat_pvalue_manual(stat.test,label = "p",hide.ns = T,
                         tip.length = 0.02)
  }
  ggsave(filename = filename,width = 6,height = 5)
}

## Total GSRS score
gsrs_violin_plot(sheet = "Total",ylab = "Total GSRS score",filename = "Total_GSRS_violin_plot.pdf")

## Constipation syndrome score
gsrs_violin_plot(sheet = "Constipation",ylab = "Constipation syndrome score",filename = "Constipation_violin_plot.pdf")

## Reflux syndrome score
gsrs_violin_plot(sheet = "Reflux",ylab = "Reflux syndrome score",filename = "Reflux_violin_plot.pdf")

## Abdominal pain score
gsrs_violin_plot(sheet = "Abdominal_pain",ylab = "Abdominal pain score",filename = "AbdominalPain_violin_plot.pdf")

## Diarrhoea syndrome score
gsrs_violin_plot(sheet = "Diarrhea", ylab = "Diarrhoea syndrome score",filename = "Diarrhoea_violin_plot.pdf")

## Indigestion syndrome score
gsrs_violin_plot(sheet = "Dyspepsia",ylab = "Indigestion syndrome score", filename = "Indigestion_violin_plot.pdf")
