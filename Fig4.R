rm(list = ls())
###############################################################################
library(tidyverse)
library(rstatix)
library(ggplot2)
library(ggpubr)

# Module abundance heatmap
rt <- read.table("diff_module_abundance.txt",header = T,sep = "\t",)
data <- rt %>% column_to_rownames(var = "ID") %>%
  select(4:ncol(.)) %>% t() %>% as.data.frame() %>%
  rename(.,group=ID) %>%
  mutate_at(vars(starts_with("CQM")),as.numeric) %>%
  group_by(group) %>% summarise_if(is.numeric,mean,na.rm=T)

info <- rt[-1,1:4]

data2 <- data %>% mutate_at(vars(-group),scale) %>%
  pivot_longer(cols = -group,names_to = "module",
               values_to = "abundance") %>%
  mutate(type=str_split(group,"_",simplify=TRUE)[,2])

data2$group <- factor(data2$group,
                      levels = c("pla_0w","pro_0w","pla_4w","pro_4w",
                                 "pla_8w","pro_8w","pla_12w","pro_12w"))
data2$type <- factor(data2$type,
                     levels = c("0w","4w","8w","12w"))
data2$module <- factor(data2$module,
                       levels = info[order(info$module,info$ID),"ID"])

df3<-data.frame(
  x = seq(0.5,2.5,1),
  xend = seq(0.5,2.5,1),
  y = -Inf,
  yend = Inf
)
df3
df4<-data.frame(
  x = -Inf,
  xend = Inf,
  y = seq(0.5,14.5,1),
  yend = seq(0.5,14.5,1)
)

p1 <- ggplot(data=data2,aes(x=group,y=module,color=abundance))+
  geom_point(aes(size=abundance,
  ),shape=15) + scale_size(range = c(2,12))+
  scale_color_gradientn(colours = c("#313695","#4575b4","#74add1","#abd9e9",
                                    "white","#fee090","#fdae61","#f46d43",
                                    "#d73027","#a50026"))+
  scale_x_discrete(expand = c(0,0),name='',labels=c("pla","pro"))+
  facet_grid(.~type,scales="free",space = "free_x")+
  geom_segment(data=df3,aes(x=x,xend=xend,y=y,yend=yend),
               color="grey",linewidth=3)+
  geom_segment(data=df4,aes(x=x,xend=xend,y=y,yend=yend),
               color="grey",linewidth=3)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.margin = margin(),
        panel.border = element_rect(color="grey"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank())

color_code <- c("#ee716b","#c99e0d","steelblue","#5ab033","#37b28f","#a08ec3")
library(ggforce)
info$ID <- factor(info$ID,levels = info[order(info$module,info$ID),"ID"])
info$Percent <- as.numeric(info$Percent)
p2 <- ggplot(data = info[order(info$module,info$ID),],aes(x=1,y=ID))+
  geom_tile(aes(x=2,fill=module),color="grey")+
  geom_point(aes(x=1,size=Percent),color="orange")+
  scale_fill_manual(values = color_code)+
  scale_x_discrete(expand = expansion(mult = c(0.3,0)))+
  scale_y_discrete(expand = expansion(mult = c(0,0)))+
  theme_bw()+
  theme(panel.spacing = unit(0,"lines"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(),
        strip.background = element_rect(fill=c("white"),color="gray"),
        axis.title = element_blank(),
        axis.ticks = element_blank()
  )

library(patchwork)
p2+p1+plot_layout(widths = c(1,7),guides = 'collect')

ggsave(filename = "SGB_module_heatmap.pdf",width = 7.5,height = 7.5)

################################################################################
# differential metabolites from melonnpan
melon_metab <- readxl::read_xlsx("diff_melonnpan_metabolite.xlsx")

melon_metab2 <- melon_metab %>%
  separate(group,into = c("Group","Time"),sep = "_")

melon_metab2$Time <- factor(melon_metab2$Time,
                            levels = c("0w","4w","8w","12w"))

boxDotPlot_sc <- function(data,xvar,yvar,group,alpha,color,dotshp=21,yp){
  fo <- as.formula(paste(yvar,'~',group))
  stat_test <- data %>% group_by(!!sym(xvar)) %>%
    wilcox_test(.,formula=fo) %>%
    add_xy_position(x=xvar)
  
  p <- ggplot(data,aes(x=!!sym(xvar),y=!!sym(yvar)))+
    geom_boxplot(aes(fill=!!sym(group)),position = position_dodge(width = 0.9),
                 outlier.shape = NA,alpha=alpha)+
    geom_jitter(aes(fill=!!sym(group)),shape=dotshp,show.legend = FALSE,
                color="white",size=2.5,alpha=alpha,stroke=0.8,
                position = position_jitterdodge(0.1,dodge.width = 0.9)
    )+#scale_fill_jco()+
    scale_fill_manual(values = color)+
    scale_x_discrete(name = "")+
    scale_y_continuous(labels = scales::scientific)+
    theme_bw() + theme(panel.grid = element_blank())+
    stat_pvalue_manual(stat_test,label = "p",
                       y.position = yp,tip.length = 0.005)
  
  return(p)
}

#1
boxDotPlot_sc(data=melon_metab2,xvar="Time",yvar="chenodeoxycholate",
              group = "Group",alpha=0.7,color = color_code,yp=0.0012)
ggsave(filename = "chenodeoxycholate_boxplot.pdf",width = 6,height = 5.5)

#2
boxDotPlot_sc(data=melon_metab2,xvar="Time",yvar="phytosphingosine",
              group = "Group",alpha=0.7,color = color_code,yp=0.00013)
ggsave(filename = "phytosphingosine_boxplot.pdf",width = 6,height = 5.5)

#3
boxDotPlot_sc(data=melon_metab2,xvar="Time",yvar="citrulline",
              group = "Group",alpha=0.7,color = color_code,yp=0.0011)
ggsave(filename = "citrulline_boxplot.pdf",width = 6,height = 5.5)

#4
boxDotPlot_sc(data=melon_metab2,xvar="Time",yvar="glutamate",
              group = "Group",alpha=0.7,color = color_code,yp=0.0056)
ggsave(filename = "glutamate_boxplot.pdf",width = 6,height = 5.5)

###############################################################################
# NMDS
library(vegan)

mp <- read.table("MelonnPan_Predicted_Metabolites.txt",header = T,sep = "\t",row.names = 1,check.names = F)
meta <- readxl::read_xlsx("grouping_info.xlsx")
meta <- meta %>% separate(col = group,into = c("Treat","Time"),sep = "_",remove = F)

# BiocManager::install('mixOmics')
library(mixOmics)

samples_0w <- rownames(group[which(group$period=="0w"),])
samples_4w <- rownames(group[which(group$period=="4w"),])
samples_8w <- rownames(group[which(group$period=="8w"),])
samples_12w <- rownames(group[which(group$period=="12w"),])

group_0w <- group[samples_0w,]
group_4w <- group[samples_4w,]
group_8w <- group[samples_8w,]
group_12w <- group[samples_12w,]

library(vegan)
bray = vegdist(mp,method = "bray")
bray <- as.matrix(bray)

nmds_plot <- function(time,metadata=meta,bray=bray,filename){
  subdata <- metadata %>% filter(Time==time)
  nmds <- metaMDS(bray[subdata$sample,subdata$sample],k=2)
  stress <- nmds$stress
  stress_text <- paste("Stress =",round(stress,4))
  poi <- as.data.frame(nmds$points) %>%
    rownames_to_column(var = "sample") %>%
    left_join(metadata,by="sample")
  
  ggplot(data = poi,aes(x=MDS1,y=MDS2)) +
    geom_point(aes(color=group,shape=group),size=3) +
    scale_color_manual(values = c("blue","violetred"))+
    labs(x="NMDS1",y="NMDS2")+
    theme_bw() + theme(panel.grid = element_blank()) +
    geom_text(aes(x=max(MDS1),y=min(MDS2)),
              label=stress_text,size=3,hjust=1,vjust=0)
  
  ggsave(filename = filename,width = 6.5,height = 5)
}


nmds_plot(time = "0w",metadata = meta,bray = bray,filename = "NMDS_0w.pdf")
nmds_plot(time = "4w",metadata = meta,bray = bray,filename = "NMDS_4w.pdf")
nmds_plot(time = "8w",metadata = meta,bray = bray,filename = "NMDS_8w.pdf")
nmds_plot(time = "12w",metadata = meta,bray = bray,filename = "NMDS_12w.pdf")