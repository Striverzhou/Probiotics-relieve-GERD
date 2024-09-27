rm(list = ls())
#############################################################################
# vOTU sankey network
votu_info1 <- read.table("vOTU_family_tax.txt",header = T,sep = "\t")
votu_info2 <- read.table("annotated_vOTU_info.txt",header = T,sep = "\t")

votu_info <- votu_info1 %>% mutate(
  group = ifelse(ictv_family=="NULL","unclassfied","classfied")
) %>% left_join(.,votu_info2[,c("qname","checkv_quality")],by="qname")

votu_info$ictv_family <- gsub("NULL",NA,votu_info$ictv_family)

votu_info_tb <- data.table::setDT(votu_info[,c("group","checkv_quality","ictv_family")])[,list(Count=.N),names(votu_info[,c("group","checkv_quality","ictv_family")])]

votu_info_data2 <- data.table::setDT(votu_info[,c("checkv_quality","ictv_family")])[,list(count=.N),names(votu_info[,c("checkv_quality","ictv_family")])]
votu_info_data2 <- na.omit(votu_info_data2) %>%
  set_colnames(c("source","target","value"))

nodes2 <- data.frame(names=c(as.character(votu_info_data2$source),
                            as.character(votu_info_data2$target)) %>%
                      unique()
)
# nodes2 <- data.frame(names=c(
#   "Medium-quality","High-quality","Complete","Siphoviridae",
#   "Podoviridae","Papillomaviridae","Myoviridae",
#   "Microviridae","Metaviridae","Inoviridae","Herelleviridae",
#   "crAss-phage","Circoviridae","Adenoviridae","Ackermannviridae"
# ))

votu_info_data2$IDsource <- match(votu_info_data2$source,nodes2$names)-1
votu_info_data2$IDtarget <- match(votu_info_data2$target,nodes2$names)-1

votu_info_data2 <- votu_info_data2[order(votu_info_data2$IDsource,
                                         votu_info_data2$IDtarget,
                                         decreasing = F),]
library(ggalluvial)
library(networkD3)
library(ggforce)
library(ggplot2)
library(ggpubr)
library(patchwork)

sankeyNetwork(Links = votu_info_data2,Nodes = nodes2,
                   Source = "IDsource",Target = "IDtarget",
                   Value = "value", NodeID = "names",#LinkGroup = 'source',
                   nodeWidth = 60,fontSize = 15,sinksRight = F)


votu_info_tb2 <- to_lodes_form(votu_info_data2,key = "x",value = "stratum",
                               id = "alluvium",axes = 1:2)

p1 <- ggplot(votu_info_tb2,aes(x=x,y=value,fill=stratum, label=stratum,
                               stratum = stratum, alluvium  = alluvium)) +
  geom_flow(width = 0.3,#linkline width
            curve_type = "sine",
            #curve shape，linear、cubic、quintic、sine、arctangent、sigmoid
            alpha = 0.5,
            color = 'white') +#gap color
  geom_stratum(width = 0.28, alpha=0.7,colour="white")+
  theme_void()

classfied_data <- as.data.frame(table(votu_info$group))
classfied_data$prop <- paste0(round(prop.table(table(votu_info$group)),4)*100,"%")

color_code <- c("#ff7f00","#1f78b4")
p2 <- ggplot(data = classfied_data,
             aes(x=1,y=Freq,fill=Var1)) +
  geom_tile() + scale_fill_manual(values = color_code)+
  theme_bw()+
  geom_text(aes(x=1,y=Freq,label=prop),color='black')+
  theme(panel.spacing = unit(0,"lines"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(fill=c("white"),color="gray"),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()
  )

p2+p1+plot_layout(widths = c(1,7),guides = 'collect')

ggsave(filename = "vOTU_composition_sankey.pdf",width = 8,height = 7)

##########################################################################
# Abundance stack diagram
#RPKM2Abundance
library(rstatix)
fp <- read.table("vOTU_family_abundance.txt",header = T,sep = "\t",row.names = 1)
metadata <- readxl::read_xlsx("grouping_info.xlsx",sheet=1)

fp <- scale(fp,center = F,scale = colSums(fp))*1000000
hp_data <- t(fp) %>% as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  left_join(metadata,by="sample")
hp_data$group <- factor(hp_data$group,levels = c("pla_0w","pro_0w","pla_4w","pro_4w",
                                                 "pla_8w","pro_8w","pla_12w","pro_12w"))

stack_input <- hp_data %>% group_by(group) %>%
  get_summary_stats(type = "mean")

# creat stack graph
library(ggsci)
col=pal_d3("category10")(10)

ggplot(stack_input) +
  geom_bar(aes(x=group,weight=mean,fill=reorder(variable,mean)),
           position = 'stack',width = 0.7) +
  scale_fill_manual(values = rev(col),name=NULL) +
  theme_bw() + theme(panel.grid = element_blank())+
  labs(x="",y="Relative abundance")+
  scale_y_continuous(breaks = c(0,250000,500000,750000,1000000),
                     labels = c("0%","25%","50%","75%","100%"))

ggsave(file = "vOTU_relative_abundance_stack_plot.pdf",width = 8,height = 6) 

###############################################################
# Differential vOTU abundance heatmap
votu <- read.table("diff_vOTU_abundance.txt",header = T,sep = "\t")

votu2 <- votu %>% column_to_rownames(var = "Votu") %>%
  select(is.numeric)
ann <- votu %>% select(Name,Votu) %>% column_to_rownames(var = "Votu") %>%
  arrange(Name) %>% set_colnames("family")
votu2 <- votu2[rownames(ann),]

library(pheatmap)
pheatmap(t(votu2),scale = "column",
         color = colorRampPalette(c("#313695","#4575b4","#74add1","#abd9e9",
                                    "white","#fee090","#fdae61","#f46d43",
                                    "#d73027"))(100),
         border_color = "white",
         annotation_col = ann,
         cellwidth = 20,cellheight = 20,
         cluster_rows = F,cluster_cols = F,
         # angle_col = 45,fontsize = 6,
         width = 18,height = 6,
         filename = "vOTU_abund_heatmap.pdf"
)



