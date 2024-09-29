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
# procruste plot
bins_table <- read.table("GERD_intestinal_flora__abundance.txt",header = T,sep = "\t",row.names = 1)
votu_table <- read.table("all_votu_rpkm.txt",header = T,sep = "\t",row.names = 1)
metadata <- readxl::read_xlsx("grouping_info.xlsx",sheet=1)
group <- metadata %>% separate(col = "group",into = c("Treat","Time"),sep = "_",remove = F)

spe.dist <- vegdist(t(bins_table)) # Bray-Curtis
env.dist <- vegdist(scale(t(votu_table),center = F,scale = colSums(t(votu_table))))

# Dimensionality reduction analysis
# NMDS
mds.s <- monoMDS(spe.dist)
mds.e <- monoMDS(env.dist)

# symmetric pattern
pro.s.e <- procrustes(mds.s,mds.e, symmetric = TRUE)
summary(pro.s.e)


# statistic test
pro.s.e_t <- protest(mds.s,mds.e, permutations = 999)
pro.s.e_t
#pro.s.e_t$ss

# p value
pro.s.e_t$signif


library(ggplot2)
# get the coordination
Pro_Y <- cbind(data.frame(pro.s.e$Yrot), data.frame(pro.s.e$X)) %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  left_join(group,by="sample")
Pro_X <- data.frame(pro.s.e$rotation)

ggplot(data=Pro_Y,) +
  geom_segment(aes(x = X1, y = X2, color = group,
                   xend = (X1 + MDS1)/2, yend = (X2 + MDS2)/2),
  ) +
  geom_segment(aes(x = (X1 + MDS1)/2, y = (X2 + MDS2)/2,
                   xend = MDS1, yend = MDS2,color = group)
  ) +
  geom_point(aes(X1, X2,color = group),shape=16,position="jitter") +
  geom_point(aes(MDS1, MDS2, color = group),shape=17,position="jitter") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = 'black',
                                        fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black', size=14),
        axis.title.y=element_text(colour='black', size=14),
        axis.text=element_text(colour='black',size=12)) +
  labs(x = 'MDS1', y = 'MDS2', color = '') +
  labs(title="Correlation between microbial community and metabolites") +
  theme(plot.title = element_text(size=8,colour = "black",
                                  hjust = 0.5,face = "bold"))

ggsave(filename = "procruste_plot2.pdf",width = 6,height = 5)

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



