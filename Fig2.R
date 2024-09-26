rm(list=ls())
################################################################
library(tidyverse)
library(reshape2)
library(vegan)
library(ggplot2)
library(ggpubr)

# read the intestinal flora abundance abundance data of GERD
bins_table <- read.table("GERD_intestinal_flora__abundance.txt",header = T,sep = "\t",row.names = 1)

# grouping information
meta <- readxl::read_xlsx("grouping_info.xlsx")
meta <- meta %>% separate(group,into = c("Treat","Time"),sep = "_")

# Alpha diversity
## Shannon index
A.shannon=diversity(t(bins_table),index="shannon")

## Simpson index
A.Simpson=diversity(t(bins_table),index="simpson")

A.div = data.frame(A.shannon,A.Simpson) %>%
  rownames_to_column(var = "sample") %>%
  left_join(meta,by="sample") %>%
  separate(group,into = c("Treat","Time"),remove = F)

A.div$Time <- factor(A.div$Time,levels = c("0w","4w","8w","12w"))

# Creating boxplot
color_code <- c("#0072b5ff","orange")
## Shannon index
stat.test <- A.div %>% group_by(Time) %>%
  wilcox_test(A.shannon~Treat) %>%
  add_xy_position(x = "Time", dodge = 0.85)

ggplot(A.div,aes(x=Time,y=A.shannon,color=Treat)) +
  geom_boxplot(alpha=0.7,outlier.shape = NA,
               position = position_dodge(0.85))+
  geom_jitter(shape=16,show.legend = FALSE,
              position = position_jitterdodge(0.2,dodge.width = 0.85))+
  scale_color_manual(values = color_code) +
  scale_y_continuous(name="Shannon index") +
  scale_x_discrete(name = "") +
  theme_classic() +
  theme(axis.text.x = element_text(vjust = 0.5)) +
  stat_pvalue_manual(stat.test,label = "p",hide.ns = T,
                     tip.length = 0.02)
ggsave(filename = "Gut_flora_Shannon_boxplot.pdf",width = 6,height = 5)


# Simpson index
stat.test <- A.div %>% group_by(Time) %>%
  wilcox_test(A.Simpson~Treat) %>%
  add_xy_position(x = "Time", dodge = 0.85)

ggplot(A.div,aes(x=Time,y=A.Simpson,color=Treat)) +
  geom_boxplot(alpha=0.7,outlier.shape = NA,
               position = position_dodge(0.85))+
  geom_jitter(shape=16,show.legend = FALSE,
              position = position_jitterdodge(0.2,dodge.width = 0.85))+
  scale_color_manual(values = color_code) +
  scale_y_continuous(name="Simpson index") +
  scale_x_discrete(name = "") +
  theme_classic() +
  theme(axis.text.x = element_text(vjust = 0.5)) +
  stat_pvalue_manual(stat.test,label = "p",hide.ns = T,
                     tip.length = 0.02)
ggsave(filename = "Gut_flora_Simpson_boxplot.pdf",width = 6,height = 5)


# Beta diversity
bray = vegdist(t(bins_table),method = "bray")
bray <- as.matrix(bray)

pcoa_plot <- function(metadata=meta,time,bray=bray,bins_table=bins_table,
                      colors=color_code,filename){
  subdata <- metadata %>% filter(.,endsWith(group,time)) %>%
    separate(col = group,into = c("Treat","Time"),sep = "_",remove = F)
  
  pcoa <- cmdscale(bray[subdata$sample,subdata$sample],k=3,eig = T)
  poi <- as.data.frame(pcoa$points) %>%
    set_colnames(paste0("PCoA",1:3)) %>%
    rownames_to_column(var = "sample") %>%
    left_join(subdata,by="sample")

  eig_percent <- round(pcoa$eig/sum(pcoa$eig)*100,1)
  aspe <- t(bins_table[,subdata$sample])
  div_adonis <- adonis2(aspe ~ Treat,data = subdata,permutations = 999,method = "bray")
  adonis_R2 <- round(div_adonis$R2[1],2)
  adonis_p <- div_adonis$`Pr(>F)`[1]

  ggplot(data = poi,aes(x=PCoA1,y=PCoA2))+
    geom_point(aes(color=group,shape=group),size=6)+
    scale_colour_manual(values = colors,name=NULL)+
    labs(x=paste("PCoA 1 (",eig_percent[1],"%)",sep = ""),
         y=paste("PCoA 2 (",eig_percent[2],"%)",sep = ""),
         colour=NULL,shape=NULL
    )+
    theme_bw()+
    geom_text(aes(x=max(PCoA1),y=min(PCoA2),
                  label=paste0("adonis==",adonis_R2)),
              hjust=1,vjust=-0.25,
              size=8,parse = T)+
    geom_text(aes(x=max(PCoA1),y=min(PCoA2),
                  label=paste0("italic(P)==",adonis_p)),
              hjust=1,vjust=1,
              size=8,parse = T) +
    theme(axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          panel.grid = element_blank(),
          legend.text = element_text(size = 25),
          plot.margin=unit(c(0.2,0.2,2,2),'lines'),
          axis.text = element_text(size = 20),
          legend.position = c(0.9,0.9))

  ggsave(filename = filename,width = 12,height = 12)
}

## 0w
pcoa_plot(metadata=meta,time="0w",bray=bray,bins_table=bins_table,
          colors=color_code,filename="pcoa_0w.pdf")
## 4w
pcoa_plot(metadata=meta,time="4w",bray=bray,bins_table=bins_table,
          colors=color_code,filename="pcoa_4w.pdf")

## 8w
pcoa_plot(metadata=meta,time="8w",bray=bray,bins_table=bins_table,
          colors=color_code,filename="pcoa_8w.pdf")

## 12w
pcoa_plot(metadata=meta,time="12w",bray=bray,bins_table=bins_table,
          colors=color_code,filename="pcoa_12w.pdf")


#########################################################################################
# circle heatmap
all_df_bins <- read.table("diff_bins.txt",header = T,sep = "\t")
df_bins_table <- t(bins_table) %>% as.data.frame() %>%
  select(all_df_bins$ID) %>%
  rownames_to_column(var = "sample") %>%
  left_join(meta,by="sample")

# calculate average value
df_bins_mean <- df_bins_table %>%
  group_by(group) %>%
  summarise_if(is.numeric, mean) %>%
  column_to_rownames(var = "group")

#SGB tax information
sgb_tax <- read.table("bins_tax_info.txt",sep = "\t",row.names = 1,header = F)
colnames(sgb_tax) <- c("Kingdom", "phylum", "class", "order", "family", "genus", "species")
df_bins_names <- gsub("^s__","",sgb_tax[colnames(df_bins_mean),"species"])

df_bins_tax <- sgb_tax %>% select(family,species) %>%
  filter(rownames(.) %in% all_df_bins$ID) %>%
  separate(col = family,into = c(NA,"family"),sep="__") %>%
  separate(col = species,into = c(NA,"species"),sep = "__") %>%
  rownames_to_column(var = "ID") %>%
  arrange(family)

# set family level color
library(RColorBrewer)
colors_ann <- c(brewer.pal(12,"Set3"),brewer.pal(9,"Set1")[3],brewer.pal(9,"Set1")[9])
family_ann <- colors_ann[match(df_bins_tax[,"family"],unique(df_bins_tax[,"family"]))]
names(family_ann) <- df_bins_tax$family
# scale
mat_mean <- t(df_bins_mean) %>% as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  left_join(df_bins_tax,by="ID")


mat_mean2 <- mat_mean %>% column_to_rownames(var = "species") %>%
  select(is.numeric) %>% 
  t() %>% as.data.frame() %>% scale()
  
  
cirle_inner <- mat_mean2 %>% t() %>%
  as.data.frame() %>%
  select(pla_0w,pro_0w,pla_4w,pro_4w,pla_8w,pro_8w)

cirle_outer <- mat_mean2 %>% t() %>%
  as.data.frame() %>%
  select(pla_12w,pro_12w)

#plot circle heatmap
library(circlize)
library(dendextend)
library(gridBase)
library(ComplexHeatmap)

col_fun1 = colorRamp2(c(-2.5,0,2.5),c("#74B9F7","white","#FF931E")) #for circle_inner 
col_fun2 = colorRamp2(c(-2.5,0,2.5),c("steelblue","white","#fcaf17")) #for circle_outer

####################################################################################
pdf("DiffBins_circle_heatmap.pdf",width = 8,height = 6)
circos.clear()
circos.par(gap.after = c(30))##gap for add colnames
# from the outside in
circos.heatmap(cirle_outer[df_bins_tax$species,],col=col_fun2,bg.lwd = 0.5,cluster = F,
               bg.border = "black",rownames.side = "outside",track.height=0.08)

circos.track(track.index = get.current.track.index(),
             
             panel.fun = function(x, y) {
               
               if(CELL_META$sector.numeric.index == 1) { # the last sector
                 
                 cn = colnames(cirle_outer)
                 
                 n = length(cn)
                 
                 circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), ##x-axis coordinate
                             
                             1:n - convert_y(0.5, "mm"), ##y-axis coordinate
                             
                             cn, ##colnames
                             
                             cex = 0.4, ##size of colnames
                             
                             adj = c(0, 0.5),
                             
                             facing = "inside")
                 
               }
               
             }, bg.border = NA)

circos.heatmap(cirle_inner[df_bins_tax$species,],
               track.margin=c(0.003,0.01),
               col = col_fun1, ## color
               track.height=0.21,
               #cell_width = 0.2,
               cluster = F,
               dend.side = "none",##tree
               bg.border = "black",bg.lwd = 0.5,
               #rownames.side = "outside",
               dend.track.height = 0.2,
               dend.callback = function(dend, m, si) {
                 # when k = 1, it renders one same color for the whole dendrogram
                 color_branches(dend, k = 4, col = 2:5)
               }
)
circos.track(track.index = get.current.track.index(),
             
             panel.fun = function(x, y) {
               
               if(CELL_META$sector.numeric.index == 1) { # the last sector
                 
                 cn = colnames(cirle_inner)
                 
                 n = length(cn)
                 
                 circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                             
                             1:n - convert_y(0.5, "mm"), 
                             
                             cn, 
                             
                             cex = 0.4,
                             
                             adj = c(0, 0.5),
                             
                             facing = "inside")
                 
               }
               
             }, bg.border = NA)

circos.heatmap(df_bins_tax[,"family"],col = family_ann,track.height=0.05,track.margin=c(0,0),cluster = F)

# For legend
lg1 <- Legend(title = "Relative abundance",col_fun = col_fun1,
              #size = unit(1.5,"mm"),
              title_gp = gpar(fontsize=8),
              direction = c("horizontal"),title_position = "topcenter")
lg2 <- Legend(title = "Relative abundance",col_fun = col_fun2,
              #size = unit(1.5,"mm"),
              title_gp = gpar(fontsize=8),
              direction = c("horizontal"),title_position = "topcenter")
lg3 <- Legend(title = "Family", at = unique(names(family_ann)),
              labels_gp = gpar(fontsize=8),
              legend_gp = gpar(fill = unique(family_ann)))

#grid.draw(lg)
h=dev.size()
gd_list= packLegend(lg1,lg2,lg3,max_height = unit(2*h, "inch"))

circle.size = unit(0.02,"snpc")
draw(gd_list,x=circle.size,just="left")

dev.off()

###############################################################################
# Relationship between flora abundance and clinical score
id_map <- read.table("id_mapping.txt",header = T,sep = "\t")
colnames(id_map) <- c("randomID","drugID","GID","group")

rdq <- read.table("RDQ_mean_score.txt",header = F,sep = "\t")
rdq2 <- rdq %>% select(V1,V3:V6) %>%
  set_colnames(c("sample","T0w","T4w","T8w","T12w"))

sgb64wk4 <- bins_table["SGB.64",which(endsWith(colnames(bins_table),"B"))]

sgb64wk4 <- t(sgb64wk4) %>% as.data.frame() %>%
  rownames_to_column(var = "GID") %>%
  mutate_at(vars(GID),~str_split(.,"B",simplify=TRUE)[,1])
sgb4w <- merge(id_map,sgb64wk4,by="GID") %>% rename(.,sample=drugID)

#Merge bins and RDQ
sgb_rdq4w_pro <- merge(sgb4w,rdq2,by="sample") %>% filter(group=="pro")

ggplot(data = sgb_rdq4w_pro,aes(x = SGB.64,y = T4w))+
  geom_point(colour="#FDBA6B")+
  scale_y_continuous(name = "RDQ score") +
  scale_x_continuous(name="Butyricimonas virosa")+
  theme_bw()+geom_smooth(
    method = "lm",se=T,
    color='cornflowerblue',
    size=1.5,fill="lightgray")+
  theme(legend.position = "none",panel.grid = element_blank(),
        axis.title=element_text(size=13,face="plain",color="black"),
        axis.text = element_text(size=12,face="plain",color="black"),
        plot.title = element_text(size = 10,hjust = 0.5))+
  stat_cor(method = "pearson",digits = 3,size=4)
ggsave(filename = "SGB64_RDQ_4w_relative_scatter.pdf",width = 5,height = 5)



