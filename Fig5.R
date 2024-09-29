rm(list = ls())
################################################################################
library(ggplot2)
library(tidyr)
library(FactoMineR)
library(factoextra)
library(ropls)

# oplsda
tt <- readxl::read_xlsx("metobalite_raw.data.xlsx")

data <- tt %>% column_to_rownames(var = "Sample") %>%
  t() %>% as.data.frame() %>%
  mutate_at(vars(-Label),as.numeric)
oplsda_plot <- function(data=data,time,filename){
  set.seed(1234)
  data0 <- data %>% filter(.,endsWith(Label,time))
  subdata <- data0 %>% select(is.numeric)
  lv <- as.factor(data0$Label)
  oplsda <- opls(x=subdata,y=lv,predI = 1,orthoI = NA,fig.pdfC = "none")

  data2 <- as.data.frame(oplsda@scoreMN)
  o1 <- oplsda@orthoScoreMN[,1]
  data2$o1 <- o1
  data2$group <- lv
  data2$samples <- rownames(data2)
  
  x_lab <- oplsda@modelDF[,"R2X"]*100
  
  ggplot(data2,aes(x=p1,y=o1,color=group))+
    theme_bw()+
    geom_point(size=3)+
    theme(panel.grid = element_blank())+
    geom_vline(xintercept = 0,lty="dashed",color="grey")+
    geom_hline(yintercept = 0,lty="dashed",color="grey")+
    labs(x=paste0("T score[1] (",x_lab[1],"%)"),
         y=paste0("Orthogonal T score[1] (",x_lab[2],"%)"),
         title = "OPLS-DA") +
    stat_ellipse(linetype = 2)+
    scale_color_manual(values = c("#1b72b3","orange")) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid=element_blank())
  
  ggsave(filename = filename,width = 5,height = 5)
}

oplsda_plot(data = data,time = "0w",filename = "oplsda_0w.pdf")

oplsda_plot(data = data,time = "4w",filename = "oplsda_4w.pdf")

oplsda_plot(data = data,time = "8w",filename = "oplsda_8w.pdf")

oplsda_plot(data = data,time = "12w",filename = "oplsda_12w.pdf")

# VIP bubble plot
library(tidyverse)
library(magrittr)

rt <- readxl::read_xlsx("Metabolites_VIP.xlsx",col_names = F)

rt1 <- rt %>% set_colnames(c("Time","ID","VIP","pla","pro","name","M/Z")) %>%
  mutate(
    Time = case_when(
      Time == "4w_8w" ~ "4w & 8w",
      Time == "8w_12w" ~ "8w & 12w",
      Time == "4w_8w_12w" ~ "4w & 8w & 12w"
    ),Foldchange = pro/pla
  )

rt1$Time <- factor(rt1$Time,levels = c("4w","8w","12w","4w & 8w","8w & 12w","4w & 8w & 12w"))

p1 <- ggplot(rt1,aes(x=VIP, y=reorder(name,VIP),size=Foldchange,color=Time)) +
  geom_point(alpha=0.5) +
  scale_size(range = c(1, 8), name="Foldchange") +
  theme_bw() +
  scale_y_discrete(name=NULL,labels=NULL)+
  theme(panel.grid = element_blank(),
        panel.spacing = unit(0,"lines"),
        #panel.border = element_blank(),
        plot.margin = margin(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())

tt <- rt1 %>% arrange(VIP) %>%
  select(pla,pro,name) %>%
  column_to_rownames(var = "name")

#install.packages("ggheatmap")
library(ggheatmap)
p2 <- ggheatmap(tt,text_position_rows = "left",
                text_position_cols = "top",
                levels_rows = rownames(tt))+
  theme(panel.spacing = unit(0,"lines"),
        plot.margin = margin())

library(patchwork)
p2+p1+plot_layout(widths = c(1,5),guides = 'collect')

ggsave("untarget_diff_metabolite_VIP_bubble.pdf",width = 8,height = 8)  

