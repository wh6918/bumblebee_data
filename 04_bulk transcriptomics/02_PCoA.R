#=====================================================================================
#  Load package
#=====================================================================================
library(vegan);cat("future:",as.character(packageVersion("vegan")),"\n")
library(grid);cat("future:",as.character(packageVersion("grid")),"\n")
library(tidyr);cat("future:",as.character(packageVersion("tidyr")),"\n")
library(ggplot2);cat("future:",as.character(packageVersion("ggplot2")),"\n")
library(readr);cat("future:",as.character(packageVersion("readr")),"\n")
library(dplyr);cat("future:",as.character(packageVersion("dplyr")),"\n")

#=====================================================================================
#  Import data
#=====================================================================================
setwd("/Users")
data <- read.csv('./expr.csv',header = T)
head(data)
gene <- data$ID
row.names(data) <- gene
data <- select(data,-"ID")
data <- t(data)
ma <- as.matrix(data)
sample <- row.names(ma)
group <- c("C", "A", "B",
           "A", "C", "B", 
           "C", "B", "A",
           "B", "C", "A" )
df_group <- data.frame(sample,group)

#=====================================================================================
#   PCOA analysis
#=====================================================================================
distance <- vegdist(ma,method='bray')


pcoa <- cmdscale(distance,eig = T)
pcoa.data <- data.frame(pcoa$points)[,1:2]
pcoa.data$sample <- row.names(pcoa.data)
pcoa.data <- merge(pcoa.data,df_group,by = 'sample',all.x = T)
names(pcoa.data)[2:3] <- c('PCOA1','PCOA2')

pcoa1 <- pcoa$eig[1]/sum(pcoa$eig)
pcoa2 <- pcoa$eig[2]/sum(pcoa$eig)

#=====================================================================================
#   Drawing
#=====================================================================================
mytheme <- 
  theme_bw()+
  theme(#axis.text.x = element_blank(),
    #axis.title =element_blank(),
    #axis.ticks.x = element_blank(),
    axis.line = element_blank()
  )+
  theme(legend.text = element_text(color = 'gray40',size = rel(0.9)),
        legend.title = element_blank()
        #legend.position = 'none'
  )+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.border = element_rect(color='gray20'),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(color = 'gray20', fill = 'gray90'))

p<-ggplot(pcoa.data,aes(x=PCOA1,y=PCOA2))+
  stat_ellipse(aes(color=group,fill = group), level = 0.95)+
  geom_point(aes(color=group),size = 2.5,alpha = 0.8)+
  labs(x=paste0('PCOA1(',round(pcoa1*100,1),'%)'),y=paste0('PCOA2(',round(pcoa2*100,1),'%)'))+
  #ggtitle(sdp)+
  scale_color_brewer(palette = 'Set2')+
  geom_vline(xintercept=0,size=0.2,alpha=0.2)+
  geom_hline(yintercept=0,size=0.2,alpha=0.2)+
  mytheme+
  scale_color_manual(values = c("#eda3a7", "#ddd6d0","#a5c5d3"))+
  scale_fill_manual(values = c("#eda3a7", "#ddd6d0","#a5c5d3"))


p <- p + geom_text(aes(label=sample), vjust=-0.5, hjust=0.5)
pdf('./pcoa.pdf',width = 5,height = 4)
print(p)
dev.off()


#=====================================================================================
#   PERMANOVA analysis
#=====================================================================================

permanova_results <- adonis(distance ~ group, data = df_group)
print(permanova_results)

F_value <- permanova_results$aov.tab$F.Model[1]  
P_value <- permanova_results$aov.tab$`Pr(>F)`[1]  



