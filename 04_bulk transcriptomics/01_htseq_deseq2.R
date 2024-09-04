#=====================================================================================
#  Load package
#=====================================================================================

library(DESeq2)
library(ggplot2)
library(plotly)
library(grid)

#=====================================================================================
#  Import data
#=====================================================================================
path <- '/Users'
count_data <- read.table(paste0(path,"/merge_b7.htcount.txt"),header = T, sep = "\t", row.names = 1)
colnames(count_data) <- c("C1", "A1", "B1",
                          "A2", "C2", "B2", 
                          "C5", "B5", "A5",
                          "B6", "C6", "A6")

count_data <- as.matrix(count_data)
count_data[is.na(count_data)] <- 0

#=====================================================================================
#  DESeq2 analysis
#=====================================================================================
condition <- factor(c("c", "a", "b",
                      "a", "c", "b",  
                      "c", "b", "a",
                      "b", "c", "a" ))


dds <- DESeqDataSetFromMatrix(count_data, DataFrame(condition), design= ~condition )
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

#normlization matrix
rld <- rlogTransformation(dds)
expr_data <- assay(rld) 
write.csv(expr_data,paste0(path,'/expr.csv'))

control <- 'b'
test <- 'a'

res <- results(dds,contrast=c("condition",test,control))
diff_stat <- as.data.frame(res)
diff_stat[which(diff_stat$padj < 0.05 & diff_stat$log2FoldChange >= 1),'diff'] <- 'up'
diff_stat[which(diff_stat$padj < 0.05 & diff_stat$log2FoldChange <= -1),'diff'] <- 'down'
diff_stat[!(diff_stat$diff %in% c('up', 'down')),'diff'] <- 'no'
diff_stat$id <- row.names(diff_stat)
write.csv(diff_stat,paste0(path,'/a-b.csv'))

#=====================================================================================
#  Drawing volcano map
#=====================================================================================

deseq2_plot <- function(diff_stat){
  mytheme <- 
    theme_bw()+
    theme(legend.text = element_text(color = 'gray40',size = rel(0.9)),
          legend.title = element_text(color = 'gray20',size = rel(1)))+
    theme(panel.border = element_rect(color='gray20'),
          panel.grid.major.y = element_blank(), 
          panel.grid.minor.y = element_blank(), 
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          panel.background = element_rect(color = 'gray',fill = 'transparent'))
  
  ggplot(diff_stat, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = diff),size = 2) +
    scale_colour_manual(limits = c('up', 'down', 'no'), values = c("#B688FB","#f6d250",'gray')) +
    labs(x = paste0('log2 (',control,'/',test,')'), y = '-log10 padj') +
    geom_vline(xintercept = c(-1, 1), color = 'gray70', linewidth = 0.5,linetype = "dashed") + 
    geom_hline(yintercept = -log10(0.05), color = 'gray70', linewidth = 0.5,linetype = "dashed")+
    geom_text(data = subset(diff_stat, diff!= 'no'), aes(label = id), vjust = -1, hjust = 0.5, size = 3, color = 'black') +
    mytheme+
    xlim(c(-5, 5))
}

{
  viewport(layout.pos.row=x,layout.pos.col=y)
}

p <- deseq2_plot(diff_stat)
p <- p+ggtitle(paste0(test,'-',control))
p <- ggplotly(p)
p
htmlwidgets::saveWidget(as.widget(p), paste0(path,'\\',test,'-',control,'.html'))

pdf(paste0(path,'\\',test,'-',control,'.pdf'),width = 20,height = 20)
print(p)


dev.off()

