#Fig6a cMET high volcano -----------------------------------------------------------------
load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/cMET_vs_wild_two_res.RData')
cMET_vs_wild_two_res$gene <- cMET_vs_wild_two_res$X
cMET_vs_wild_two_res$sig  <- ifelse(cMET_vs_wild_two_res$padj<0.05 & cMET_vs_wild_two_res$log2FoldChange > 1, 'up',
                                    ifelse(cMET_vs_wild_two_res$padj <0.05 & cMET_vs_wild_two_res$log2FoldChange < -1, 'dn', 'no'))
de_genes                  <- c('EP300', 'CREBBP', 'BMPR2','TGFBR2', 'JAK2','MET','KLK2','KLK4','KLK3')

ggplot(cMET_vs_wild_two_res, aes(x = log2FoldChange, y = -log10(padj),color=sig))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "black")+
  geom_point()+
  scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
  scale_size_continuous(range = c(0,1))+
  theme(legend.title = element_blank() )+
  geom_label_repel(data = subset(cMET_vs_wild_two_res, gene %in% de_genes),
                   max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                   aes(label = gene),
                   size = 5, 
                   color = 'black',fontface = "bold.italic",
                   segment.color = "black",   
                   segment.linetype = "solid", 
                   segment.size = 0.5 ,
                   box.padding = 3,        
                   point.padding = 0.2 ) +
  xlab("log2FC")+
  ylab("-log10 p-adj")+
  theme_bw()+
  ggtitle('')+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),     # 坐标轴标题加粗
    axis.text = element_text(size = 14, face = "bold"),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.2, "cm"),         # 增加刻度线长度
    axis.ticks = element_line(size = 1.2)        # 增加刻度线粗细
  )

ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/cMET.volcano.pdf', width = 8, height = 8)



# Fig.6b ------------------------------------------------------------------


load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/jak.score.RData')
jak.score %>% ggplot(aes(x= type, y = value, color = type))+
  geom_point(size = 7)+
  scale_color_manual(values =c(wild_type="Blue",cMET_high="Red"))+
  ylab('JAK/STAT score')+
  xlab('')+
  theme_bw()+
  scale_x_discrete(labels=(c('wild_type'='Control',
                             c('cMET_high'='VCAP-MET-OE'))))+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )


# Fig.6c ------------------------------------------------------------------

load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/MSPC.stem.cor.RData')
MSPC.stem.cor %>% ggplot(aes(x = engineered.cell.line, y = value))+
  
  geom_violin(aes(fill = engineered.cell.line, colour = engineered.cell.line), alpha = 0.5) +
  # alpha控制不透明度
  geom_boxplot(aes(colour = engineered.cell.line), width = 0.2)+
  scale_color_manual(values = c("PC3_PROSTATE" = "firebrick", 'MET_OE2' = '#FBCE6A', 'MET_OE1'='#FBCE6A', 'LNCAPKO2'='#974F9F',
                                'LNCAPKO1'='#974F9F','LNCaP.enza.1'='#354898','LNCaP.enza.2'='#354898','LNCaP.enza.3'='#354898'))+
  scale_fill_manual(values = c("PC3_PROSTATE" = "firebrick", 'MET_OE2' = '#FBCE6A', 'MET_OE1'='#FBCE6A', 'LNCAPKO2'='#974F9F',
                               'LNCAPKO1'='#974F9F','LNCaP.enza.1'='#354898','LNCaP.enza.2'='#354898','LNCaP.enza.3'='#354898'))+
  scale_x_discrete(labels = c('PC3_PROSTATE' = paste('PC3', sep = '\n', '(n = 2785)'),
                              'MET_OE2' = paste('VCaP-MET-OE1', sep = '\n', '(n = 2785)'),
                              'MET_OE1' = paste('VCaP-MET-OE2', sep = '\n', '(n = 2785)'),
                              'LNCAPKO2'=paste('LNCaP-TP53/RB1-KO1', sep = '\n','(n= 2785)'),
                              'LNCAPKO1'=paste('LNCaP-TP53/RB1-KO2', sep = '\n','(n= 2785)'),
                              'LNCaP.enza.1'=paste('LNCaP-enza1',sep = '\n','(n = 2785)'),
                              'LNCaP.enza.2'=paste('LNCaP-enza2',sep = '\n','(n = 2785)'),
                              'LNCaP.enza.3'=paste('LNCaP-enza3',sep = '\n','(n = 2785)')))+
  
  theme_niwot()+
  xlab('')+
  ylab('Transcriptome similarity')+
  RotatedAxis()+
  theme(legend.position = 'none',
        axis.title = element_text(size = 30, color = 'black'),
        axis.text.x  = element_text(size = 30, color = 'black'),
        axis.text.y  = element_text(size = 30, color = 'black'))+
  stat_compare_means(paired = T, comparisons = list(c('PC3_PROSTATE','MET_OE2')))

ggsave('./pcdata/paper/plot/section5/engineered.pdf', width = 28, height = 12)




# Fig.6d ------------------------------------------------------------------

col.fun    <- circlize::colorRamp2(seq(0.3,0.7, length.out = 101),colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(101))
p1 <- Heatmap(
  mat.en.cor,
  name = "Correlation",
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = T,
  row_names_side = "left",
  col = col.fun,
  column_names_rot = 315
  
)


pdf('/home/liuxueying/Proj_MetPCCelline/output/plot/lineage.annotation.pdf', width = 8, height = 8)
p1
dev.off()



# Fig.6e ------------------------------------------------------------------

load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/MSPC.stem.cor.RData')

MSPC.stem.cor %>% ggplot(aes(x = organoid, y = value))+
  
  geom_violin(aes(fill = organoid, colour = organoid), alpha = 0.5) +
  geom_boxplot(aes(colour = organoid), width = 0.2)+
  scale_color_manual(values = color_map) +
  scale_fill_manual(values = color_map)+
  scale_x_discrete(labels = c('PC3_PROSTATE' = paste('PC3', sep = '\n', '(n = 2785)'),
                              'MSKPCa12' = paste('MSKPCa12', sep = '\n', '(n = 2785)'),
                              'MSKPCa9' = paste('MSKPCa9', sep = '\n', '(n = 2785)'),
                              'MSKPCa13'=paste('MSKPCa13', sep = '\n','(n= 2785)'),
                              'MSKPCa18'=paste('MSKPCa18', sep = '\n','(n= 2785)'),
                              'MSKPCa8'=paste('MSKPCa8',sep = '\n','(n = 2785)'),
                              'MSKPCa3'=paste('MSKPCa3',sep = '\n','(n = 2785)'),
                              'MSKPCa11'=paste('MSKPCa11',sep = '\n','(n = 2785)'),
                              'MSKPCa17'=paste('MSKPCa17',sep = '\n','(n = 2785)'),
                              'MSKPCa15'=paste('MSKPCa15',sep = '\n','(n = 2785)'),
                              'MSKPCa20'=paste('MSKPCa20',sep = '\n','(n = 2785)')))+
  
  xlab('')+
  ylab('Transcriptome similarity')+
  RotatedAxis()+
  theme(legend.position = 'none',
        axis.title = element_text(size = 30, color = 'black'),
        axis.text.x  = element_text(size = 30, color = 'black'),
        axis.text.y  = element_text(size = 30, color = 'black'))+
  stat_compare_means(paired = T, comparisons = list(c('PC3_PROSTATE','MSKPCa12'),
                                                    c('PC3_PROSTATE','MSKPCa9'),
                                                    c('PC3_PROSTATE','MSKPCa13'),
                                                    c('PC3_PROSTATE','MSKPCa18')))



# Fig.6f ------------------------------------------------------------------

load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/mat.organoid.cor.RData')
col.fun    <- circlize::colorRamp2(seq(0.3,0.7, length.out = 101),colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(101))

p1 <- Heatmap(
  mat.organoid.cor,
  name = "Correlation",
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = T,
  right_annotation = row_anno,
  row_names_side = "left",
  col = col.fun,
  column_names_rot = 315
  
)

pdf('/home/liuxueying/Proj_MetPCCelline/output/plot/lineage.annotation.pdf', width = 8, height = 8)
p1
dev.off()


# Fig.6g-h ----------------------------------------------------------------

load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/stem.pc3.data.log2cpm.symbol.RData')
df            <- stem.pc3.data.log2cpm.symbol['TP63',] %>% sort() %>% as.data.frame()
colnames(df)  <- 'expression'
df$sample     <- rownames(df)

df            <- df[order(df$expression, decreasing = T),]
df$sample     <- factor(df$sample, levels = rev(df$sample))

ggplot(df, aes(x = sample, y = expression)) +
  geom_bar(stat = "identity", fill = "#D8CDEB") +
  
  coord_flip() +
  labs(title = "",
       x = "",
       y = "")+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 20, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 20, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.25, "cm"),         # 增加刻度线长度
    
  )
ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/TP63_expr.pdf', width = 8, height = 6)

df              <- stem.pc3.data.log2cpm.symbol['KRT5',] %>% sort() %>% as.data.frame()
colnames(df)    <- 'expression'
df$sample       <- rownames(df)

df <- df[order(df$expression, decreasing = T),]
df$sample <- factor(df$sample, levels = rev(df$sample))

ggplot(df, aes(x = sample, y = expression)) +
  geom_bar(stat = "identity", fill = "#D8CDEB") +
  
  coord_flip() +
  labs(title = "",
       x = "",
       y = "")+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 20, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 20, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.25, "cm"),         # 增加刻度线长度
    
  )
ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/KRT5_expr.pdf', width = 8, height = 6)














