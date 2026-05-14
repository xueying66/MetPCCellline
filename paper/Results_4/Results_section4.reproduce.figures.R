
# Fig.4a ------------------------------------------------------------------
load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/pca.rna.RData')

ggplot(pca.rna,aes(x= PC1, y = PC2, color = cell.line))+
  geom_point(size = 8)+
  geom_text_repel(data = subset(pca), 
                  aes(label = gsub("VCAP", "VCaP", 
                                   gsub("NCIH660", "NCI-H660", 
                                        gsub("MDAPCA2B", "MDA-PCa-2b", 
                                             gsub("LNCAPCLONEFGC", "LNCaP", cell.line))))),
                  nudge_y = 0.02,  
                  color = "black",
                  size = 6,  # Increase the font size
                  # Set the angle of the text (e.g., 45 degrees)
                  arrow = arrow(length = unit(0.02, "npc")))+
  theme(legend.position = 'none')

ggsave(file = '/home/liuxueying/Proj_MetPCCelline/output/plot/pca.pdf',width = 8, height = 6)


# Fig. 4b -----------------------------------------------------------------

load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/msigdb.enrich.result.RData')
ggplot(msigdb.enrich.result, aes(x = gene_ratio, y = reorder(Term, count,decreasing = F))) +
  geom_point(aes(size = count, fill = P.value), shape = 21) +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(
    
    size = "count",
    x = "Gene Ratio",
    y = "",
    title = "MsigDB Hallmark"
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 38))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color="black", linetype="solid", size = 1),
        axis.text.x = element_text( size = 14,color = 'black'),
        axis.text.y = element_text(size = 14,color = 'black'),
        legend.text = element_text(color = "black",size = 13),
        legend.title = element_text(color = "black",size = 13),
        axis.ticks.x = element_line(size=1),
        axis.ticks.y = element_line(size=1),
        plot.title = element_text(color = "black",face = "bold", size = 17,hjust = 0.5))

ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/MSigDB.pdf', height = 6.5, width = 9)



# Fig.4c-e ----------------------------------------------------------------

load( file = '/home/liuxueying/Proj_MetPCCelline/output/data/expr.marker.RData')

p1 <- ggplot(expr.marker, aes(x = reorder(cell.line, CD44, decreasing = T), y = CD44))+
  geom_col(width = 0.6, fill = "#D8CDEB")+
  scale_x_discrete(labels = c('NCIH660'='NCI-H660',
                              'VCAP'='VCaP',
                              'MDAPCA2B'='MDA-PCa-2b',
                              'LNCAPCLONEFGC'='LNCaP'))+
  theme_bw()+
  RotatedAxis()+
  ylab('CD44 expression')+xlab('')

p2 <- ggplot(expr, aes(x = reorder(cell.line, AR, decreasing = T), y = AR))+
  geom_col(width = 0.6, fill = "#D8CDEB")+
  scale_x_discrete(labels = c('NCIH660'='NCI-H660',
                              'VCAP'='VCaP',
                              'MDAPCA2B'='MDA-PCa-2b',
                              'LNCAPCLONEFGC'='LNCaP'))+
  
  
  RotatedAxis()+
  ylab('AR expression')+xlab('')


p3 <- ggplot(expr, aes(x = reorder(cell.line, SYP, decreasing = T), y = SYP))+
  geom_col(width = 0.6, fill = "#D8CDEB")+
  scale_x_discrete(labels = c('NCIH660'='NCI-H660',
                              'VCAP'='VCaP',
                              'MDAPCA2B'='MDA-PCa-2b',
                              'LNCAPCLONEFGC'='LNCaP'))+
  
  RotatedAxis()+
  ylab('SYP expression')+xlab('')

pdf(file = '/home/liuxueying/Proj_MetPCCelline/output/plot/CD44.AR.SYP.pdf', width = 15, height = 5)
p1|p2|p3
dev.off()



# Fig.4f ------------------------------------------------------------------

load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/crpc.seob.RData')
CRPC.tumor.seob <- subset(crpc.seob, tumor_subtype != 'normal')


p1 <- FeaturePlot(CRPC.tumor.seob, features = c('rna_AR'))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('AR')+
  NoLegend()
p2 <- FeaturePlot(CRPC.tumor.seob, features = c('SYP'))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('SYP')+
  NoLegend()




p3 <- FeaturePlot(CRPC.tumor.seob, features = c('rna_CD44'))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('CD44')+
  NoLegend()

p4 <- DimPlot(CRPC.tumor.seob, group.by = 'tumor_subtype',label = T, 
              cols = c('MSPC'='#D65813', 'ARPC' = '#2CAD3F', 'NEPC'='#6792CD'))+
  ggtitle('')+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('Subtype')+
  NoLegend()

pdf('/home/liuxueying/Proj_MetPCCelline/output/plot/marker.expression.pdf',width = 20, height = 5)
p4|p1|p2|p3
dev.off()


# Fig.4g -----------------------------------------------------------------

load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/mspc.cor.df.RData')

ggplot(mspc.cor.df, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 3) +  
  labs(
    
    x = "Rank",
    y = "Transcriptome similarity"
  ) +
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +
  

  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),    
    axis.text = element_text(size = 25, color = 'black'),      
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), 
    axis.ticks.length = unit(0.4, "cm"),        
    
  )+ggtitle('MSPC')
ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/MSPC.pdf', width = 9, height = 8)


# Fig.4h ------------------------------------------------------------------
load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/lineage.cor.res.RData')

ggplot(lineage.cor.res, aes(x = variable, y = value)) +
  geom_violin(aes(fill = variable, colour = variable), alpha = 0.5) +
  # alpha控制不透明度
  geom_boxplot(aes(colour = variable), width = 0.2)+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )+
  scale_color_manual(values = c("BE" = '#D65813', "LE" = '#2CAD3F', 'NE'='#6792CD'))+
  scale_fill_manual(values = c("BE" = '#D65813', "LE" = '#2CAD3F', 'NE'='#6792CD'))+
  theme_niwot()+
  stat_compare_means(paired = T, comparisons = list(c('BE','LE')))+
  
  scale_x_discrete(labels = c('BE' = paste('Basal', sep = '\n', '(n = 2785)'),
                              'LE' = paste('Luminal', sep = '\n', '(n = 2785)'),
                              'NE' = paste('Neuroendocrine', sep = '\n', '(n = 2785)')))+
  
  xlab('')+
  ylab('Correlation')

ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/lineage.cor.res.pdf', width = 12, height = 6)


# Fig.5a ------------------------------------------------------------------

load(file = '/home/liuxueying/Proj_MetPCCelline/output/pca.atac.RData')

ggplot(pca.atac,aes(x= PC1, y = PC2, color = cell.line))+
  geom_point(size = 8)+
  geom_text_repel(data = subset(pca), 
                  aes(label = gsub("VCAP", "VCaP", 
                                   gsub("NCIH660", "NCI-H660", 
                                        gsub("MDAPCA2B", "MDA-PCa-2b", 
                                             gsub("LNCAPCLONEFGC", "LNCaP", cell.line))))),
                  nudge_y = 0.02,  
                  color = "black",
                  size = 6,  # Increase the font size
                  # Set the angle of the text (e.g., 45 degrees)
                  arrow = arrow(length = unit(0.02, "npc")))+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )+
  theme(legend.position = 'none')

ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/atac.pca.pdf', width = 10, height = 10)




# Fig. 5b -----------------------------------------------------------------

load(file = '/home/liuxueying/met_pc_cell_line/pcdata/paper/data/section55/feature_annotation.RData')

### FUNCTION: getGRanges() ---
getGRanges <- function(dat){
  dat <- subset(dat, select=c("seqnames","start","end"))
  gr <- GenomicRanges::makeGRangesFromDataFrame(dat, keep.extra.columns=FALSE)
  return(gr)
}


### FUNCTION: processEnrichment() ---
processEnrichment <- function(dat, pval){
  dat <- subset(dat, select=c("id","fold_enrichment","p_value","p_adjust"))
  dat <- dat[which(dat$p_adjust <= pval),]
  dat <- dat[order(dat$p_adjust, decreasing=FALSE),]
  dat$nlogp <- -log10(dat$p_adjust)
  
  return(dat)
}

### FUNCTION: plotBar() ----
plotBar <- function(dat, analysis_id){
  # PREPARE DATA ---
  dat <- dat[order(dat$nlogp, decreasing=FALSE),]
  dat$id <- factor(dat$id, levels=dat$id)
  
  # PLOT ---
  p <- ggplot(dat, aes(x=nlogp, y=id)) +
    geom_bar(stat="identity", fill="#006094", width=0.8) +
    coord_cartesian(xlim=c(0,12)) +
    scale_x_continuous(breaks=seq(0,12,by=4)) +
    theme(
      axis.text.x = element_text(size = 11, color="#000000"),
      axis.text.y = element_text(size = 11, color="#000000", ),
      axis.title = element_text(size = 11, color="#000000"),
      plot.title = element_text(size = 11, color="#000000", hjust=0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(size=0.2, color="#000000"), 
      strip.text = element_text(size=11, color="#000000"),
      strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
      panel.background = element_rect(fill="#FFFFFF", color="#000000"),
      legend.text = element_text(size = 11, color="#000000"),
      legend.title = element_blank(),
      legend.key.size = unit(0.3, "cm"),
      legend.position = "none") +
    ylab("") +
    xlab("-log10(FDR)") + 
    ggtitle(analysis_id) 
  
  return(p)
}

gr_all             <- getGRanges(dat = feature_annotation)

### SET SEED ---
set.seed(12345)

### GREAT ENRICHMENT ---

obj_all        <- rGREAT::great(gr=gr_all, gene_sets="MSigDB:H", tss_source="txdb:hg38", 
                                biomart_dataset = NULL, min_gene_set_size = 5, 
                                mode = "basalPlusExt", basal_upstream = 5000, basal_downstream = 1000, extension = 1000000, 
                                extended_tss = NULL, background = NULL, exclude = "gap", cores = 50, verbose = TRUE)
### PARSE ENRICHMENT DATA ---

dat_all        <- obj_all@table

df_all         <- processEnrichment(dat=dat_all, pval=0.05)

p1          <- plotBar(dat=df_all, analysis_id="MSigDB Hallmark")

pdf('/home/liuxueying/Proj_MetPCCelline/output/plot/enrich.pdf',width = 8 ,height = 6)
p1
dev.off()


# Fig.5b ------------------------------------------------------------------

load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/atac.organoid.RData')

cor.stem.res             <- cor(atac.three.subtype.mcrpc[rownames(mat_features),stem.id], atac.organoid[rownames(mat_features),], method = 'spearman')
cor.stem.res %>% reshape2::melt() %>% 
  ggplot(aes(x = reorder(Var2,value, median), y = value , fill = Var2, color =Var2))+
  geom_boxplot(alpha = 0.5)+
  geom_jitter(
    aes(color = Var2),              
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
    size = 1.2,
  )+
  stat_compare_means(comparisons = list(c('PC3','MSK-PCa12')
  ))+xlab('')+ylab('Correlation')+theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )+RotatedAxis()
ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/organoid.atac.pdf', width = 12, height = 6)



# Fig.5c-e-----------------------------


load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/patient.samples.RData')

cor.AR.res                <- cor(atac.three.subtype.mcrpc[rownames(mat_features),AR.id], atac.count[rownames(mat_features),], method = 'spearman')
boxplot(cor.AR.res)
cor.AR.res %>% reshape2::melt() %>% 
  ggplot(aes(x = reorder(Var2,value, median), y = value, fill = Var2, colour = Var2))+
  
  geom_boxplot(alpha = 0.5)+
  scale_fill_manual(values = my.cols)+
  scale_colour_manual(values = my.cols)+
  geom_jitter(
    aes(color = Var2),              
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
    size = 1.2,
   
   
  )+
  stat_compare_means(comparisons = list(
    c('VCaP','LNCaP')), method = 'wilcox.test', paired = T)+
  xlab('')+
  ylab('')+
  ylab('Correlation')+
  theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 35),
    axis.title = element_text(size = 35, hjust = 0.5, color = 'black'),
    axis.text = element_text(size = 35, color = 'black'),
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA),
    axis.ticks.length = unit(0.6, "cm")
    
  )+
  theme(legend.position = "none")+
  ggtitle('ARPC')
ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/ARPC.atac.pdf', width = 15, height = 8)

cor.NE.res                <- cor(atac.three.subtype.mcrpc[rownames(mat_features),neuroendocrine.id], atac.count[rownames(mat_features),], method = 'spearman')

cor.NE.res %>% reshape2::melt() %>% 
  ggplot(aes(x = reorder(Var2,value, median), y = value, fill = Var2, colour = Var2))+
  geom_boxplot(alpha = 0.5)+
  scale_fill_manual(values = my.cols)+
  scale_colour_manual(values = my.cols)+
  geom_jitter(
    aes(color = Var2),                
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
    size = 1.2,
  )+
  stat_compare_means(comparisons = list(c('NCI-H660','DU145'),c('NCI-H660','PC3')), method = 'wilcox.test',paired = T)+
  scale_x_discrete(labels = c('NCI-H660' = 'NCI-H660',
                              '22RV1' = '22RV1',
                              'DU145' = 'DU145',
                              'LNCaP'='LNCaP',
                              
                              'PC3'='PC3',
                              'VCaP'='VCaP')
  )+
  ylab('Correlation')+
  theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 35),
    axis.title = element_text(size = 35, hjust = 0.5, color = 'black'),
    axis.text = element_text(size = 35, color = 'black'),
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA),
    axis.ticks.length = unit(0.6, "cm")
  )+
  theme(legend.position = "none")+
  xlab('')+
  
  ggtitle('NEPC')
ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/NEPC.atac.pdf', width = 15, height = 8)




cor.stem.res              <- cor(atac.three.subtype.mcrpc[rownames(mat_features),stem.id], atac.count[rownames(mat_features),], method = 'spearman')
cor.stem.res %>% reshape2::melt() %>% 
  ggplot(aes(x = reorder(Var2,value, median), y = value, fill = Var2, colour = Var2))+
  geom_boxplot(alpha = 0.5)+
  scale_fill_manual(values = my.cols)+
  scale_colour_manual(values = my.cols)+
  geom_jitter(
    aes(color = Var2),                
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
    size = 1.2,
  )+
  stat_compare_means(comparisons = list(c('NCI-H660','DU145'),c('NCI-H660','PC3')), method = 'wilcox.test',paired = T)+
  ylab('Correlation')+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )+
  theme(legend.position = "none")+
  
  xlab('')+
  ggtitle('MSPC')
ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/MSPC.atac.pdf', width = 15, height = 8)









