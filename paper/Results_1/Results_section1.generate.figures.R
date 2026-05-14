setwd('/home/liuxueying/met_pc_cell_line/')
# Fig.1a ---------------------------------------------------------------------
load(file = './data/CCLE.SU2C.TCGA.mutation.matrix.RData')
load(file = './data/compare.somatic.mutation.RData')
full.gene               <- c(MET.highly.mutated.gene,MET.dm.gene) %>% unique  

CCLE.mutation.fequency  <- apply(CCLE.mutation.matrix[full.gene,], 1, sum) / ncol(CCLE.mutation.matrix)



gene.order           <- CCLE.mutation.fequency %>% sort(decreasing = T) %>% names()
TCGA.sample.order    <- apply(TCGA.mutation.matrix[full.gene,], 2, sum) %>% sort(decreasing = T) %>% names()
TCGA.mutation.matrix <- TCGA.mutation.matrix[, TCGA.sample.order]
SU2C.sample.order    <- apply(SU2C.mutation.matrix[full.gene,], 2, sum) %>% sort(decreasing = T) %>% names()
SU2C.mutation.matrix <- SU2C.mutation.matrix[, SU2C.sample.order]
CCLE.sample.order    <- apply(CCLE.mutation.matrix[full.gene,], 2, sum) %>% sort(decreasing = T) %>% names()
CCLE.mutation.matrix <- CCLE.mutation.matrix[, CCLE.sample.order]


mutation.profile <- cbind(TCGA.mutation.matrix[gene.order,],SU2C.mutation.matrix[gene.order,])


col.df.pro <- c( rep(times= colnames(TCGA.mutation.matrix)  %>% length, x='TCGA'),
                 rep(times= colnames(SU2C.mutation.matrix)   %>% length, x='SU2C')
)

names(col.df.pro)          <- c(colnames(TCGA.mutation.matrix), colnames(SU2C.mutation.matrix))
col.df.pro                 <- as.data.frame(col.df.pro)
colnames(col.df.pro)       <- 'type'

col.annotation    <-  HeatmapAnnotation(type=col.df.pro$type, 
                                        col = list(type = c("SU2C" =  "#4DAE48", "TCGA" = "#974F9F")) ,
                                        show_annotation_name = F,## 不显示type
                                        show_legend = F
)


# 定义 alter_fun
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h, gp = gpar(fill = '#F0F0F0', col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.8, gp = gpar(fill = "black", col = NA))
  }
)


col <- c("MUT" = "black", "background" = '#F0F0F0')
sample_order <- rownames(col.df.pro)
identical(colnames(mutation.profile), rownames(col.df.pro))

p1 <- oncoPrint(
  mutation.profile,
  get_type = function(x) {ifelse(x == 1, "MUT", "background")},
  row_names_side = "left",          
  show_heatmap_legend = FALSE,      
  alter_fun = alter_fun,           
  col = col,                       
  show_column_names = FALSE,        
  show_row_names = F, # 不显示行名
  top_annotation = col.annotation,  
  right_annotation = NULL,         
  show_pct = FALSE,                
  row_names_gp = gpar(fontsize = 40, fontface = 'bold'),
  column_order = sample_order,
  row_order = gene.order)




###### CCLE


mutation.profile <- CCLE.mutation.matrix[full.gene,]
sample_order     <- mutation.profile %>% colSums() %>% sort(decreasing = T) %>% names()
mutation.profile <- mutation.profile[gene.order,sample_order]


col.df.pro <- c( rep(times= colnames(mutation.profile)  %>% length, x='CCLE')
)


names(col.df.pro)          <- c(colnames(mutation.profile))
col.df.pro                 <- as.data.frame(col.df.pro)
colnames(col.df.pro)       <- 'type'

col.annotation    <-  HeatmapAnnotation(type=col.df.pro$type, 
                                        col = list(type = c("CCLE" =  "#377EB9")) ,
                                        show_annotation_name = F,## 不显示type
                                        show_legend = F
)



alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h, gp = gpar(fill = '#F0F0F0', col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.8, gp = gpar(fill = "black", col = NA))
  }
)
'#D3D3D3'
'#DCDCDC'

'#F0F0F0'


'#F5F5F5'


col            <- c("MUT" = "black", "background" = '#F0F0F0')
sample_order   <- rownames(col.df.pro)
identical(colnames(mutation.profile), rownames(col.df.pro))
# 绘制 OncoPrint
p2 <- oncoPrint(
  mutation.profile,
  get_type = function(x) {ifelse(x == 1, "MUT", "background")},
  #row_names_side = "right",          
  show_heatmap_legend = FALSE,     
  alter_fun = alter_fun,           
  col = col,                        
  # show_column_names = T,        
  top_annotation = col.annotation,  
  right_annotation = NULL,         
  show_pct = FALSE,                 
  row_names_gp = gpar(fontsize = 40, fontface = 'bold'),
  show_row_names = F,
  
  row_order = gene.order,
  column_order = sample_order,
)



p3 <- grid.grabExpr({
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2, widths = unit(c(9, 2), "null"))))
  
  pushViewport(viewport(layout.pos.col = 1))
  draw(p1, newpage = FALSE)
  upViewport()
  
  pushViewport(viewport(layout.pos.col = 2))
  draw(p2, newpage = FALSE)
  upViewport(2)
})
## 画图
pdf('aa.pdf', width = 90, height = 40)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2, widths = unit(c(9, 2), "null"))))

pushViewport(viewport(layout.pos.col = 1))
draw(p1, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.col = 2))
draw(p2, newpage = FALSE)
upViewport(2)
dev.off()

#CCLE  mutation frequency 
CCLE.mutation.fequency  <- apply(CCLE.mutation.matrix[full.gene,], 1, sum) / ncol(CCLE.mutation.matrix)
mutation.frequency      <- CCLE.mutation.fequency
mutation.frequency      <- as.data.frame(mutation.frequency)
colnames(mutation.frequency) <- 'CCLE.mutation.fequency'
col.df                  <- c('CCLE')
gene.order              <- CCLE.mutation.fequency %>% sort(decreasing = T) %>% names()
names(col.df)           <- colnames(mutation.frequency)
col.df                  <- as.data.frame(col.df)
colnames(col.df)        <- 'type'
col.ha <- HeatmapAnnotation(type = col.df$type,
                            col = list(type = c("CCLE" = "#377EB9")),
                            show_legend = T,
                            show_annotation_name = F,### rm type 
                            annotation_name_gp = gpar(fontface = "bold"),
                            annotation_legend_param = list(
                              title_gp = gpar(fontsize = 20, fontface = 'bold'), ## 图例字体
                              labels_gp = gpar(fontsize = 20),
                              legend_height = unit(10, "cm"), ###调整col annotation 图例大小
                              legend_width = unit(5, "cm")
                            ))
mutation.frequency <- as.matrix(mutation.frequency)
p4 <- Heatmap(mutation.frequency,
              show_column_dend = FALSE, 
              col=colorRamp2(c(0,1),c('blue','red')),
              row_names_gp = gpar(fontsize = 40,fontface='bold'),
              width = unit(1, "cm"),
              show_heatmap_legend = T,
              top_annotation=col.ha,
              show_column_names = F,
              cluster_columns = F,
              cluster_rows = F,
              column_order = c("CCLE.mutation.fequency"),
              row_order = gene.order,
              name = 'Mutation \nfrequency',
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 20, fontface = "bold"),
                labels_gp = gpar(fontsize = 20),                  
                legend_height = unit(10, "cm"),                    
                legend_width = unit(5, "cm")                      
              ))


grob_p4 <- grid.grabExpr({
  draw(p4, newpage = FALSE)
})



pdf("./plot/Fig.1a.pdf", width = 90, height = 40)


grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2, widths = unit(c(4, 0.2), "null"))))


pushViewport(viewport(layout.pos.col = 1))
grid.draw(p3)
upViewport()


pushViewport(viewport(layout.pos.col = 2))
grid.draw(grob_p4)
upViewport(2)
dev.off()

#Fig. 1b Number of mutated genes -------------------------------------------------
load(file = './data/compare.somatic.mutation.RData')
load(file= './data/CCLE.mutation.burden.RData')
full.gene                                      <- c(MET.highly.mutated.gene,MET.dm.gene) %>% unique 
resemble.gene.cnt                              <- apply(CCLE.mutation.matrix[full.gene,],2,sum)
resemble.gene.mutation.burden <- resemble.gene.cnt %>% sort() %>% as.data.frame()
colnames(resemble.gene.mutation.burden)        <- 'mutation_burden'
resemble.gene.mutation.burden$rank             <- rank(resemble.gene.mutation.burden$mutation_burden)
resemble.gene.mutation.burden$cell.line        <- rownames(resemble.gene.mutation.burden)

df                                             <- resemble.gene.mutation.burden

df$name <- gsub(pattern = '_PROSTATE', '', df$cell.line)
df$name <- gsub(df$name, pattern = 'LNCAPCLONEFGC', replacement = 'LNCaP')
df$name <- gsub(df$name, pattern = 'MDAPCA2B', replacement = 'MDA-PCa-2b')
df$name <- gsub(df$name, pattern = 'NCIH660', replacement = 'NCI-H660')
df$name <- gsub(df$name, pattern = 'VCAP', replacement = 'VCaP')

#write.csv(df, file = "/home/liuxueying/met_pc_cell_line/pcdata/paper/data/reproduce/section1/Fig.1b.csv", row.names = TRUE)


df$name <- factor(df$name, levels = df$name[order(df$rank)])


ggplot(df, aes(x = name, y = mutation_burden)) +
  geom_col(width = 0.6, fill = "#D8CDEB") +
  labs(y = "Number of mutated genes") +
  xlab('')+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.title = element_text(size = 15, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 15, color = 'black'),      # 坐标轴刻度加粗
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.2, "cm"),         # 增加刻度线长度
    
  )
ggsave('./plot/Number.of.mutated.genes.pdf', width = 5, height = 5)

#Fig. 1c CCLE mutation burden ----------------------------------------------------

load(file= './pcdata/paper/data/section1/CCLE.mutation.burden.RData')
CCLE.df      <- data.frame(name = names(CCLE.mutation.burden), mutation.burden = log2(CCLE.mutation.burden))


df                                             <- CCLE.df

df$name <- gsub(pattern = '_PROSTATE', '', df$name)
df$name <- gsub(df$name, pattern = 'LNCAPCLONEFGC', replacement = 'LNCaP')
df$name <- gsub(df$name, pattern = 'MDAPCA2B', replacement = 'MDA-PCa-2b')
df$name <- gsub(df$name, pattern = 'NCIH660', replacement = 'NCI-H660')
df$name <- gsub(df$name, pattern = 'VCAP', replacement = 'VCaP')
df$rank <- rank(df$mutation.burden)
df      <- df[order(df$rank),]
ggplot(df,aes(x= name, y = mutation.burden))+
  geom_col(width = 0.6, fill = "#D8CDEB")+
  xlab('Rank')+
  ylab('Mutation burden')+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.title = element_text(size = 15, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 15, color = 'black'),      # 坐标轴刻度加粗
    
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.2, "cm"),         # 增加刻度线长度
    
  )

ggsave('./plot/section1/CCLE.mutation.burden.pdf', width = 5.5, height = 5)
#write.csv(df, file = "/home/liuxueying/met_pc_cell_line/pcdata/paper/data/reproduce/section1/Fig.1c.csv", row.names = TRUE)


# Fig.1d ------------------------------------------------------------------


load(file= './data/CCLE.mutation.burden.RData')
resemble.gene.cnt                              <- apply(CCLE.mutation.matrix[full.gene,],2,sum)
( resemble.gene.cnt[names(CCLE.mutation.burden)] / CCLE.mutation.burden ) %>% sort
resemble.gene.cnt                       <- apply(CCLE.mutation.matrix[full.gene,],2,sum)
aa                                      <- ( resemble.gene.cnt[names(CCLE.mutation.burden)] / CCLE.mutation.burden )
aa                                      <- as.data.frame(aa)
aa
aa                                      <- aa[order(aa$aa),, drop=F]
aa$rank <- rank(aa$aa)
aa$name <-  gsub(rownames(aa), pattern = '_PROSTATE',replacement = '')
aa$name <- gsub(aa$name, pattern = 'LNCAPCLONEFGC', replacement = 'LNCaP')
aa$name <- gsub(aa$name, pattern = 'MDAPCA2B', replacement = 'MDA-PCa-2b')
aa$name <- gsub(aa$name, pattern = 'NCIH660', replacement = 'NCI-H660')
aa$name <- gsub(aa$name, pattern = 'VCAP', replacement = 'VCaP')

aa$name <- factor(aa$name, levels = aa$name[order(aa$rank)])


p1 <- ggplot(aa,aes(x= name, y = aa)) + 
  geom_col(width = 0.6, fill = "#D8CDEB")+
  ylab('Ratio')+
  xlab('')+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.title = element_text(size = 15, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 15, color = 'black'),      # 坐标轴刻度加粗
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.2, "cm"),         # 增加刻度线长度
    
  )
p1
ggsave('./plot/section1/ratio.pdf', width = 5, height = 5)
#write.csv(aa, file = "/home/liuxueying/met_pc_cell_line/pcdata/paper/data/reproduce/section1/Fig.1d.csv", row.names = TRUE)

#Fig.1e hotspot maf -------------------------------------------------------------
load(file = './data/supplementary.table1.RData')


hotspot.maf = read.maf(maf = supplementary.table1)


pdf('./plot/AR.pdf', width = 15, height = 10)
lollipopPlot(
  maf = hotspot.maf,
  gene = 'AR',
  AACol = 'HGVSp_Short',
  showMutationRate = TRUE,
  labelPos =c(702,742, 875,878)
  
)
dev.off()


pdf('./plot/TP53.pdf', width = 15, height = 10)
lollipopPlot(
  maf = hotspot.maf,
  gene = 'TP53',
  AACol = 'HGVSp_Short',
  showMutationRate = TRUE,
  labelPos = c(282,273,248,245,220,179,175)
  
)
dev.off()




#write.csv(supplementary.table1, file = "/home/liuxueying/met_pc_cell_line/pcdata/paper/data/reproduce/section1/Fig.1e.csv", row.names = TRUE)







#Fig.1f hot.mutation.code.ccle --------------------------------------------------
load('./data/organize.mutation.data.RData')
#Next, let us check hotspot mutations
non.silent.PC.SU2C.maf.data  <- PC.maf.data.list$SU2C %>% filter(Variant_Classification != 'Silent')  #### 转移癌样本体细胞突变数据  
PC.SU2C.hotspot.mutation     <- table (non.silent.PC.SU2C.maf.data$mutation.code) %>% as.data.frame %>% filter(Freq >= 3)  %>%   dplyr::select(Var1) %>% unlist %>% as.character()  #####转移癌中高突变的位置
df                           <- filter(PC.maf.data.list$CCLE, mutation.code %in% PC.SU2C.hotspot.mutation) %>% as.data.frame  ###### 转移癌中高突变的位置在细胞系中出现在哪些细胞系里面
View(df) #"1@6257785@RPL22" resembled by three hyper-mutated cell lines; this mutation is usually seen in hyper-mutated samples, another piece of proof! 
#Ppaer:High Frequency of RPL22 Mutations in Microsatellite-Unstable Colorectal and Endometrial Tumors
#Quantitative Proteomics of the Cancer Cell Line Encyclopedia
df <- df[,c(38,39)]
colnames(df)[1] <- 'CELL_LINE_NAME'
cell_lines <- PC.maf.data.list$CCLE$CCLE_Name %>% table() %>% names()

df$CELL_LINE_NAME <- factor(df$CELL_LINE_NAME, levels = cell_lines)

hot.mutation.code.ccle <- table(df$mutation.code, df$CELL_LINE_NAME)

hot.mutation.code.ccle <- ifelse(hot.mutation.code.ccle > 0, 1, 0)

hot.mutation.code.ccle

original_rownames <- rownames(hot.mutation.code.ccle)

split_names <- strsplit(original_rownames, "@")

new_rownames <- sapply(split_names, function(x){
  paste0(x[3], " (chr", x[1], ":", x[2], ")")
})

rownames(hot.mutation.code.ccle) <- new_rownames

hot.mutation.code.ccle
colnames(hot.mutation.code.ccle) <- gsub(colnames(hot.mutation.code.ccle), pattern = '_PROSTATE',replacement = '')
colnames(hot.mutation.code.ccle)[5] <- 'MDA-PCa-2b'
colnames(hot.mutation.code.ccle)[6] <- 'NCI-H660'
colnames(hot.mutation.code.ccle)[4] <- 'LNCaP'
colnames(hot.mutation.code.ccle)[12] <- 'VCaP'
sample.order <- apply(hot.mutation.code.ccle, 2, sum) %>% sort() %>% names()
all_samples  <- colnames(hot.mutation.code.ccle)


p1 <- oncoPrint(
  hot.mutation.code.ccle,
  get_type = function(x) {ifelse(x == 1, "MUT", "background")},
  row_names_side = "left",         
  show_heatmap_legend = F,     
  alter_fun = alter_fun,           
  col = col,                        
  show_column_names = T,        
  top_annotation = NULL,        
  right_annotation = NULL,          
  show_pct = FALSE,                 
  row_names_gp = gpar(fontsize = 20),    
  column_names_gp = gpar(fontsize = 25), 
  column_order = sample.order, 
  column_names_rot = 45)        
pdf('./plot/hot.mutation.code.ccle.pdf', width = 12, height = 9)
p1
dev.off()


#Fig. 1g all samples dm cnv ------------------------------------------------------

load(file = './data/TCGA.gene.cnv.matrix.RData')
load(file = './data/SU2C.gene.cnv.matrix.RData')


c.gene <- intersect(rownames(TCGA.gene.cnv.matrix), rownames(SU2C.gene.cnv.matrix))
rs.df  <- foreach(g = c.gene,.combine = 'rbind') %do% {
  p.value <- wilcox.test(TCGA.gene.cnv.matrix[g,],SU2C.gene.cnv.matrix[g,])$p.value
  delta   <- SU2C.gene.cnv.matrix[g, ] %>% median - TCGA.gene.cnv.matrix[g,] %>% median
  data.frame(delta = delta, p.value = p.value,gene = g)
}  
rownames(rs.df)       <- rs.df$gene
rs.df$p.adj           <- p.adjust(rs.df$p.value,method = 'bonferroni')


rs.df$sig <- ifelse(rs.df$p.adj<0.05 & rs.df$delta > 0.3, 'up',
                    ifelse(rs.df$p.adj <0.05 & rs.df$delta < -0.3, 'dn', 'no'))


ggplot(rs.df, aes(x = delta, y = -log10(p.adj),color=sig))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-0.3,0.3), linetype = "dashed", color = "black")+
  
  scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
  scale_size_continuous(range = c(0,1))+
  theme(legend.title = element_blank() )+
  geom_label_repel(data = subset(rs.df, gene %in% 'AR'),
                   max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                   aes(label = gene),
                   size = 5, 
                   color = 'black',fontface = "bold.italic",
                   segment.color = "black",   
                   segment.linetype = "solid", 
                   segment.size = 0.5 ,
                   box.padding = 3,        
                   point.padding = 0.2 ) +
  xlab("Difference")+
  ylab("-log10 p-adj")+
  xlim(-0.9, 0.9) +
  ggtitle('')+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.title = element_text(size = 15, hjust = 0.5, color = 'black'),     
    axis.text = element_text(size = 15, color = 'black'),      
    axis.text.x = element_text( hjust = 1, vjust = 1),
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), 
    axis.ticks.length = unit(0.2, "cm"),        
    
  )


ggsave('./pcdata/paper/plot/section1/all.dm.CNV.pdf', width = 5, height = 5)
#write.csv(rs.df, file = "/home/liuxueying/met_pc_cell_line/pcdata/paper/data/reproduce/section1/Fig.1g.csv", row.names = TRUE)




#Fig.1h AR CNV  -----------------------------------------------------------------
load(file = './data/CCLE.PC.gene.cnv.matrix.RData')


df <- CCLE.PC.gene.cnv.matrix['AR', ] %>% sort() %>% as.data.frame()
colnames(df) <- 'CNV'
df$Rank <- rank(df$CNV)
df$name <- gsub(pattern = '_PROSTATE', x = rownames(df), replacement = '')

df$name <- gsub(df$name, pattern = 'LNCAPCLONEFGC', replacement = 'LNCaP')
df$name <- gsub(df$name, pattern = 'MDAPCA2B', replacement = 'MDA-PCa-2b')
df$name <- gsub(df$name, pattern = 'NCIH660', replacement = 'NCI-H660')
df$name <- gsub(df$name, pattern = 'VCAP', replacement = 'VCaP')

df$name <- factor(df$name, levels = df$name[order(df$Rank)])
ggplot(df, aes(x = name, y = CNV)) +
  geom_col(width = 0.6, fill = "#D8CDEB") +
  labs(y = "AR CNV") +
  xlab('')+
  ggplot.style
ggsave('./plot/AR.CNV.pdf', width = 4, height = 4)

#write.csv(df, file = "/home/liuxueying/met_pc_cell_line/pcdata/paper/data/reproduce/section1/Fig.1h.csv", row.names = TRUE)
