library(data.table)
library(Seurat)
library(dplyr)
library(doParallel)
library(iterators)
library(foreach)
library(DESeq2)

source('/home/liuxueying/met_pc_cell_line/code/BioKLab.util.R')
load('/home/liuxueying/met_pc_cell_line/processed.data/CCLE.transcriptome.RData')
hgnc.data                    <- read.delim("./raw.data/HGNC/hgnc_complete_set.txt", stringsAsFactors=FALSE)

flag                         <- grepl(colnames(CCLE.log2.tpm.matrix), pattern ='PROSTATE')
pc.cell.line.log2.tpm.matrix <- CCLE.log2.tpm.matrix[,flag]
flag                         <- colnames(pc.cell.line.log2.tpm.matrix) != 'PRECLH_PROSTATE' ## remove no cancer 
pc.cell.line.log2.tpm.matrix <- pc.cell.line.log2.tpm.matrix[,flag]

pca.rs                       <- prcomp(pc.cell.line.log2.tpm.matrix %>% t)    
plot(pca.rs$x[,1:2])

pc1.tail.gene                <- pca.rs$rotation[,'PC1'] %>% sort %>% tail(200) %>% names
pc1.head.gene                <- pca.rs$rotation[,'PC1'] %>% sort %>% head(200) %>% names


# pc1 tail 200 gene enrichment analysis --------------------------------------------------
`MSigDB_Hallmark_2020_table.(2)` <- read.delim("/home/liuxueying/met_pc_cell_line/pcdata/explore_pc3/input/MSigDB_Hallmark_2020_table (2).txt")
msigdb.enrich.result <- `MSigDB_Hallmark_2020_table.(2)`

msigdb.enrich.result <- msigdb.enrich.result %>% 
  separate(col = 'Overlap',into = c('count','gene_set_count'),sep = '/',remove = F)

msigdb.enrich.result$count <- as.numeric(msigdb.enrich.result$count)
msigdb.enrich.result <- msigdb.enrich.result %>% mutate(gene_ratio=msigdb.enrich.result$count/200)


msigdb.enrich.result <- msigdb.enrich.result[order(-msigdb.enrich.result$count), ]
msigdb.enrich.result <- msigdb.enrich.result %>% filter(Adjusted.P.value < 0.05)


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
  theme_bw()
  
  

# AR CD44 SYP expression --------------------------------------------------
pc.cell.line.log2.tpm.matrix[symbol.to.ensemble('CD44'),] %>% sort()
pc.cell.line.log2.tpm.matrix[symbol.to.ensemble('AR'),] %>% sort()
pc.cell.line.log2.tpm.matrix[symbol.to.ensemble('SYP'),] %>% sort()

##scRNA pc3 validation
pc3.seob                    <- readRDS("/home/liuxueying/met_pc_cell_line/pcdata/met_pc/pc3.seob.rds")

p1 <- FeaturePlot(pc3.seob, features = c('AR'))
p2 <- FeaturePlot(pc3.seob, features = c('SYP'))
p3 <- FeaturePlot(pc3.seob, features = c('CD44'))


# CRPC expression ---------------------------------------------------------
load(file = '/home/liuxueying/met_pc_cell_line/pcdata/paper/data/section3/CRPC.confident.seob.RData')
CRPC.tumor.seob <- subset(CRPC.confident.seob, tumor_subtype != 'normal')


p1 <- FeaturePlot(CRPC.tumor.seob, features = c('AR'))
p2 <- FeaturePlot(CRPC.tumor.seob, features = c('SYP'))
p3 <- FeaturePlot(CRPC.tumor.seob, features = c('CD44'))



#  MSPC  TC analysis--------------------------------------------------------------


####Fig 4 MSPC 

load(file = './pcdata/paper/data/section3/CRPC.confident.seob.RData')
CRPC.tumor.seob <- subset(CRPC.confident.seob, tumor_subtype != 'normal')


MSPC_seob <- subset(CRPC.tumor.seob, tumor_subtype == 'MSPC')


MSPC.cor.res <- TC.rs.list(MSPC_seob)

df <- MSPC.cor.res$median.cor
colnames(df)[1] <- 'correlation'


df$name      <- rownames(df)
df$prostate  <- ifelse(rownames(df) %in% pc.cell.line.name, 
                       "prostate", 
                       "others")
colnames(df)[1] <- "correlation"

df$rank     <- rank(df$correlation)
df_ranked   <- df[order(df$rank), ]

#write.xlsx(df_ranked, file = 'mspc.tc.xlsx')
#--- B. plot ---
highlight_points <- df_ranked[df_ranked$rank %in% c(1018, 1019,1017), ]
ggplot(df_ranked, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 3) +  
  labs(
    # 标题使用该列名，或你想替换的名字
    x = "Rank",
    y = "Transcriptome similarity"
  ) +
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +
  # 在图中高亮 pc.cell.line.name 这几行
  geom_point(
    data = subset(df, name %in% pc.cell.line.name),
    color = "red",
    size = 4
  )+
  geom_text_repel(data = subset(df, name %in% pc.cell.line.name), 
                  aes(label = gsub(x= name, pattern = '_PROSTATE','')), 
                  nudge_y = 0.02,  # 标签在纵坐标上的偏移量
                  color = "black",
                  size = 4,
                  arrow = arrow(length = unit(0.02, "npc")))+
  
 
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )+ggtitle('MSPC')

ggsave('./pcdata/paper/plot/section4/MSPC.pdf', width = 8, height = 8)





# Fig4 ARPC TC analysis ---------------------------------------------------------------


####Fig 4 ARPC  

load(file = './pcdata/paper/data/section3/CRPC.confident.seob.RData')
CRPC.tumor.seob <- subset(CRPC.confident.seob, tumor_subtype != 'normal')


ARPC_seob <- subset(CRPC.tumor.seob, tumor_subtype == 'ARPC')


ARPC.cor.res <- TC.rs.list(ARPC_seob)

df <- ARPC.cor.res$median.cor
colnames(df)[1] <- 'correlation'


df$name      <- rownames(df)
df$prostate  <- ifelse(rownames(df) %in% pc.cell.line.name, 
                       "prostate", 
                       "others")
colnames(df)[1] <- "correlation"

df$rank     <- rank(df$correlation)
df_ranked   <- df[order(df$rank), ]

#--- B. plot ---

ggplot(df_ranked, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 3) +  
  labs(
    # 标题使用该列名，或你想替换的名字
    x = "Rank",
    y = "Transcriptome similarity"
  ) +
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +
  # 在图中高亮 pc.cell.line.name 这几行
  geom_point(
    data = subset(df, name %in% pc.cell.line.name),
    color = "red",
    size = 4
  )+
  geom_text_repel(data = subset(df, name %in% pc.cell.line.name), 
                  aes(label = gsub(x= name, pattern = '_PROSTATE','')), 
                  nudge_y = 0.02,  # 标签在纵坐标上的偏移量
                  color = "black",
                  size = 4,
                  arrow = arrow(length = unit(0.02, "npc")))+
  
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )+ggtitle('ARPC')

ggsave('./pcdata/paper/plot/section4/ARPC.pdf', width = 8, height = 8)


# Fig NEPC TC analysis----------------------------------------------------------------


####Fig 4 NEPC  

load(file = './pcdata/paper/data/section3/CRPC.confident.seob.RData')
CRPC.tumor.seob <- subset(CRPC.confident.seob, tumor_subtype != 'normal')


NEPC_seob <- subset(CRPC.tumor.seob, tumor_subtype == 'NEPC')


NEPC.cor.res <- TC.rs.list(NEPC_seob)

df <- NEPC.cor.res$median.cor
colnames(df)[1] <- 'correlation'


df$name      <- rownames(df)
df$prostate  <- ifelse(rownames(df) %in% pc.cell.line.name, 
                       "prostate", 
                       "others")
colnames(df)[1] <- "correlation"

df$rank     <- rank(df$correlation)
df_ranked   <- df[order(df$rank), ]

#--- B. plot ---

ggplot(df_ranked, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 3) +  
  labs(
    # 标题使用该列名，或你想替换的名字
    x = "Rank",
    y = "Transcriptome similarity"
  ) +
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +
  # 在图中高亮 pc.cell.line.name 这几行
  geom_point(
    data = subset(df, name %in% pc.cell.line.name),
    color = "red",
    size = 4
  )+
  geom_text_repel(data = subset(df, name %in% pc.cell.line.name), 
                  aes(label = gsub(x= name, pattern = '_PROSTATE','')), 
                  nudge_y = 0.02,  # 标签在纵坐标上的偏移量
                  color = "black",
                  size = 4,
                  arrow = arrow(length = unit(0.02, "npc")))+
  
  
  
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )+ggtitle('NEPC')

ggsave('./pcdata/paper/plot/section4/NEPC.pdf', width = 8, height = 8)



# ARPC-vcap ---------------------------------------------------------------
load(file = '/home/liuxueying/met_pc_cell_line/pcdata/paper/data/section1/CCLE.log2.tpm.matrix.symbol.RData')
load("/home/liuxueying/met_pc_cell_line/pcdata/CCLE.rna.seq.marker.gene.1000.symbol.RData")
ARPC_seob <- subset(CRPC.confident.seob, tumor_subtype =='ARPC')

gene.id <- intersect(rownames(ARPC_seob), rownames(CCLE.log2.tpm.matrix.symbol))
gene.id <- intersect(gene.id, CCLE.rna.seq.marker.gene.1000.symbol)


ARPC.vcap.cor <- cor(ARPC_seob@assays$RNA$data[gene.id, ] %>% as.matrix() ,CCLE.log2.tpm.matrix.symbol[gene.id, 'VCAP_PROSTATE', drop=F], method = 'spearman')

ARPC.vcap.cor <- ARPC.vcap.cor %>% melt() %>% as.data.frame()
colnames(ARPC.vcap.cor) <- c('barcode', 'type','value')

# MSPC-PC3 ----------------------------------------------------------------


MSPC_seob     <- subset(CRPC.confident.seob, tumor_subtype =='MSPC')
gene.id       <- intersect(rownames(MSPC_seob), rownames(CCLE.log2.tpm.matrix.symbol))
gene.id       <- intersect(gene.id, CCLE.rna.seq.marker.gene.1000.symbol)

MSPC.pc3.cor  <- cor(MSPC_seob@assays$RNA$data[gene.id, ] %>% as.matrix() ,CCLE.log2.tpm.matrix.symbol[gene.id, 'PC3_PROSTATE', drop=F], method = 'spearman')


MSPC.pc3.cor            <- MSPC.pc3.cor %>% melt %>% as.data.frame()
colnames(MSPC.pc3.cor)  <- c('barcode', 'type','value')
MSPC.ARPC.cor           <- rbind(MSPC.pc3.cor, ARPC.vcap.cor)

MSPC.ARPC.cor$type <- factor(MSPC.ARPC.cor$type, levels = c('VCAP_PROSTATE','PC3_PROSTATE'))
MSPC.ARPC.cor %>% ggboxplot(x = 'type', y = 'value', color = 'type', add = 'point', size = 1, 
                            add.params = list(size = 2),
                            order = c('VCAP_PROSTATE', 'PC3_PROSTATE'))+
  scale_x_discrete(labels = c('VCAP_PROSTATE' = paste('VCAP_PROSTATE', sep = '\n', '(n = 686)'),
                              'PC3_PROSTATE' = paste('PC3_PROSTATE', sep = '\n', '(n = 2785)')))+
  stat_compare_means()


# Lineage analysis --------------------------------------------------------


load(file = '/home/liuxueying/met_pc_cell_line/pcdata/paper/data/section3/CRPC.confident.seob.RData')
MSPC.seob                   <- subset(CRPC.confident.seob, tumor_subtype == 'MSPC')
norm.seob                   <- readRDS("/home/liuxueying/met_pc_cell_line/met_pc_cell_line/pcdata/norm.seob.rds")
three.cell.type.seob        <- subset(norm.seob, cell_type_new =='BE'| cell_type_new=='LE'|cell_type_new=='NE')

stem.cell.pesudobulk <- AggregateExpression(
  three.cell.type.seob,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "cell_type_new",
  add.ident = NULL,
  slot = "counts",
  verbose = TRUE
)

stem.cell.pesudobulk                     <- stem.cell.pesudobulk$RNA
stem.cell.pesudobulk.cpm                 <- apply(stem.cell.pesudobulk, 2, function(x) { x/sum(x)*1000000 })
stem.cell.pesudobulk.log2cpm             <- log2(stem.cell.pesudobulk.cpm +1)

three.cell.type.seob  <- FindVariableFeatures(three.cell.type.seob, selection.method = "vst", nfeatures = 5000)
norm.5000.marker.gene <- head(VariableFeatures(three.cell.type.seob), 5000)

nor.data                                                        <- as.matrix(GetAssayData(MSPC.seob,layer='data'))
gg                 <- intersect(rownames(nor.data),norm.5000.marker.gene)
it.data            <- nor.data[gg,]  
norm.data          <- three.cell.type.aver.matrix[gg,]
cor.matrix <- foreach(it = iterators::iter(obj = it.data,by = 'column',chunksize= 100) ,.combine = 'rbind') %dopar% {
  rs         <- cor(it,norm.data,method = 'spearman')
  rs  
}
cor.matrix 
cor.max                             <- colnames(cor.matrix)[apply(cor.matrix, 1, function(x) which(x == max(x))[1])]
table(cor.max)
cor.anno.matrix                     <- cbind(cor.matrix, cor.max)
pheatmap(cor.matrix, show_rownames = F, cluster_rows = F, cluster_cols = F)



