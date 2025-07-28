library(DESeq2)
library(GSVA)

################### Yes, pc1.tail.gene are cMET co-expressed genes! ######################
cor.rs    <- cor(CCLE.log2.tpm.matrix[MET,] %>% c,CCLE.log2.tpm.matrix[pc1.tail.gene,] %>% t , method = 'spearman')
cor.rd.rs <- cor(CCLE.log2.tpm.matrix[MET,] %>% c,CCLE.log2.tpm.matrix[sample(size = 200, x = 1:nrow(CCLE.log2.tpm.matrix)),] %>% t , method = 'spearman')
cbind(cor.rd.rs %>% c, cor.rs %>% c) %>% boxplot() # Compared with random 200 genes, correlation significantly higher!!!


with.PC3.corr         <- cor(CCLE.log2.tpm.matrix[CCLE.rna.seq.marker.gene.1000,'PC3_PROSTATE'], CCLE.log2.tpm.matrix[CCLE.rna.seq.marker.gene.1000,],method = 'spearman')
plot(x = CCLE.log2.tpm.matrix[MET,], y = with.PC3.corr %>% c) # wow, correlation with PC3 is correlated with  cMET expression
cor.test(x = CCLE.log2.tpm.matrix[MET,], y = with.PC3.corr %>% c,method = "spearman")

cor(x = CCLE.log2.tpm.matrix[MET,], y = with.PC3.corr %>% c)




################### Yes, pc1.tail.gene are cMET co-expressed genes! ######################
cor.rs    <- cor(CCLE.log2.tpm.matrix[MET,] %>% c,CCLE.log2.tpm.matrix[pc1.tail.gene,] %>% t , method = 'spearman')
cor.rd.rs <- cor(CCLE.log2.tpm.matrix[MET,] %>% c,CCLE.log2.tpm.matrix[sample(size = 200, x = 1:nrow(CCLE.log2.tpm.matrix)),] %>% t , method = 'spearman')
cbind(cor.rd.rs %>% c, cor.rs %>% c) %>% boxplot() # Compared with random 200 genes, correlation significantly higher!!!


with.PC3.corr         <- cor(CCLE.log2.tpm.matrix[CCLE.rna.seq.marker.gene.1000,'PC3_PROSTATE'], CCLE.log2.tpm.matrix[CCLE.rna.seq.marker.gene.1000,],method = 'spearman')
plot(x = CCLE.log2.tpm.matrix[MET,], y = with.PC3.corr %>% c) # wow, correlation with PC3 is correlated with  cMET expression
cor.test(x = CCLE.log2.tpm.matrix[MET,], y = with.PC3.corr %>% c,method = "spearman")

cor(x = CCLE.log2.tpm.matrix[MET,], y = with.PC3.corr %>% c)

# VCAP overexpression vs control DE analysis -----------------------------------------


load("~/met_pc_cell_line/pcdata/vcap_cMET_high_bulk_data/gene.expression.RData")

cts                                            <- 2^log2.read.count.matrix-1
cts <- cts[,-6]
cts <- cts[,-6]
coldata                                        <- data.frame(condition = factor(c(rep('wild_type',4),rep('cMET_high' ,2))),row.names = colnames(cts))

class(coldata)
rownames(coldata)                              <- colnames(cts)
any(colnames(cts)==rownames(coldata))
cts                                            <- round(cts)
dds                                            <- DESeqDataSetFromMatrix(countData = cts,
                                                                         colData = coldata,
                                                                         design = ~ condition)


dds$condition                                  <- relevel(dds$condition, ref = "wild_type")
dds                                            <- DESeq(dds)
res                                            <- results(dds)
res.df                                         <- as.data.frame(res)
res.df                                         <- res.df[complete.cases(res.df),]
up.gene                                     <- filter(res.df, log2FoldChange > 1 & padj<0.05)
dn.gene                                     <- filter(res.df, log2FoldChange < -1 & padj<0.05)

load("/home/liuxueying/met_pc_cell_line/pcdata/vcap_cMET_high_bulk_data/gene.expression.RData")
colnames(log2.tpm.matrix) <- c(rep('wild_type',4), rep('cMET_high',4))
vcap_bulk                 <- log2.tpm.matrix

gseascore                                     <- gsva(vcap_bulk,
                                                    sig.list.en, 
                                                    method="ssgsea",
                                                    verbose=FALSE,
                                                    ssgsea.norm=T) 





# engineered correlation with MSPC ------------------------------------------------------
##LNCAP knockout
load(file = './pcdata/GSE175975_LNCAP_scrna/lncap.six.sample.ave.cpm.RData')
lncap.six.sample.ave.cpm 

lncap.ko  <- lncap.six.sample.ave.cpm[, 5:6]
lncap.ko.id <- colnames(lncap.ko)




### LNCAP enza 
load("~/met_pc_cell_line/pcdata/GSE162225/RData/gene.expression.RData")
head(log2.tpm.matrix)
lncap.enza <- log2.tpm.matrix
colnames(lncap.enza) <- c(rep('enza_d0', 3), rep('enza_d14',3))


lncap.enza <- lncap.enza[,c(4,5,6)]
colnames(lncap.enza) <- c('LNCaP.enza.1','LNCaP.enza.2','LNCaP.enza.3')
lncap.enza.id <- colnames(lncap.enza)


####VCAP MET OE 
load("~/met_pc_cell_line/pcdata/vcap_cMET_high_bulk_data/gene.expression.RData")
cMET_high_expr <- log2.tpm.matrix
cMET_high_expr <- cMET_high_expr[,-6]
cMET_high_expr <- cMET_high_expr[,-6]
cMET_high_expr <- cMET_high_expr[, 5:6]
colnames(cMET_high_expr) <- c('MET_OE1', 'MET_OE2')

vcap.oe    <- cMET_high_expr
vcap.oe.id <- colnames(vcap.oe)


gene.id <- intersect(rownames(CCLE.log2.tpm.matrix), rownames(lncap.ko))
gene.id <- intersect(rownames(lncap.enza), gene.id)
gene.id <- intersect(rownames(vcap.oe), gene.id)
lncap.ko      <- lncap.ko[gene.id,]
lncap.enza    <- lncap.enza[gene.id,]
vcap.oe       <- vcap.oe[gene.id,]
ccle.tpm      <- CCLE.log2.tpm.matrix[gene.id,]

all.sample.tpm <- cbind(lncap.ko, lncap.enza)
all.sample.tpm <- cbind(all.sample.tpm, vcap.oe)
all.sample.tpm <- cbind(all.sample.tpm, ccle.tpm)









# GSE181374 organoid stem---------------------------------------------------------------

load("/home/liuxueying/met_pc_cell_line/met_pc_cell_line/pcdata/stem.organoid/RData/gene.expression.RData")

stem.bulk.organoid.tpm  <- log2.tpm.matrix

# correlation with MSPC cell 
stem.bulk.organoid.tpm[1:4, 1:4]
CCLE.log2.tpm.matrix[1:4, 1:4]
gene.id <- intersect(rownames(stem.bulk.organoid.tpm), rownames(CCLE.log2.tpm.matrix))
stem.tpm <- stem.bulk.organoid.tpm[gene.id,]
ccle.tpm <- CCLE.log2.tpm.matrix[gene.id,]
identical(rownames(stem.tpm), rownames(ccle.tpm))
stem.organoid.ccle <- cbind(stem.tpm, ccle.tpm)
dim(stem.organoid.ccle)

load(file = './pcdata/paper/data/section3/CRPC.confident.seob.RData')
CRPC.tumor.seob <- subset(CRPC.confident.seob, tumor_subtype != 'normal')

stem.cell.pesudobulk <- AggregateExpression(
  CRPC.tumor.seob,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "tumor_subtype",
  add.ident = NULL,
  slot = "counts",
  verbose = TRUE
)

stem.cell.pesudobulk                     <- stem.cell.pesudobulk$RNA
stem.cell.pesudobulk.cpm                 <- apply(stem.cell.pesudobulk, 2, function(x) { x/sum(x)*1000000 })
stem.cell.pesudobulk.log2cpm             <- log2(stem.cell.pesudobulk.cpm +1)

stem.cell.pesudobulk.log2cpm             <- change_rownames_symbol_to_ensemble(stem.cell.pesudobulk.log2cpm)

marker.gene           <- intersect(rownames(stem.cell.pesudobulk.log2cpm),CCLE.rna.seq.marker.gene.1000)  
marker.gene           <- intersect(rownames(stem.organoid.ccle),marker.gene) 
ccle.stem.cor         <- cor(stem.cell.pesudobulk.log2cpm[marker.gene,],stem.organoid.ccle[marker.gene,],method='spearman')
ccle.stem.cor         <- ccle.stem.cor %>% t %>% as.data.frame()
stem.id      <- colnames(stem.tpm)
ccle.stem.cor[pc.cell.line.name,]
ccle.stem.cor[stem.id,]
df           <- ccle.stem.cor[, 'MSPC', drop = FALSE]
df$name      <- rownames(df)

df$prostate  <- ifelse(rownames(df) %in% pc.cell.line.name, 
                       "prostate", 
                       ifelse(rownames(df) %in% stem.id, 'stem_organoid', 'others'))

colnames(df)[1] <- "correlation"
df_ranked <- df %>% arrange(correlation)
df_ranked$rank <- 1:nrow(df_ranked)
ggplot(df_ranked, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 2.5) +
  geom_point(
    data = subset(df_ranked, name %in% pc.cell.line.name),
    color = "firebrick",
    size = 4
  ) +
  geom_point(
    data = subset(df_ranked, name %in% stem.id),
    color = "steelblue",
    size = 4
  )



# MSPC organoid lineage annotation -------------------------------------------------

load( file = './pcdata/paper/data/section5/three.cell.type.seob.RData')
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

three.cell.type.seob <- FindVariableFeatures(three.cell.type.seob, selection.method = "vst", nfeatures = 5000)

norm.5000.marker.gene <- head(VariableFeatures(three.cell.type.seob), 5000)


nor.data                                                        <- stem.bulk.organoid.tpm.symbol
gg                 <- intersect(rownames(nor.data),norm.5000.marker.gene)
it.data            <- nor.data[gg,]  ##scrna expression matrix
norm.data          <- stem.cell.pesudobulk.log2cpm[gg,]
cor.matrix <- foreach(it = iterators::iter(obj = it.data,by = 'column',chunksize= 100) ,.combine = 'rbind') %dopar% {
  rs         <- cor(it,norm.data,method = 'spearman')
  rs  
}
cor.matrix  ## result get  one cell - one cell line  correlation
cor.max                             <- colnames(cor.matrix)[apply(cor.matrix, 1, function(x) which(x == max(x))[1])]
cor.anno.matrix                     <- cbind(cor.matrix, cor.max)

pheatmap(cor.matrix, show_rownames = F, cluster_rows = F, cluster_cols = F)
