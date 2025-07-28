###########################  Organize SU2C CNV data ###############################
data_cna_hg19         <- read.delim("/home/liuxueying/met_pc_cell_line/raw.data/cBioPortal/prad_su2c_2019/data_cna_hg19.seg")
SU2C.meta             <- fread(input = "/home/liuxueying/met_pc_cell_line/raw.data/cBioPortal/prad_su2c_2019/data_clinical_sample.txt",  stringsAsFactors=FALSE,skip = 4,header=TRUE) %>%
  unique %>%
  filter(PATHOLOGY_CLASSIFICATION == 'Adenocarcinoma' & TISSUE_SITE %in% c('Prostate','Unknown')  == FALSE) %>%   # I only want metastatic samples
  as.data.frame

metastasis.sample.cnv.data  <- filter(data_cna_hg19, ID %in% unique(SU2C.meta$SAMPLE_ID)) %>% dplyr::select(c('ID','chrom','loc.start','loc.end','seg.mean'))
load('/home/liuxueying/met_pc_cell_line/output/hg19.gene.info.R.output/hg19.gene.info.RData')
SU2C.obj                    <- CNSeg(metastasis.sample.cnv.data)
SU2C.gene.cnv               <- getRS(object = SU2C.obj,by = 'gene',imput = FALSE,XY=TRUE,geneMap = hg19.gene.info,what = "mean")@rs
SU2C.gene.cnv               <- SU2C.gene.cnv[,5:ncol(SU2C.gene.cnv)]
rownames(SU2C.gene.cnv)     <- SU2C.gene.cnv$genename
SU2C.gene.cnv$genename      <- NULL
SU2C.gene.cnv.matrix        <- as.matrix(SU2C.gene.cnv)

save(SU2C.gene.cnv.matrix, file = 'home/liuxueying/met_pc_cell_line/pcdata/paper/data/section1/SU2C.gene.cnv.matrix.RData')
########################### Organize TCGA CNV data ##################################
data_cna_hg19               <- read.delim("./raw.data/cBioPortal/prad_tcga_pan_can_atlas_2018/data_cna_hg19.seg")
data_cna_hg19               <- dplyr::select(data_cna_hg19,c('ID','chrom','loc.start','loc.end','seg.mean'))
TCGA.obj                    <- CNSeg(data_cna_hg19)
TCGA.gene.cnv               <- getRS(object = TCGA.obj,by = 'gene',imput = FALSE,XY=TRUE,geneMap = hg19.gene.info,what = "mean")@rs
TCGA.gene.cnv               <- TCGA.gene.cnv[,5:ncol(TCGA.gene.cnv)]
rownames(TCGA.gene.cnv)     <- TCGA.gene.cnv$genename
TCGA.gene.cnv$genename      <- NULL
TCGA.gene.cnv.matrix        <- as.matrix(TCGA.gene.cnv)

save(TCGA.gene.cnv.matrix, file = '/home/liuxueying/met_pc_cell_line/pcdata/paper/data/section1/TCGA.gene.cnv.matrix.RData')
########################### Organize CCLE CNV data   WES##################################

data_cna_hg19               <- read.csv("./pcdata/paper/data/section1/CCLE_wes_segment_cn.csv", header = T)
data_cna_hg19               <- dplyr::select(data_cna_hg19,c('DepMap_ID','Chromosome','Start','End','Segment_Mean'))
colnames(data_cna_hg19)     <- c('ID','chrom','loc.start','loc.end','seg.mean')
 
sample_info          <- read.csv('/home/liuxueying/met_pc_cell_line/pcdata/paper/data/section1/sample_info.csv', header = T)
flag                 <- grepl(pattern = '_PROSTATE', x = sample_info$CCLE_Name)
pc.ccle.info         <- sample_info[flag,]
pc.ccle.depmap.id    <- pc.ccle.info$DepMap_ID

PC.data_cna_hg19               <- data_cna_hg19[data_cna_hg19$ID %in% pc.ccle.depmap.id,]
load('/home/liuxueying/met_pc_cell_line/output/hg19.gene.info.R.output/hg19.gene.info.RData')
CCLE.PC.obj                    <- CNSeg(PC.data_cna_hg19)
CCLE.PC.gene.cnv               <- getRS(object = CCLE.PC.obj,by = 'gene',imput = FALSE,XY=TRUE,geneMap = hg19.gene.info,what = "mean")@rs
CCLE.PC.gene.cnv               <- CCLE.PC.gene.cnv[,5:ncol(CCLE.PC.gene.cnv)]
rownames(CCLE.PC.gene.cnv)     <- CCLE.PC.gene.cnv$genename
CCLE.PC.gene.cnv$genename      <- NULL
CCLE.PC.gene.cnv.matrix        <- as.matrix(CCLE.PC.gene.cnv)


tmp                                 <- pc.ccle.info[pc.ccle.info$DepMap_ID %in% colnames(CCLE.PC.gene.cnv.matrix),]
CCLE.PC.gene.cnv.matrix             <- CCLE.PC.gene.cnv.matrix[,tmp$DepMap_ID]
identical(tmp$DepMap_ID, colnames(CCLE.PC.gene.cnv.matrix))
colnames(CCLE.PC.gene.cnv.matrix)   <- tmp$CCLE_Name

save('CCLE.PC.gene.cnv.matrix', file = '/home/liuxueying/met_pc_cell_line/pcdata/paper/data/section1/CCLE.PC.gene.cnv.matrix.RData')


# metastasis and primary samples differential cnv ------------------------------------------------------

load(file = '/home/liuxueying/met_pc_cell_line/pcdata/paper/data/section1/TCGA.gene.cnv.matrix.RData')
load(file = '/home/liuxueying/met_pc_cell_line/pcdata/paper/data/section1/SU2C.gene.cnv.matrix.RData')


c.gene <- intersect(rownames(TCGA.gene.cnv.matrix), rownames(SU2C.gene.cnv.matrix))
rs.df  <- foreach(g = c.gene,.combine = 'rbind') %do% {
  p.value <- wilcox.test(TCGA.gene.cnv.matrix[g,],SU2C.gene.cnv.matrix[g,])$p.value
  delta   <- SU2C.gene.cnv.matrix[g, bone.id] %>% median - TCGA.gene.cnv.matrix[g,] %>% median
  data.frame(delta = delta, p.value = p.value,gene = g)
}  ####计算TCGA SU2C 差异cnv 基因 (使用wilcox 方法)
rownames(rs.df)       <- rs.df$gene
rs.df$p.adj           <- p.adjust(rs.df$p.value,method = 'bonferroni')


rs.df$sig <- ifelse(rs.df$p.adj<0.05 & rs.df$delta > 0.3, 'up',
                    ifelse(rs.df$p.adj <0.05 & rs.df$delta < -0.3, 'dn', 'no'))


# select differential CNV gene site specific ----------------------------------------
load(file = './pcdata/paper/data/section1/TCGA.gene.cnv.matrix.RData')
load(file = './pcdata/paper/data/section1/SU2C.gene.cnv.matrix.RData')


c.gene <- intersect(rownames(TCGA.gene.cnv.matrix), rownames(SU2C.gene.cnv.matrix))
rs.df.bone  <- foreach(g = c.gene,.combine = 'rbind') %do% {
  p.value <- wilcox.test(TCGA.gene.cnv.matrix[g,],SU2C.gene.cnv.matrix[g,bone.id])$p.value
  delta   <- SU2C.gene.cnv.matrix[g, bone.id] %>% median - TCGA.gene.cnv.matrix[g,] %>% median
  data.frame(delta = delta, p.value = p.value,gene = g)
}  
rownames(rs.df.bone)       <- rs.df.bone$gene
rs.df.bone$p.adj           <- p.adjust(rs.df.bone$p.value,method = 'bonferroni')



c.gene <- intersect(rownames(TCGA.gene.cnv.matrix), rownames(SU2C.gene.cnv.matrix))
rs.df.liver  <- foreach(g = c.gene,.combine = 'rbind') %do% {
  p.value <- wilcox.test(TCGA.gene.cnv.matrix[g,],SU2C.gene.cnv.matrix[g,liver.id])$p.value
  delta   <- SU2C.gene.cnv.matrix[g, liver.id] %>% median - TCGA.gene.cnv.matrix[g,] %>% median
  data.frame(delta = delta,p.value = p.value,gene = g)
}  
rownames(rs.df.liver)       <- rs.df.liver$gene
rs.df.liver$p.adj           <- p.adjust(rs.df.liver$p.value,method = 'bonferroni')


c.gene <- intersect(rownames(TCGA.gene.cnv.matrix), rownames(SU2C.gene.cnv.matrix))
rs.df.LN  <- foreach(g = c.gene,.combine = 'rbind') %do% {
  p.value <- wilcox.test(TCGA.gene.cnv.matrix[g,],SU2C.gene.cnv.matrix[g,LN.id])$p.value
  delta   <- SU2C.gene.cnv.matrix[g, LN.id] %>% median - TCGA.gene.cnv.matrix[g,] %>% median
  data.frame(delta = delta,p.value = p.value,gene = g)
} 
rownames(rs.df.LN)       <- rs.df.LN$gene
rs.df.LN$p.adj           <- p.adjust(rs.df.LN$p.value,method = 'bonferroni')


