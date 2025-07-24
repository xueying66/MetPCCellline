library(data.table)
library(tidyverse)
library(ggpubr)
library(Seurat)
library(ggrepel)
library(foreach)
library(data.table)

load ('/home/liuxueying/met_pc_cell_line/processed.data/CCLE.transcriptome.RData')  
############################################################################################################################
################## Compute transcriptomic correlation between cell lines and MET500 samples ########################################
############################################################################################################################
load('/home/liuxueying/met_pc_cell_line/processed.data/BioklabData/liuke/KeData/MET500/MET500.RData')
flag                                     <- grepl(x=MET500.sample.meta$cancer.type,pattern='Prostate Adenocarcinoma') 
MET500.sample                            <- MET500.sample.meta$Run[flag]
tmp                                      <- MET500.sample.meta[MET500.sample,]

pc.MET.bone.polyA.sample  <- tmp$Run[tmp$biopsy.site == 'BONE'       & tmp$LibrarySelection == 'PolyA']
pc.MET.liver.polyA.sample <- tmp$Run[tmp$biopsy.site == 'LIVER'      & tmp$LibrarySelection == 'PolyA']
pc.MET.lymph.polyA.sample <- tmp$Run[tmp$biopsy.site == 'LYMPH_NODE' & tmp$LibrarySelection == 'PolyA']

pc.MET.bone.hs.sample  <- tmp$Run[tmp$biopsy.site == 'BONE'       & tmp$LibrarySelection == 'Hybrid Selection']
pc.MET.liver.hs.sample <- tmp$Run[tmp$biopsy.site == 'LIVER'      & tmp$LibrarySelection == 'Hybrid Selection']
pc.MET.lymph.hs.sample <- tmp$Run[tmp$biopsy.site == 'LYMPH_NODE' & tmp$LibrarySelection == 'Hybrid Selection']


pc.MET.sample.list <- list(pc.MET.bone.polyA.sample= pc.MET.bone.polyA.sample,
                           pc.MET.liver.polyA.sample = pc.MET.liver.polyA.sample,
                           pc.MET.lymph.polyA.sample = pc.MET.lymph.polyA.sample,
                           pc.MET.bone.hs.sample=pc.MET.bone.hs.sample,
                           pc.MET.liver.hs.sample = pc.MET.liver.hs.sample,
                           pc.MET.lymph.hs.sample = pc.MET.lymph.hs.sample
)



#  Compute PCa cell line correlation with all MET500 samples
source('/home/liuxueying/met_pc_cell_line/code/BioKLab.util.R')
MET500.TC.rs   <- pick.out.cell.line(expr.of.samples = MET500.log2.tpm.matrix[,pc.MET.sample.list %>% unlist],expr.of.cell.lines = CCLE.log2.tpm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
MET500.TC.rs$cell.line.median.cor %>% as.data.frame() %>% View()

################## Do metastasis site and platform matter? Actually NO ##########################  
MET500.site.TC.rs.list <- foreach (item = pc.MET.sample.list) %do% {
  rs   <- pick.out.cell.line(expr.of.samples = MET500.log2.tpm.matrix[,item],expr.of.cell.lines = CCLE.log2.tpm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
  rs
}
names(MET500.site.TC.rs.list)      <-  names(pc.MET.sample.list)
all.cell.line                      <-  colnames(CCLE.log2.tpm.matrix)
pc.cell.line                       <-  colnames(CCLE.log2.tpm.matrix)[colnames(CCLE.log2.tpm.matrix) %>% grepl(pattern='PROSTATE')]

MET500.site.cell.line.correlation.matrix <- foreach(item = MET500.site.TC.rs.list,.combine='cbind')  %do% {
  item$cell.line.median.cor[all.cell.line]
}
colnames(MET500.site.cell.line.correlation.matrix) <- names(MET500.site.TC.rs.list)


save(list = c('MET500.site.cell.line.correlation.matrix','MET500.site.TC.rs.list','MET500.TC.rs','MET500.site.cell.line.pc.correlation.matrix'),file = '/home/liuxueying/met_pc_cell_line/pcdata/paper/data/section2/MET500.RData')





# scRNA adenocacinoma validation-----------------------------------------------------------

source('/home/liuxueying/met_pc_cell_line/pcdata/method.R')
source('/home/liuxueying/met_pc_cell_line/code/Manuscript.Figure/ggplot.style.R')
#pc.cell.line.name
flag                         <- grepl(colnames(CCLE.log2.tpm.matrix), pattern ='PROSTATE')
pc.cell.line.log2.tpm.matrix <- CCLE.log2.tpm.matrix[,flag]
flag                         <- colnames(pc.cell.line.log2.tpm.matrix) != 'PRECLH_PROSTATE' ## remove no cancer 
pc.cell.line.log2.tpm.matrix <- pc.cell.line.log2.tpm.matrix[,flag]
pc.cell.line.name            <- colnames(pc.cell.line.log2.tpm.matrix) 

# TC rs function ----------------------------------------------------------

source('/home/liuxueying/met_pc_cell_line/code/BioKLab.util.R')
source('/home/liuxueying/met_pc_cell_line/pcdata/my_function.R')
load("/home/liuxueying/met_pc_cell_line/pcdata/CCLE.rna.seq.marker.gene.1000.symbol.RData")
load('/home/liuxueying/processed.data/CCLE.transcriptome.RData')
CCLE.log2.tpm.matrix.symbol <- change_rownames_ensemble_to_symbol(CCLE.log2.tpm.matrix)
TC.rs.list <- function(seob){
  TC.rs             <- pick.out.cell.line(expr.of.samples = seob@assays$RNA$data %>% as.matrix(), expr.of.cell.lines = CCLE.log2.tpm.matrix.symbol,marker.gene = CCLE.rna.seq.marker.gene.1000.symbol)
  median.cor        <- TC.rs$cell.line.median.cor %>% as.data.frame()
  pc.cell.line.name <- grep(pattern = '_PROSTATE', x = rownames(median.cor), value = T)
  pc.median.cor     <- median.cor[pc.cell.line.name,, drop=F]
  TC.rs.list        <-  list(median.cor = median.cor,
                             pc.median.cor = pc.median.cor,
                             correlation.matrix = TC.rs[["correlation.matrix"]])
}

load('/home/liuxueying/met_pc_cell_line/pcdata/SCP1244_14mCRPC/scp1224.seob.RData')

scp1224.malignant.seob <- subset(scp1224.seob, cell_type == 'prostate cancer cell')
adeno.seob             <- subset(scp1224.malignant.seob, biosample_id != '09171135') # remove small cell cacinoma
adeno.res              <- TC.rs.list(adeno.seob)
adeno.res$median.cor %>% View()

# plot --------------------------------------------------------------------

df                    <- adeno.res$median.cor %>% as.data.frame()
df$name               <- rownames(df)
df$prostate           <- ifelse(rownames(df) %in% pc.cell.line.name,'prostate','others')
colnames(df)[1]       <- 'correlation'
df_ordered            <- df[order(df$correlation, decreasing = T),]
df$rank               <- rank(df$correlation) 
df_ranked             <- df[order(df$rank),]


ggplot(df_ranked, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 3) 


#scRNA-seq site-specific correlation ---------------------------------------------------------- 
# bone
adeno.seob$organ__ontology_label %>% table()
bone.seob <- subset(adeno.seob, organ__ontology_label == 'bone tissue')
bone.res <- TC.rs.list(bone.seob)
bone.res$median.cor %>% View()
# lymph node
ln.seob   <- subset(adeno.seob, organ__ontology_label == 'lymph node')
ln.res    <- TC.rs.list(ln.seob)
ln.res$median.cor %>% View()



# bone ln scRNA spearman correlation
bone.median.cor                    <- bone.res$pc.median.cor %>% as.data.frame()
colnames(bone.median.cor)          <- 'bone.cor.median'
flag                               <- rownames(bone.median.cor) != 'PRECLH_PROSTATE'
bone.pc.median.cor                 <- bone.median.cor[flag,,drop = F ]

ln.median.cor                      <- ln.res$pc.median.cor %>% as.data.frame()
colnames(ln.median.cor)            <- 'ln.cor.median'
flag                <- rownames(ln.median.cor) !='PRECLH_PROSTATE'
ln.pc.median.cor    <- ln.median.cor[flag,,drop = F]
ln.pc.median.cor    <- ln.pc.median.cor[rownames(bone.pc.median.cor),,drop = F]
c.median.cor        <- cbind(ln.pc.median.cor, bone.pc.median.cor)
bone.ln.cor         <- cor(c.median.cor[,1],c.median.cor[,2], method = 'spearman')
df                  <- c.median.cor

ggplot(df, aes(x= ln.cor.median, y=bone.cor.median))+
  geom_point(size=8)+
  theme_bw()+
  xlab('Lymph node')+
  ylab('Bone')+
  ggtitle('')+
  xlim(c(0,0.3))+ylim(c(0,0.3))+geom_abline(slope = 1,intercept = 0)+
  annotate(
    'text',
    x = 0.2,            
    y = 0.15,           
    label = paste0("Spearman-rank correlation = ", round(bone.ln.cor, digits = 2)), 
    size = 5,          
    color = "black",     
    hjust = 0           
  )+theme_bw()
 













