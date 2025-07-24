load('/home/liuxueying/met_pc_cell_line/processed.data/organize.mutation.data.RData')


############################################################################################################################
################## Compute genes that are differentialy mutated between metastatic and primary samples, and those that are highly mutated in metastatic samples
############################################################################################################################
TCGA.mutation.matrix                           <- PC.mutation.matrix.list$TCGA$SNP + PC.mutation.matrix.list$TCGA$INS + PC.mutation.matrix.list$TCGA$DEL
TCGA.mutation.matrix[TCGA.mutation.matrix > 0] <- 1
TCGA.mutation.freq                             <- apply(TCGA.mutation.matrix,1,function(x) sum(x > 0) / length(x))

SU2C.mutation.matrix                           <- PC.mutation.matrix.list$SU2C$SNP + PC.mutation.matrix.list$SU2C$INS + PC.mutation.matrix.list$SU2C$DEL
SU2C.mutation.matrix[SU2C.mutation.matrix > 0] <- 1

CCLE.mutation.matrix                           <- PC.mutation.matrix.list$CCLE$SNP + PC.mutation.matrix.list$CCLE$INS + PC.mutation.matrix.list$CCLE$DEL
CCLE.mutation.matrix[CCLE.mutation.matrix > 0] <- 1


liver.sample.idx            <- which(grepl(x=colnames(SU2C.mutation.matrix), pattern = '@Liver'))
bone.sample.idx             <- which(grepl(x=colnames(SU2C.mutation.matrix), pattern = '@Bone'))
ln.sample.idx               <- which(grepl(x=colnames(SU2C.mutation.matrix), pattern = '@LN'))
SU2C.mutation.freq.liver    <- apply(SU2C.mutation.matrix[,liver.sample.idx],1, function(x) sum(x > 0) / length(x))
SU2C.mutation.freq.bone     <- apply(SU2C.mutation.matrix[,bone.sample.idx], 1, function(x) sum(x > 0) / length(x))
SU2C.mutation.freq.ln       <- apply(SU2C.mutation.matrix[,ln.sample.idx],   1, function(x) sum(x > 0) / length(x))
tmp                         <- as.integer(SU2C.mutation.freq.liver > 0.05) + as.integer(SU2C.mutation.freq.bone > 0.05) + as.integer(SU2C.mutation.freq.ln > 0.05) ##why choose 0.05
MET.highly.mutated.gene     <- names(SU2C.mutation.freq.liver)[tmp >=2]  #Let us consider three tissues   select MET highly mutated genes


dm.analysis  <- function(PRI.mutation.freq, MET.mutation.freq, PRI.mutation.matrix,MET.mutation.matrix) {
  g.vec <- c( names(PRI.mutation.freq)[PRI.mutation.freq > 0.01 ],names(MET.mutation.freq)[MET.mutation.freq > 0.01 ]) %>% unique
  rs.df <- foreach(g = g.vec,.combine = 'rbind') %do% {
    x       <- matrix(0,nrow = 2,ncol = 2) 
    x[1,1]  <- sum(PRI.mutation.matrix[g,] == 0)
    x[1,2]  <- sum(PRI.mutation.matrix[g,] >0)
    x[2,1]  <- sum(MET.mutation.matrix[g,] == 0)
    x[2,2]  <- sum(MET.mutation.matrix[g,] >0)
    p.value <- fisher.test(x)$p.value  
    
    #p.value <- 1- pbinom(q = sum(MET.mutation.matrix[g,] >0),size = ncol(MET.mutation.matrix),prob = PRI.mutation.freq[g] )
    delta   <- MET.mutation.freq[g]  - PRI.mutation.freq[g] 
    data.frame(delta = delta, p.value = p.value)
  }
  rownames(rs.df) <- g.vec
  rs.df$fdr       <- p.adjust(rs.df$p.value,method = 'fdr')
  rs.df$gene      <- rownames(rs.df)
  rs.df
}

bone.dm.rs   <- dm.analysis(PRI.mutation.freq = TCGA.mutation.freq, MET.mutation.freq = SU2C.mutation.freq.bone,  PRI.mutation.matrix = TCGA.mutation.matrix, MET.mutation.matrix = SU2C.mutation.matrix[,bone.sample.idx])
liver.dm.rs  <- dm.analysis(PRI.mutation.freq = TCGA.mutation.freq, MET.mutation.freq = SU2C.mutation.freq.liver, PRI.mutation.matrix = TCGA.mutation.matrix, MET.mutation.matrix = SU2C.mutation.matrix[,liver.sample.idx])
ln.dm.rs     <- dm.analysis(PRI.mutation.freq = TCGA.mutation.freq, MET.mutation.freq = SU2C.mutation.freq.ln,    PRI.mutation.matrix = TCGA.mutation.matrix, MET.mutation.matrix = SU2C.mutation.matrix[,ln.sample.idx])
dm.gene.df   <- rbind(bone.dm.rs,liver.dm.rs)
dm.gene.df   <- rbind(dm.gene.df,ln.dm.rs)
tmp          <- filter(dm.gene.df,fdr < 0.05)  
MET.dm.gene  <- table(tmp$gene) %>% as.data.frame() %>% filter(Freq >=2) %>% dplyr::select(Var1)  %>% unlist %>% as.character()##filter(Freq >=2) 在两种转移类型以上都显著差异的基因
full.gene                                      <- c(MET.highly.mutated.gene,MET.dm.gene) %>% unique ## select 73 metastasis-related genes


########## Compute mutation burden ########################
tmp                         <- ddply(PC.maf.data.list$CCLE,.(dcast.id), function(x) sum(x$Variant_Type == 'SNP' & x$Variant_Classification != 'silent'))
CCLE.mutation.burden        <- tmp$V1
names(CCLE.mutation.burden) <- tmp$dcast.id    ####为什么只选择SNP 不选择 del ins

tmp                         <- ddply(PC.maf.data.list$SU2C,.(dcast.id), function(x) sum(x$Variant_Type == 'SNP' & x$Variant_Classification != 'silent'))
SU2C.mutation.burden        <- tmp$V1
names(SU2C.mutation.burden) <- tmp$dcast.id

tmp                         <- ddply(PC.maf.data.list$TCGA,.(dcast.id), function(x) sum(x$Variant_Type == 'SNP' & x$Variant_Classification != 'silent'))
TCGA.mutation.burden        <- tmp$V1
names(TCGA.mutation.burden) <- tmp$dcast.id

resemble.gene.cnt                              <- apply(CCLE.mutation.matrix[full.gene,],2,sum)
plot(x = log2(CCLE.mutation.burden), y= resemble.gene.cnt[names(CCLE.mutation.burden)]) #! X and Y are highly correlated ####
sort(log2(CCLE.mutation.burden)) %>% plot # LNCAPCLONEFGC,DU145,22RV1,MDAPCA2B are hyper-mutated  

( resemble.gene.cnt[names(CCLE.mutation.burden)] / CCLE.mutation.burden ) %>% sort # good, Vcap rank is high    

# We use 9 as the cut-off to pick out hyper-mutated samples, the fisher's exact test is to prove the cutoff makes sense
boxplot(log2(TCGA.mutation.burden),log2(SU2C.mutation.burden),log2(CCLE.mutation.burden))
log2(TCGA.mutation.burden) %>% sort %>% plot
log2(SU2C.mutation.burden) %>% sort %>% plot
log2(CCLE.mutation.burden) %>% sort %>% plot

MMR.gene        <- c('MLH1','MLH3','MSH2','MSH6','MSH3','PMS1','PMS2','POLE','POLD1') # https://www.nature.com/articles/ncomms15180/
TCGA.MSI.sample <- names(TCGA.mutation.burden)[log2(TCGA.mutation.burden) > 9]
SU2C.MSI.sample <- names(SU2C.mutation.burden)[log2(SU2C.mutation.burden) > 9]
CCLE.MSI.sample <- names(CCLE.mutation.burden)[log2(CCLE.mutation.burden) > 9]
MSI.sample      <- c(TCGA.MSI.sample,SU2C.MSI.sample,CCLE.MSI.sample)    

MMR.mutation.matrix <- cbind(TCGA.mutation.matrix[MMR.gene,],SU2C.mutation.matrix[MMR.gene,],CCLE.mutation.matrix[MMR.gene,])
MMR.mut.sample      <- colnames(MMR.mutation.matrix)[apply(MMR.mutation.matrix,2,sum) >0 ]  

MSI.no.sample       <- setdiff(colnames(MMR.mutation.matrix),MSI.sample)  
MMR.no.mut.sample   <- setdiff(colnames(MMR.mutation.matrix),MMR.mut.sample)   

x      <- matrix(0,nrow = 2, ncol = 2) #row: MMR mutation status; col:MSI mutation status
x[1,1] <- intersect(MMR.mut.sample,MSI.sample)       %>% length
x[1,2] <- intersect(MMR.mut.sample,MSI.no.sample)    %>% length
x[2,1] <- intersect(MMR.no.mut.sample,MSI.sample)    %>% length
x[2,2] <- intersect(MMR.no.mut.sample,MSI.no.sample) %>% length
fisher.test(x)  


#Next, let us check hotspot mutations
non.silent.PC.SU2C.maf.data  <- PC.maf.data.list$SU2C %>% filter(Variant_Classification != 'Silent')  #### 转移癌样本体细胞突变数据  
PC.SU2C.hotspot.mutation     <- table (non.silent.PC.SU2C.maf.data$mutation.code) %>% as.data.frame %>% filter(Freq >= 3)  %>%   dplyr::select(Var1) %>% unlist %>% as.character()  
df                           <- filter(PC.maf.data.list$CCLE, mutation.code %in% PC.SU2C.hotspot.mutation) %>% as.data.frame  
View(df) #"1@6257785@RPL22" resembled by three hyper-mutated cell lines; this mutation is usually seen in hyper-mutated samples, another piece of proof! 


# DE analysis between MSI sample and non-MSI samples, a piece of proof we should use hyper-mutated cell line to study the metastasis mechanism of hyper-mutated prostate cancer 

load('/home/liuxueying/met_pc_cell_line/TCGA/Prostate Adenocarcinoma.RData')  
TCGA.MSI.patient    <- names(TCGA.mutation.burden)[log2(TCGA.mutation.burden) >= 9]
TCGA.MSI.sample     <- PC.maf.data.list$TCGA$Tumor_Sample_Barcode[match(x= TCGA.MSI.patient,table = PC.maf.data.list$TCGA$dcast.id)]
TCGA.MSI.sample     <- intersect(TCGA.MSI.sample, colnames(log2.fpkm.matrix)) 
TCGA.non.MSI.sample <- setdiff(colnames(log2.fpkm.matrix),TCGA.MSI.sample )    

source('./code/BioKLab.util.R')
de.res <- perform.DE.analysis.between.TRE.and.CON(CON.log2.read.count.matrix = log2.read.count.matrix[,TCGA.non.MSI.sample],
                                                  TRE.log2.read.count.matrix = log2.read.count.matrix[,TCGA.MSI.sample],
                                                  CON.log2.tpm.matrix        = log2.tpm.matrix[,TCGA.non.MSI.sample],
                                                  TRE.log2.tpm.matrix        = log2.tpm.matrix[,TCGA.MSI.sample]
)   ####得到超突变样本和非超突变样本差异表达基因

up.gene <- rownames(de.res)[de.res$log2FoldChange > 1  & de.res$padj < 0.05] # GO analysis shows cell cycle - related process elevated, and T cell abundance elevated 
dn.gene <- rownames(de.res)[de.res$log2FoldChange < -1 & de.res$padj < 0.05] # Gamma-delta T cell marker TRGV9 down-regulated













