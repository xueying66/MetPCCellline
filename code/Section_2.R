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



