require(plyr)
require(dplyr)
require(reshape2)
require(foreach)
require(data.table)



##PROCESS CCLE transcriptom data

if(file.exists('./processed.data/CCLE.transcriptome.RData') == FALSE){
  tmp     <- fread(input = '/home/liuxueying/met_pc_cell_line/raw.data/CCLE/CCLE_RNAseq_rsem_genes_tpm_20180929.txt',header=TRUE)
  gene.id <- tmp$gene_id
  data    <- tmp[,3:ncol(tmp)] %>% as.matrix
  remove.ensemble.version.id <- function(x){
    v <- strsplit(x = x,split = "\\.") %>% unlist
    v[1]
  }
  clean.gene.id        <- sapply(gene.id,remove.ensemble.version.id)
  names(clean.gene.id) <- NULL
  CCLE.log2.tpm.matrix <- log2(data + 1)
  rownames(CCLE.log2.tpm.matrix) <- clean.gene.id
  CCLE.median                    <- apply(CCLE.log2.tpm.matrix,1,median)
  CCLE.expressed.gene            <- names(CCLE.median)[CCLE.median > 1]
  tmp                            <- CCLE.log2.tpm.matrix[CCLE.expressed.gene,]
  tmp.rank                       <- apply(tmp,2,rank)
  rank.mean                      <- apply(tmp.rank,1,mean)
  rank.sd                        <- apply(tmp.rank,1,sd)
  plot(x=rank.mean,y=rank.sd)
  lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
  CCLE.rna.seq.marker.gene.1000   <- names(sort(rank.sd,decreasing =TRUE))[1:1000]
  save(file = '/home/liuxueying/met_pc_cell_line/processed.data/CCLE.transcriptome.RData',list = c('CCLE.log2.tpm.matrix','CCLE.rna.seq.marker.gene.1000'))
}


#### PROCESS organize mutation data
ENST2HGNC                 <- read.csv("/home/liuxueying/met_pc_cell_line/raw.data/ensemble/ENST2HGNC.csv") %>% unique
protein.coding.symbol.vec <- filter(ENST2HGNC,Gene.type == 'protein_coding') %>% dplyr::select(HGNC.symbol) %>% unlist %>% unique %>% as.character()
protein.coding.symbol.vec <- protein.coding.symbol.vec[protein.coding.symbol.vec != '']

##########################  Organize CCLE  somatic mutation data ##########################  
PC.sample_info <- read.csv("/home/liuxueying/met_pc_cell_line/raw.data/CCLE/sample_info.csv") %>%
  filter(primary_disease == 'Prostate Cancer')
tmp            <- fread(input = "./raw.data/CCLE/CCLE_mutations.csv",stringsAsFactors = FALSE) %>%
  filter(Variant_Classification %in% c("5'Flank","3'Flank") == FALSE & Hugo_Symbol != 'Unknown' & DepMap_ID %in% PC.sample_info$DepMap_ID) %>%
  unique        %>%
  as.data.frame 

PC.CCLE.maf.data  <- merge(x= tmp, y = PC.sample_info[,c('DepMap_ID','CCLE_Name')], by = 'DepMap_ID') #Map DepMap id to CCLE name

remove.ensemble.version.id <- function(x){
  v <- strsplit(x = x,split = "\\.") %>% unlist
  v[1]
}

clean.tr.id                  <- sapply(PC.CCLE.maf.data$Annotation_Transcript,remove.ensemble.version.id)
names(clean.tr.id)           <- NULL
PC.CCLE.maf.data$clean.tr.id <- clean.tr.id
tmp                          <- merge(x= PC.CCLE.maf.data, y= ENST2HGNC, by.x = 'clean.tr.id', by.y = 'Transcript.stable.ID') # Let me re-assign gene symbols based on ENSEMBLE annotation
PC.CCLE.maf.data             <- tmp
PC.CCLE.maf.data$dcast.id    <- PC.CCLE.maf.data$CCLE_Name


##########################  Organize SU2C somatic mutation data (metastatic prostate cancer, there are patients overlapping with MET500) ##################
SU2C.meta             <- fread(input = "/home/liuxueying/met_pc_cell_line/raw.data/cBioPortal/prad_su2c_2019/data_clinical_sample.txt",  stringsAsFactors=FALSE,skip = 4,header=TRUE) %>%
  unique %>%
  filter(PATHOLOGY_CLASSIFICATION == 'Adenocarcinoma' & TISSUE_SITE %in% c('Prostate','Unknown')  == FALSE) %>%   # Here, we only want metastatic samples and Adenocarcinoma samples
  as.data.frame
tmp                   <- fread(file = "./raw.data/cBioPortal/prad_su2c_2019//data_mutations_extended.txt", stringsAsFactors=FALSE) %>% 
  unique   %>% 
  filter(Variant_Classification %in% c("5'Flank","3'Flank") == FALSE & Hugo_Symbol != 'Unknown') %>%
  as.data.frame
PC.SU2C.maf.data      <- merge(x= tmp, y = SU2C.meta,  by.x='Tumor_Sample_Barcode',by.y = 'SAMPLE_ID') # I only want metastatic samples

tmp                          <- merge(x= PC.SU2C.maf.data, y= ENST2HGNC, by.x = 'Transcript_ID', by.y = 'Transcript.stable.ID')  # Let me re-assign gene symbols based on ENSEMBLE annotation
PC.SU2C.maf.data             <- tmp
PC.SU2C.maf.data$dcast.id    <- paste(PC.SU2C.maf.data$PATIENT_ID,PC.SU2C.maf.data$TISSUE_SITE,sep ='@'  )


########################## Organize TCGA  somatic mutation data ##########################
PC.TCGA.maf.data          <- read.delim("/home/liuxueying/met_pc_cell_line/raw.data/cBioPortal/prad_tcga_pan_can_atlas_2018/data_mutations_extended.txt", stringsAsFactors=FALSE) %>% 
  unique   %>% 
  filter(Variant_Classification %in% c("5'Flank","3'Flank") == FALSE & Hugo_Symbol != 'Unknown' ) 

tmp                       <- merge(x= PC.TCGA.maf.data, y= ENST2HGNC, by.x = 'Transcript_ID', by.y = 'Transcript.stable.ID') # Let me re-assign gene symbols based on ENSEMBLE annotation
PC.TCGA.maf.data          <- tmp


get.TCGA.patient.id <- function(x) {
  paste(strsplit(x, split ='-') %>% unlist  %>% head(3),collapse = '-')
}
PC.TCGA.maf.data$dcast.id   <- sapply(PC.TCGA.maf.data$Tumor_Sample_Barcode,get.TCGA.patient.id)



get.mutation.matrix.list <- function(maf.data){
  tmp             <- filter(maf.data,HGNC.symbol %in% protein.coding.symbol.vec) # Here, I only consider protein-coding genes
  dcast.id.vec    <- maf.data$dcast.id %>% as.character() %>% unique
  
  mutation.matrix.list <- foreach(variant_type = c('SNP','INS','DEL')) %do% {
    VAR.matrix.df              <- reshape2::dcast(tmp[tmp$Variant_Type == variant_type,],HGNC.symbol ~ dcast.id,fun.aggregate = function(x) ifelse(length(x) >= 1, 1, 0)   )
    rownames(VAR.matrix.df)    <- VAR.matrix.df$HGNC.symbol
    VAR.matrix.df$HGNC.symbol  <- NULL
    VAR.matrix                 <- as.matrix(VAR.matrix.df)
    
    delta.dcast.id         <- setdiff(dcast.id.vec,colnames(VAR.matrix)) 
    if( length(delta.dcast.id) > 0 ){
      col.delta.matrix           <- matrix(0,nrow = nrow(VAR.matrix), ncol = length(delta.dcast.id))
      rownames(col.delta.matrix) <- rownames(VAR.matrix)
      colnames(col.delta.matrix) <- delta.dcast.id
      VAR.matrix                 <- cbind(VAR.matrix, col.delta.matrix)
    }
    
    delta.symbol           <- setdiff(protein.coding.symbol.vec,rownames(VAR.matrix))
    if( length(delta.symbol) > 0 ){
      row.delta.matrix           <- matrix(0,nrow = length(delta.symbol), ncol = ncol(VAR.matrix))
      rownames(row.delta.matrix) <- delta.symbol
      colnames(row.delta.matrix) <- colnames(VAR.matrix)
      VAR.matrix                 <- rbind(VAR.matrix,row.delta.matrix)
    }
    
    VAR.matrix[protein.coding.symbol.vec,dcast.id.vec]
  }
  names(mutation.matrix.list) <- c('SNP','INS','DEL')
  mutation.matrix.list
}


PC.mutation.matrix.list <- list(CCLE = get.mutation.matrix.list(PC.CCLE.maf.data),
                                SU2C = get.mutation.matrix.list(PC.SU2C.maf.data),
                                TCGA = get.mutation.matrix.list(PC.TCGA.maf.data)
)

PC.CCLE.maf.data$mutation.code <- paste(PC.CCLE.maf.data$Chromosome,PC.CCLE.maf.data$Start_position,PC.CCLE.maf.data$HGNC.symbol,sep ='@')
PC.SU2C.maf.data$mutation.code <- paste(PC.SU2C.maf.data$Chromosome,PC.SU2C.maf.data$Start_Position,PC.SU2C.maf.data$HGNC.symbol,sep ='@')
PC.TCGA.maf.data$mutation.code <- paste(PC.TCGA.maf.data$Chromosome,PC.TCGA.maf.data$Start_Position,PC.TCGA.maf.data$HGNC.symbol,sep ='@')
PC.maf.data.list <- list(CCLE = PC.CCLE.maf.data,
                         SU2C = PC.SU2C.maf.data,
                         TCGA = PC.TCGA.maf.data
)

save(file ='/home/liuxueying/met_pc_cell_line/processed.data/organize.mutation.data.RData', 
     list= c('PC.mutation.matrix.list','PC.maf.data.list')
)





