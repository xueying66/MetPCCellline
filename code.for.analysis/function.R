#####Function to pick out cell line##########
pick.out.cell.line <- function(expr.of.samples,expr.of.cell.lines,marker.gene){
  marker.gene           <- intersect(rownames(expr.of.samples),(marker.gene))  
  marker.gene           <- intersect(rownames(expr.of.cell.lines),(marker.gene)) 
  correlation.matrix    <- cor(expr.of.samples[marker.gene,],expr.of.cell.lines[marker.gene,],method='spearman')
  cell.line.median.cor  <- apply(correlation.matrix,2,median) %>% sort(decreasing = TRUE)
  best.cell.line        <- names(cell.line.median.cor)[1]
  p.value.vec           <- foreach(cell.line= setdiff(names(cell.line.median.cor),best.cell.line),.combine='c') %do% {
    v                 <- correlation.matrix[,cell.line]
    p.value           <- wilcox.test(correlation.matrix[,best.cell.line],v,alternative = 'greater',paired = TRUE)$p.value
  }
  names(p.value.vec) <- setdiff(names(cell.line.median.cor),best.cell.line)
  fdr.vec            <- p.adjust(p.value.vec,method='fdr')
  list(cell.line.median.cor=cell.line.median.cor,best.cell.line=best.cell.line,compare.fdr.vec=fdr.vec,correlation.matrix = correlation.matrix )
}



##### ensemble to symbol #####
hgnc.data                 <- read.delim("./external.data/HGNC/hgnc_complete_set.txt", stringsAsFactors=FALSE)

symbol.to.ensemble <- function(symbol) {
  mapping <- hgnc.data[match(x=symbol,table = hgnc.data$symbol),]
  mapping <- mapping$ensembl_gene_id
  return(mapping)
}
symbol.to.ensemble('S100A6')

ensembl.to.symbol <- function(ensembl){
  mapping <- hgnc.data[match(x=ensembl, table = hgnc.data$ensembl_gene_id),]
  mapping <- mapping$symbol
  return(mapping)
}
ensembl.to.symbol('ENSG00000089685')





# exchange rownames symbol to ensemble ------------------------------------

change_rownames_symbol_to_ensemble <- function(expression.matrix){
  match.pos          <- match(x= rownames(expression.matrix), table = hgnc.data$symbol)
  ense.id            <- hgnc.data[match.pos,]$ensembl_gene_id
  df                 <- data.frame(ense.id = ense.id, symbol.id = rownames(expression.matrix))
  df                 <- df[complete.cases(df),]
  expression.matrix  <- expression.matrix[df$symbol.id, ,drop=F]
  if(identical(rownames(expression.matrix), df$symbol.id)){
    rownames(expression.matrix) <- df$ense.id
  }else{
    expression.matrix
  }
  return(expression.matrix)
}


# exchange rownames ensemble to symbol ------------------------------------


change_rownames_ensemble_to_symbol <- function(expression.matrix){
  match.pos <- match(x=rownames(expression.matrix), table = hgnc.data$ensembl_gene_id)
  symbol.id <- hgnc.data[match.pos,]$symbol
  df <- data.frame(ense.id = rownames(expression.matrix), symbol.id=symbol.id)
  df <- df[complete.cases(df),]
  expression.matrix <- expression.matrix[df$ense.id,, drop=F]
  if(identical(rownames(expression.matrix), df$ense.id)){
    rownames(expression.matrix) <- df$symbol.id
  }else{
    expression.matrix
  }
  return(expression.matrix)
}


# TC rs function ----------------------------------------------------------

source('/home/liuxueying/met_pc_cell_line/code/BioKLab.util.R')
source('/home/liuxueying/met_pc_cell_line/pcdata/my_function.R')
load("~/met_pc_cell_line/pcdata/CCLE.rna.seq.marker.gene.1000.symbol.RData")
load ('./processed.data/CCLE.transcriptome.RData')
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





