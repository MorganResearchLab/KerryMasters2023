library("DESeq2")
library("IHW")
library("ggplot2")
library('ggrepel')
library('org.Hs.eg.db')
library('enrichR')
library('stringr')

#Setting up all the RNA sequencing data from all the files so that it's all merged into one file. Merged by chromosome then by Geneid column.
setwd("~/Documents/Uni/Immunology/Project/Data/quant.dir")

chromes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
             "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM")

chrome_list <- list()

wd <- getwd()

for(x in seq_along(chromes)){   
  chrome_list[[chromes[x]]] <- list.files(wd, pattern=paste0("*", chromes[x], "-counts.txt"), full.names=TRUE)
}

data_list <- list()      

for(x in seq_along(chrome_list)){
  file_list <- chrome_list[[chromes[x]]]
  merged_data <- NULL
  for(file in file_list){
    data <- read.table(file, header = TRUE) 
    if (is.null(merged_data)) {
      merged_data <- data[,!colnames(data) %in% c("Chr", 'Start', 'End', 'Strand', 'Length')]
    } else {
      merged_data <- merge(merged_data, data[,!colnames(data) %in% c("Chr", 'Start', 'End', 'Strand', 'Length')], by = 'Geneid', all = TRUE)
    }
  }
  data_list[[x]] <- merged_data
}


final <- NULL

for (x in seq_along(data_list)) {
  if (is.null(final)) {
    final <- data_list[[x]]
  } else {
    final <- rbind.data.frame(final, data_list[[x]]) 
  }
}


write.csv(final, 'final_merged_data_file')

#Reading metadata and making DESeq2 object.

data <- read.csv('final_merged_data_file')

metaData_mutations <- read.table('../final_data_meta/SBA_meta-mutations_merged.txt', sep = '\t', header = TRUE)

rownames(metaData_mutations) <- metaData_mutations$Sample

colnames(data) <- gsub('X.nfs.research.marioni.mdmorgan.KRAS_JulieStoeckus.bam.dir.', '', colnames(data))
colnames(data) <- gsub('_STARAligned.sortedByCoord.out.bam', '', colnames(data))

#Make matrix so metadata can be aligned properly 
data.mtx <- as.matrix(data[,3:ncol(data)]) 

data.mtx <- data.mtx[,which(colnames(data.mtx) %in% metaData_mutations$Sample)] 

#Check row names of metadata and column names of count data are the same
identical(colnames(data.mtx), rownames(metaData_mutations)) #TRUE

ddseqset <- DESeqDataSetFromMatrix(countData = data.mtx, colData = metaData_mutations, design = ~ X12_25245350_T_C + Disease)

#Running DESeq and collecting results, IHW and LFC then merged into table

dds <- DESeq(ddseqset) 

resultsNames(dds)

res <- results(dds, name = 'X12_25245350_T_C', filterFun = ihw, alpha = 0.01)
res$geneid <- data$Geneid
res
#run results with alpha = 0.01, and IHW first, then run lfcshrink

resLFC <- lfcShrink(dds, coef = 'X12_25245350_T_C', type = 'apeglm') 
resLFC$geneid <- data$Geneid
resLFC 

res.table <- merge(as.data.frame(resLFC[, c("log2FoldChange", "lfcSE", "geneid")]), 
                   as.data.frame(res[, c("baseMean", "stat", "pvalue", "padj", "weight", "geneid")]), by='geneid')
#added LFC and LFC standard error columns to the results table that includes the IHW, replaces LFC and LFCse columns from res table

head(res.table)

#MA plot of the base mean vs the log2 fold change

res.table.rmNA <- res.table[!is.na(res.table$padj),]

ggplot(data = res.table.rmNA, mapping = aes(x=baseMean, y=log2FoldChange)) +
  geom_point(data=res.table.rmNA[res.table.rmNA$padj>=0.05, ], color= 'grey') +
  geom_point(data=res.table.rmNA[res.table.rmNA$padj<0.05, ], color="red") + scale_x_log10()


#Setting colours for Volcano plot

mycolors <- c("blue", "red", "black")

names(mycolors) <- c("DOWN", "UP", "NO")

#Setting column for determining differentially expressed genes

res.table$diffexpressed <- "NO"

res.table$diffexpressed[res.table$log2FoldChange > 1.2 & res.table$padj < 0.05] <- "UP"

res.table$diffexpressed[res.table$log2FoldChange < -1.2 & res.table$padj < 0.05] <- "DOWN"

#Using org.Hs.eg.db package to label ensemble genes 

annots <- select(org.Hs.eg.db, keys=res.table$geneid, 
                 columns="SYMBOL", keytype="ENSEMBL")
head(annots) 

geneSymbRes <- merge(res.table, annots, by.x="geneid", by.y="ENSEMBL")  
geneSymbRes

geneSymbRes$delabel <- NA
geneSymbRes$delabel[geneSymbRes$diffexpressed != "NO"] <- geneSymbRes$SYMBOL[geneSymbRes$diffexpressed != "NO"]

#Volcano plot

ggplot(data=geneSymbRes, aes(x=log2FoldChange, y= -log10(padj), col = diffexpressed, label = delabel )) + 
  geom_point() + theme_minimal() +
  geom_vline(xintercept=c(-1.2, 1.2), col="red3") +
  geom_hline(yintercept=-log10(0.05), col="red3") +
  scale_colour_manual(values = mycolors) +
  geom_text_repel()

#Differential expressed genes headed, checks for influential immune system genes

head(geneSymbRes[order(geneSymbRes$log2FoldChange, decreasing = FALSE),]) #shows the 6 most downreg genes

head(geneSymbRes[order(geneSymbRes$log2FoldChange, decreasing = TRUE),]) #shows the 6 most upreg genes

geneSymbRes[geneSymbRes$SYMBOL %in% c('FOXP3', "IL2RB", "TIGIT", "CD4", "CD8", "IL17A"),] #seeing if important immune genes are significant



#Functional enrichment testing

listEnrichrSites()
setEnrichrSite("Enrichr") #human genes set

enr.dbs <- listEnrichrDbs() #listing the enrichment pathways 

head(enr.dbs)

enr.dbs <- c("GO_Cellular_Component_2017", "GO_Molecular_Function_2017", "GO_Biological_Process_2017")

enriched.up <- enrichr(geneSymbRes$SYMBOL[geneSymbRes$diffexpressed == 'UP'], enr.dbs)
names(enriched.up)

enriched.down <- enrichr(geneSymbRes$SYMBOL[geneSymbRes$diffexpressed == 'DOWN'], enr.dbs)
names(enriched.down)

enr.bio.pros.up <- enriched.up[["GO_Biological_Process_2017"]] 
enr.bio.pros.down <- enriched.down[["GO_Biological_Process_2017"]]

#Bar plots for functional enrichment

enr.bio.pros.up$up_y.axis <- str_wrap(enr.bio.pros.up$Term, width = 70, indent = 0, exdent = 0, whitespace_only = TRUE)
enr.bio.pros.up$up_y.axis

enr.bio.pros.down$down_y.axis <- str_wrap(enr.bio.pros.down$Term, width = 70, indent = 0, exdent = 0, whitespace_only = TRUE)
enr.bio.pros.down$down_y.axis

enr.bio.pros.up.sig <- enr.bio.pros.up[enr.bio.pros.up$'Adjusted.P.value'<0.1,]
enr.bio.pros.up.sig

enr.bio.pros.down.sig <- enr.bio.pros.down[enr.bio.pros.down$'Adjusted.P.value'<0.1,]
enr.bio.pros.up.sig


ggplot(data = enr.bio.pros.up.sig, mapping = aes(x = Odds.Ratio, y = reorder(up_y.axis, Odds.Ratio))) +
  geom_bar(stat = 'identity') +
  ylab('Biological Function Term') + xlab('Odds Ratio') +
  scale_y_discrete(limits=rev)

down_func <- c('regulation of TORC1 signaling (GO:1903432)', 'microglial cell activation (GO:0001774)', 'single strand break repair (GO:0000012)',
               'retinal ganglion cell axon guidance (GO:0031290)', 'response to amino acid (GO:0043200)', 'positive regulation of helicase activity (GO:0051096)')

ggplot(data = enr.bio.pros.down.sig, mapping = aes(x = Odds.Ratio, y = reorder(down_y.axis, Odds.Ratio))) +
  geom_bar(stat = 'identity') +
  ylab('Biological Function Term') + xlab('Odds Ratio') +
  ylim('regulation of TORC1 signaling (GO:1903432)', 'microglial cell activation (GO:0001774)', 'single strand break repair (GO:0000012)',
       'retinal ganglion cell axon guidance (GO:0031290)', 'response to amino acid (GO:0043200)', 'positive regulation of helicase activity (GO:0051096)')

