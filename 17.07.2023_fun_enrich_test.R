library("DESeq2")
library("IHW")
library("ggplot2")
#library('EnhancedVolcano')
#library('biomaRt')
library('ggrepel')
library('org.Hs.eg.db')
library('enrichR')

setwd("~/Documents/Uni/Immunology/Project/Data/final_data_meta")

data <- read.csv('final_data_file')

metaData_mutations <- read.table('SBA_meta-mutations_merged.txt', sep = '\t', header = TRUE)


rownames(metaData_mutations) <- metaData_mutations$Sample

colnames(data) <- gsub('X.nfs.research.marioni.mdmorgan.KRAS_JulieStoeckus.bam.dir.', '', colnames(data))
colnames(data) <- gsub('_STARAligned.sortedByCoord.out.bam', '', colnames(data))

#make matrix so metadata can be aligned properly 
data.mtx <- as.matrix(data[,8:ncol(data)]) 

data.mtx <- data.mtx[,which(colnames(data.mtx) %in% metaData_mutations$Sample)] 

#Check row names of metadata and column names of count data are the same
identical(colnames(data.mtx), rownames(metaData_mutations)) #did it, TRUE

ddseqset <- DESeqDataSetFromMatrix(countData = data.mtx, colData = metaData_mutations, design = ~ X12_25245350_T_C + Disease) 
#includes mutation of interest

dds <- DESeq(ddseqset) 

resultsNames(dds)

res <- results(dds, name = 'X12_25245350_T_C', filterFun = ihw, alpha = 0.05)
res$geneid <- data$Geneid
res

resLFC <- lfcShrink(dds, coef = 'X12_25245350_T_C', type = 'apeglm') 
resLFC$geneid <- data$Geneid
resLFC 

res.table <- merge(as.data.frame(resLFC[, c("log2FoldChange", "lfcSE", "geneid")]), 
                   as.data.frame(res[, c("baseMean", "stat", "pvalue", "padj", "weight", "geneid")]), by='geneid')

head(res.table) #visualise 

#setting up and down reg


res.table$diffexpressed <- "NO"

res.table$diffexpressed[res.table$log2FoldChange > 1.2 & res.table$padj < 0.05] <- "UP"

res.table$diffexpressed[res.table$log2FoldChange < -1.2 & res.table$padj < 0.05] <- "DOWN"


#make vectors for colours 

mycolors <- c("blue", "red", "black")

#name them 

names(mycolors) <- c("DOWN", "UP", "NO") 

annots <- select(org.Hs.eg.db, keys=res.table$geneid, 
                 columns="SYMBOL", keytype="ENSEMBL")
head(annots) #this is what i was looking for, but its in numerical order. 

geneSymbRes <- merge(res.table, annots, by.x="geneid", by.y="ENSEMBL") #i think i've done it 
geneSymbRes

#now to label the genes on the volcano plot 

geneSymbRes$delabel <- NA
geneSymbRes$delabel[geneSymbRes$diffexpressed != "NO"] <- geneSymbRes$SYMBOL[geneSymbRes$diffexpressed != "NO"] 

#adapt below to make final plot
ggplot(data=geneSymbRes, aes(x=log2FoldChange, y= -log10(padj), col = diffexpressed, label = delabel )) + 
  geom_point() + theme_minimal() +
  geom_vline(xintercept=c(-1.2, 1.2), col="red3") +
  geom_hline(yintercept=-log10(0.05), col="red3") +
  scale_colour_manual(values = mycolors) +
  geom_text_repel()

#make code to head 5 most upreg and 5 most downreg genes 

head(geneSymbRes[order(geneSymbRes$log2FoldChange, decreasing = FALSE),]) #shows the 6 most downreg genes

head(geneSymbRes[order(geneSymbRes$log2FoldChange, decreasing = TRUE),]) #shows the 6 most upreg genes


#order above by log fold change instead so that its not just alphabetically ordered

#foxp3 
geneSymbRes[geneSymbRes$SYMBOL %in% c('FOXP3', "IL2RB", "TIGIT", "CD4", "CD8", "IL17A"),]

#look at enriched pathways 


###

listEnrichrSites()
setEnrichrSite("Enrichr") #human genes set

enr.dbs <- listEnrichrDbs() #listing the enrichment pathways 

head(enr.dbs)

#adding the functions i was found etc, check these, just adds them as a string so not sure what the vignette meant?

enr.dbs <- c("GO_Cellular_Component_2017", "GO_Molecular_Function_2017", "GO_Biological_Process_2017")

enriched.up <- enrichr(geneSymbRes$SYMBOL[geneSymbRes$diffexpressed == 'UP'], enr.dbs)
names(enriched.up)

#upreg EGFR signalling, we know this is directly upstream from KRAS 

enriched.down <- enrichr(geneSymbRes$SYMBOL[geneSymbRes$diffexpressed == 'DOWN'], enr.dbs)
enriched.down


enr.bio.pros.up <- enriched.up[["GO_Biological_Process_2017"]] 
enr.bio.pros.down <- enriched.down[["GO_Biological_Process_2017"]]


#head(enrichr(res.table$SYMBOL, c('ST8SIA3'))[[1]]) #very unsure how to do this

####


#can plot a bar graph with x = odds ratio and y = term, 

ggplot(data = enr.bio.pros.up, mapping =x = Odds.Ratio, y = Term) +
  geom_bar(stat = 'identity')

ggplot(data = enr.bio.pros.down, mapping = aes(x = Odds.Ratio, y = Term)) +
  geom_bar(stat = 'identity') + 
  guides(guide_axis(title = Term, angle = 40, n.dodge = 210))


ggplot(data = enr.bio.pros.down, mapping = aes(x = Odds.Ratio, y = Term)) +
  geom_bar(stat = 'identity') + 
  scale_y_discrete(discrete_scale(aesthetics = NULL, scale_name = Term, palette = "Set1", breaks = waiver()))

ggplot(data = enr.bio.pros.up, mapping = aes(x = Odds.Ratio, y = Term)) +
  geom_bar(stat = 'identity') + 
  scale_y_discrete(aesthetics = NULL, scale_name = Term, palette = "Set1", breaks = waiver())

#will have to fix properly for paper
