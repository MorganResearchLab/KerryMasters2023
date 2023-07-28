library("DESeq2")
library("IHW")
library("ggplot2")
#library('EnhancedVolcano')
#library('biomaRt')
library('ggrepel')
library('org.Hs.eg.db')
library('enrichR')
library('stringr')

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
names(enriched.down)


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

#use stringr to sort strings str_wrap

#will have to fix properly for paper

#ask about KRAS mutations, not sure how to pick extra for writing about, dont cover all the ones in the new doc
#ask about the code that sorted the data into the 'final_data' and if i need it
#protein database uniprot to look at kras structure, maybe use as picture for intro

###
#Attempt with stringr

#OG ugly graphs
ggplot(data = enr.bio.pros.up, mapping = aes(x = Odds.Ratio, y = Term)) +
  geom_bar(stat = 'identity')



ggplot(data = enr.bio.pros.down, mapping = aes(x = Odds.Ratio, y = Term)) +
  geom_bar(stat = 'identity')


#Adding y axis text to str_wrap

#enr.up_yaxis <- str_wrap(colnames(enr.bio.pros.up), width = 70, indent = 0, exdent = 0, whitespace_only = TRUE)

#enr.down_yaxis <- str_wrap(colnames(enr.bio.pros.down), width = 70, indent = 0, exdent = 0, whitespace_only = TRUE)

ggplot(data = enr.bio.pros.up, mapping = aes(x = Odds.Ratio, y = enr.up_yaxis)) +
  geom_bar(stat = 'identity')

#^this doesnt work, aesthetics must be length 1 or same length as data 

ggplot(data = enr.bio.pros.up, mapping = aes(x = Odds.Ratio, y = Term)) +
  geom_bar(stat = 'identity') +
  scale_y_discrete(labels = function(y) str_wrap(y, width = 65), )
#^ this is a lot better but there are overlaps

ggplot(data = enr.bio.pros.up, mapping = aes(x = Odds.Ratio, y = Term)) +
  geom_bar(stat = 'identity') +
  scale_y_discrete(labels = function(y) str_wrap(y, width = 65), ) +
  guides(guide_axis(check.overlap = TRUE))
#^no difference, very ugly overlaps 

ggplot(data = enr.bio.pros.up, mapping = aes(x = Odds.Ratio, y = reorder(Term, Odds.Ratio))) +
  geom_bar(stat = 'identity') +
  scale_y_discrete(labels = function(y) str_wrap(y, width = 65))
#removed guides for now, ordered the y axis by odds ratio 

ggplot(data = enr.bio.pros.down, mapping = aes(x = Odds.Ratio, y = reorder(Term, Odds.Ratio))) +
  geom_bar(stat = 'identity') +
  scale_y_discrete(labels = function(y) str_wrap(y, width = 65))
#same as above but for down reg


###


#starting over but making new column in data

enr.bio.pros.up$up_y.axis <- str_wrap(colnames(enr.bio.pros.up), width = 70, indent = 0, exdent = 0, whitespace_only = TRUE)

enr.bio.pros.down$down_y.axis <- str_wrap(colnames(enr.bio.pros.down), width = 70, indent = 0, exdent = 0, whitespace_only = TRUE)

#Error in `$<-.data.frame`(`*tmp*`, up_y.axis, value = c("Term", "Overlap",  : replacement has 9 rows, data has 210

#$ for dataframes, [] for matrix

#/n is seperator for new line

#start writing on pen and paper 


#statistically significant 

 #%in% causes logical values, this selects columns, only want to select rows where value is below a specific point

#indexing rows, not column so use right hand of [,], data[data$subset<5,] (always < then = in that order)
#can use %in% if there are certain values you want to pick out. has to be specific though, not threshhold

#sig.enr.up <- enr.bio.pros.up[enr.bio.pros.up$Adjusted.P.value<0.1,]

#start from scratch if unsure/wrong
#keep seperate file for incorrect code

###

#new column for the str_wrap strings
enr.bio.pros.up$up_y.axis <- str_wrap(enr.bio.pros.up$Term, width = 70, indent = 0, exdent = 0, whitespace_only = TRUE)
enr.bio.pros.up$up_y.axis

enr.bio.pros.down$down_y.axis <- str_wrap(enr.bio.pros.down$Term, width = 70, indent = 0, exdent = 0, whitespace_only = TRUE)
enr.bio.pros.down$down_y.axis


#remake the ggplots

enr.bio.pros.up.sig <- enr.bio.pros.up[enr.bio.pros.up$'Adjusted.P.value'<0.1,]
enr.bio.pros.up.sig

enr.bio.pros.down.sig <- enr.bio.pros.down[enr.bio.pros.down$'Adjusted.P.value'<0.1,]
enr.bio.pros.up.sig


ggplot(data = enr.bio.pros.up.sig, mapping = aes(x = Odds.Ratio, y = reorder(up_y.axis, Odds.Ratio))) +
  geom_bar(stat = 'identity') +
  ylab('Function') + xlab('Odds Ratio')

  
ggplot(data = enr.bio.pros.down.sig, mapping = aes(x = Odds.Ratio, y = reorder(down_y.axis, Odds.Ratio))) +
  geom_bar(stat = 'identity') +
  ylab('Function') + xlab('Odds Ratio')

