library("DESeq2")
library("IHW")
library("ggplot2")

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

res <- results(dds, name = 'X12_25245350_T_C', filterFun = ihw, alpha = 0.1)
res$geneid <- data$Geneid
res

resLFC <- lfcShrink(dds, coef = 'X12_25245350_T_C', type = 'apeglm') 
resLFC$geneid <- data$Geneid
resLFC 

res.table <- merge(as.data.frame(resLFC[, c("log2FoldChange", "lfcSE", "geneid")]), 
                   as.data.frame(res[, c("baseMean", "stat", "pvalue", "padj", "weight", "geneid")]), by='geneid')

head(res.table) #visualise 

#trying to get the MA plot to work in ggplot

#res.table.noNA <- na.omit(res.table), this still doesnt help 'logical subscript contains NAs'
res.table.rmNA <- res.table[!is.na(res.table$padj),] #still the same issue 

ggplot(data = res.table.rmNA, mapping = aes(x=baseMean, y=log2FoldChange)) +
  geom_point(data=res.table.rmNA[res.table.rmNA$padj>=0.05, ], mapping=aes(color= 'grey')) +
  geom_point(data=res.table.rmNA[res.table.rmNA$padj<0.05, ], mapping=aes(color="red")) + scale_x_log10()

#isnt mapping the colours correctly, its putting the grey and red as a legend

ggplot(data = res.table.rmNA) +
  geom_point(data=res.table.rmNA[res.table.rmNA$padj>=0.05, ], mapping=aes(x = baseMean, y=log2FoldChange, color= 'grey')) +
  geom_point(data=res.table.rmNA[res.table.rmNA$padj<0.05, ], mapping=aes(x = baseMean, y=log2FoldChange, color="red")) + scale_x_log10()

#exact same

ggplot(data = res.table.rmNA, mapping = aes(x=baseMean, y=log2FoldChange)) +
  geom_point(data=res.table.rmNA[res.table.rmNA$padj>=0.05, ], mapping=aes(scale_colour_identity('grey'))) +
  geom_point(data=res.table.rmNA[res.table.rmNA$padj<0.05, ], mapping=aes(scale_colour_identity("red"))) + scale_x_log10()

#not valid aesthetics, 


ggplot(data = res.table.rmNA, mapping = aes(x=baseMean, y=log2FoldChange)) +
  geom_point(data=res.table.rmNA[res.table.rmNA$padj>=0.05, ], color= 'grey') +
  geom_point(data=res.table.rmNA[res.table.rmNA$padj<0.05, ], color="red") + scale_x_log10()

sum(is.na(res.table$padj)) #9176
