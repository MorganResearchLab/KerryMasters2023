


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

res <- results(dds, name = 'X12_25245350_T_C', filterFun = ihw, alpha = 0.01)
res$geneid <- data$Geneid
res
#run results with alpha = 0.01, and ihw first, then run lfcshrink

resLFC <- lfcShrink(dds, coef = 'X12_25245350_T_C', type = 'apeglm') 
resLFC$geneid <- data$Geneid
resLFC 

res.table <- merge(as.data.frame(resLFC[, c("log2FoldChange", "lfcSE", "geneid")]), 
                   as.data.frame(res[, c("stat", "pvalue", "padj", "weight", "geneid")]), by='geneid')
#added LFC and LFC standard error columns to the results table that includes the IHW, replaces LFC and LFCse columns from res table

head(res.table)

type(res.table$log2FoldChange)
type(res.table$lfcSE)
type(res.table$stat)
type(res.table$pvalue)
type(res.table$padj)
type(res.table$weight)

#try the ggplot method?

#res.table.noNA <- na.omit(res.table), this still doesnt help 'logical subscript contains NAs'
res.table.rnNA <- res.table[!is.na(res.table$padj),] #still the same issue 

ggplot(res.table.rnNA, aes(x=baseMean, y=log2FoldChange)) + scale_x_log10() +
  geom_point(data=res[res$padj>=0.01, ], mapping=aes(colour="grey")) +
geom_point(data=res[res$padj<0.01, ], mapping=aes(colour="red")) 

####

res.tableOrdered <- res[order(res$pvalue),] #can look at this to see ordered p values 
res.tableOrdered

summary(res) 
sum(res$padj < 0.05, na.rm = TRUE) #how many genes below padj > 0.05, 
sum(res$padj < 0.01, na.rm = TRUE) 
sum(res$padj < 0.001, na.rm = TRUE) 

res05 <- results(dds, alpha=0.05) #gives back results with false discovery rate adjusted for 0.05
summary(res05)

sum(res05$padj < 0.05, na.rm=TRUE) #redo this to see padj > 0.05,



#####
#applied above in results function as a whole

resIHW <- results(dds, filterFun = ihw)
summary(resIHW) #has adjp of 0.1 here
sum(resIHW$padj < 0.1, na.rm = TRUE) #3301
metadata(resIHW)$ihwResult #look into metadata to find the IHW summary

####


plotMA(resLFC, ylim=c(-2,2)) #A lot less data points that when looking only at disease design

xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(res.table, xlim=xlim, ylim=ylim, main="apeglm") #doesn't like it being in a table expects the data frame to have 3 columns, 
#two numeric ones for mean and log fold change, and a logical one for significance.


dPlotC <- plotCounts(dds, gene=which.min(res$padj), intgroup="Disease", 
                     returnData=TRUE)

ggplot(dPlotC, aes(x=Disease, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


#Still a bit unclear on order of DESeq design
#Less data points when looking setting mutation as coefficient?
#In which case should i be investigating the IHW results over the LFC results - look at deseq papers and see how info is presented, revisit papers
#Should i be setting false discovery rate with alpha = before moving on to graphs etc - use 0.01 for this
#


#comparing gene expr level with mutation count, gene exp = weighted variables +residual 
#different samples of same gene = vectors, the weighted variable applied to all of them = scailer, resodual vlues are per observation
#can plot residuals to ensure the model is viable, and look to see if theyre normaly distributed 
#can put samples into a model matrix, multiplies by vector of weights
#order matters as software sometimes takes last column in matrix, but can specify this

#volcano plots = reports log fold change on x axis and -log10 adj p value on y axis
#find top 5 upregulated genes, and top 5 downreg genes
#make volcano plot using ggplot 
#label volcano plot points with those genes, do this programatically
#ggrepel package 
#biomaRt package to map ensemble genes with gene symbols (e.g KRAS)
#Friday meeting = 11am 14th 
