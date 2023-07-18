#Making MA plot with deseq object including mutation in design, and using mutation as coef

ddseqset <- DESeqDataSetFromMatrix(countData = data.mtx, colData = metaData_mutations, design = ~ X12_25245350_T_C + Disease)

dds <- DESeq(ddseqset) 

#res <- results(dds, name = 'X12_25245350_T_C', filterFun = ihw, alpha = 0.01) #see if just this is better

#resLFC <- lfcShrink(dds, coef = 'X12_25245350_T_C', type = 'apeglm') #lfcshrink data

###

#without alpha = filtering 
res <- as.data.frame(results(dds, name = 'X12_25245350_T_C'))
res$Geneid <- data$Geneid
ggplot(res, aes(x=baseMean, y=log2FoldChange, colour=padj<0.01)) + scale_x_log10() + geom_point() + 
  scale_colour_manual(values=c("grey", "red")) 

###
#res$geneid <- data$Geneid

xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(res, xlim=xlim, ylim=ylim, main="apeglm") 

#Making MA plot with deseq object with just disease in design, and disease as coef 

ddseqset <- DESeqDataSetFromMatrix(countData = data.mtx, colData = metaData_mutations, design = ~ Disease)

dds <- DESeq(ddseqset) 

resultsNames(dds)

#res <- results(dds, name = 'Disease_celiac.disease.and.small.intestine.carcinoma_vs_celiac.disease', filterFun = ihw, alpha = 0.01)

#resLFC <- lfcShrink(dds, coef = 'Disease_celiac.disease.and.small.intestine.carcinoma_vs_celiac.disease', type = 'apeglm')

res <- as.data.frame(results(dds, name = 'Disease_celiac.disease.and.small.intestine.carcinoma_vs_celiac.disease'))
res$Geneid <- data$Geneid
ggplot(res, aes(x=baseMean, y=log2FoldChange, colour=padj<0.01)) + scale_x_log10() + geom_point() + 
  scale_colour_manual(values=c("grey", "red"))

res$geneid <- data$Geneid
#could save above as csv 

#The below doesnt work, find a way to use ggplot instead 
#xlim <- c(1,1e5); ylim <- c(-3,3)
#plotMA(res, xlim=xlim, ylim=ylim, main="apeglm") 

ggplot(res, aes(x=baseMean, y=log2FoldChange))
  geom_point(data=res[res$padj>=0.01, ], mapping=aes(colour="grey")) + 
               geom_point(data=res[res$padj<0.01, ], mapping=aes(colour="red")) +
                + scale_x_log10()            

