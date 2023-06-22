## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
 if (!requireNamespace("curatedTCGAData", quietly = TRUE)) 
    BiocManager::install("curatedTCGAData")


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
library(curatedTCGAData)


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
curatedTCGAData(diseaseCode="BRCA", assays="*",dry.run=TRUE, version="2.0.1")


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
brca <- suppressMessages(
    curatedTCGAData(
        "BRCA",
        assays="RNASeqGene",
        dry.run=FALSE, 
        version="2.0.1")
)


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
brca


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
summary(assays(brca))


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
names(assays(brca))


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
xpr <- assays(brca)[[1]]
head(xpr[1:10,1:5])


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
nrow(xpr)


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
cnames <- colnames(xpr)
hpos <- gregexpr("-",cnames)
blah2 <- sapply(1:length(cnames),function(i) {substr(cnames[i],1,hpos[[i]][3]-1)})
idx <- !duplicated(blah2)

# subset values
blah2 <- blah2[idx]
xpr <- xpr[,idx]
colnames(xpr) <- blah2


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
pheno <- colData(brca)
colnames(pheno)[1:20]
head(pheno[,1:5])


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
ids_in_both <- intersect(pheno$patientID, colnames(xpr))
pheno <- pheno[which(pheno$patientID %in% ids_in_both),]
xpr <- xpr[,which(colnames(xpr) %in% ids_in_both)]

if (all.equal(pheno$patientID, colnames(xpr))!=TRUE) {
    midx <- match(pheno$patientID, colnames(xpr))
    xpr <- xpr[,midx]
}


## Always check that samples are in the same order in the tables you are comparing, and use `match()` to reorder them if necessary.

## If they are not in the same order, you are matching data to the wrong sample, and your results will be wrong.


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
pheno <- pheno[,c("patientID","PAM50.mRNA")]
table(pheno$PAM50.mRNA)


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
idx <- which(pheno$PAM50.mRNA %in% c("Luminal A","Basal-like"))
pheno <- pheno[idx,]
xpr <- xpr[,idx] 

dim(pheno)
dim(xpr)


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
suppressMessages(require(edgeR))


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
dge <- DGEList(
    counts = xpr,
    group = factor(
        pheno$PAM50.mRNA,
        levels=c("Luminal A","Basal-like")
    ) # should be a factor
)


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
head(dge$counts)


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
cpm(dge)[1:10,1:5]


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
dim(dge) #before 
tokeep <- rowSums(cpm(dge)>100) >= 2
dge <- dge[tokeep,]
dim(dge) #after


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
dge$samples$lib.size = colSums(dge$counts) # add up per sample
head(dge$samples)


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
dge <- calcNormFactors(dge)


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
plotMDS(
    dge, 
    col=as.numeric(dge$samples$group), 
    pch=16
)
legend(
    "bottomleft", 
    as.character(unique(dge$samples$group)),
    col=c(1,2), pch=16
    )


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
d1 <- estimateCommonDisp(dge, verbose=T)
d1 <- estimateTagwiseDisp(d1)


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
diffEx <- exactTest(d1, pair=c(1,2))


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
tt <- topTags(diffEx, n=10)


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
tt <- as.data.frame(
    topTags(diffEx, n=nrow(dge)
    )
)


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
rnd_pval <- runif(nrow(tt)) 
hist(
    rnd_pval,
    n=100,
    main="Histogram for uniformly-distributed data"
)
abline(v=0.05,col="red")


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
hist(tt$PValue,n=100, main="BRCA Luminal A vs Basal-Like: nominal pvalues")
abline(v=0.05,col="red")


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
qqplot(
    tt$PValue, 
    rnd_pval, 
    xlab="p-values from real data",
    ylab="X~U[0,1]",
    pch=16,cex=0.5
)

# x=y line as reference
abline(0,1,col="red")


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
diffEx2 <- decideTestsDGE(diffEx, 
    adjust.method="BH", 
    p.value=0.05
)
summary(diffEx2)


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
midx <- match(rownames(tt), rownames(diffEx2))
diffEx2 <- diffEx2[midx,]

cols <- rep("black",nrow(diffEx2))
cols[which(diffEx2>0 )]  <- "red"
cols[which(diffEx2<0)]  <- "blue"
# volcano plot
plot(tt$logFC,-log10(tt$PValue),pch=16,
    col=cols)
abline(v=0,lty=3)


## ---- class.source="codeblock",eval=TRUE--------------------------------------------------
write.table(tt,file="diffEx.results.txt",
    sep="\t",
    col=TRUE,
    row=TRUE,
    quote=FALSE
)

