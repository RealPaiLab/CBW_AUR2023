# Module 4 - exercises

```{r, class.source="codeblock",eval=TRUE}
if (!requireNamespace("yeastRNASeq", quietly = TRUE)) 
  BiocManager::install("yeastRNASeq")

library(yeastRNASeq)

data(geneLevelData)
gene_keep <- rowMeans(geneLevelData) > 2 

# What genes pass this threshold?
filtered <- geneLevelData[gene_keep,] 
str(filtered)
```

```{r, class.source="codeblock",eval=TRUE}
#  as.matrix(filtered): the count data in the right class
# phenoData: The sample information

group <- factor(rep(c("Mut", "WT"),each=2), levels = c("WT","Mut")) 
y <- DGEList(as.matrix(filtered), 
           group = group)  

## matrix of experimental design 
mod = model.matrix(~group, y)


## Normalize data
y <- calcNormFactors(y, method = "upperquartile")
y <- estimateDisp(y, mod)

fit = glmFit(y, mod)
lrt = glmLRT(fit, coef = 2)
```

```{r, class.source="codeblock",eval=TRUE}
DEGS = topTags(lrt, n=nrow(y))$table
```

```{r, class.source="codeblock",eval=TRUE}
## check out differentially expressed genes
head(DEGS)

DEGS_sig = DEGS[DEGS$FDR < 0.05,]
head(DEGS_sig)
dim(DEGS_sig)
```