library(RColorBrewer)
library(RColorBrewer)
## Explore missing data {-}

For this let's use a little script that converts a table into black and white squares to visualize missing data. For this, install the `plotrix` package

```{r, class.source="codeblock",eval=TRUE}
if (!requireNamespace("plotrix", quietly = TRUE)) install.packages("plotrix")

suppressMessages(require(plotrix))

#' show data missingness as a chequered matrix
#' 
#' @param x (matrix) data matrix.
#' @param outFile (char) path to file for printing graph
#' @param wd (numeric) width in inches
#' @param ht (numeric) height in inches
#' @return plots missingness matrix to file
#' @import plotrix
#' @export
plotMissMat <- function(x,xlab="columns",
		ylab="rows",border=NA) {
	
	x <- !is.na(x)
	class(x) <- "numeric"
	color2D.matplot(x,show.values=FALSE,axes=FALSE,
		cs1=c(1,0),cs2=c(1,0),cs3=c(1,0),border=border,
		cex=0.8,
		xlab=xlab,ylab=ylab)
}
```

Let's construct a 100 x 10 table and plot missingness.

```{r, class.source="codeblock",eval=TRUE}
dat <- matrix(rnorm(100*10,mean=0,sd=1),nrow=100)
colnames(dat) <- c(rep("case",5), rep("control",5))
dim(dat) 
```

Now let's add some unstructured missing data:

```{r, class.source="codeblock",eval=TRUE}
dat2 <- dat; 
dat2[sample(1000,10,F)] <- NA 
plotMissMat(dat2) 
```

And this is what structured missingness looks like. It could arise due to design limitations, or systematic issues with data collection:

```{r, class.source="codeblock",eval=TRUE}
dat3 <- dat; 
dat3[80:100,6:10] <- NA 
plotMissMat(dat3) # structured missing data

```