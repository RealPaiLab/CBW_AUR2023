
# Module 3: Fitting generalized linear models
rm(list=ls())

## Merging data

## Missing data

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

require(plotrix)
dat <- matrix(rnorm(100*10,mean=0,sd=1),nrow=100)
colnames(dat) <- c(rep("case",5), rep("control",5))
dim(dat)
plotMissMat(dat) # no missing data

dat2 <- dat; 
dat2[sample(1000,10,F)] <- NA 
plotMissMat(dat2) # unstructured missing data

dat3 <- dat; 
dat3[80:100,6:10] <- NA 
plotMissMat(dat3) # structured missing data

# let's use the [state.x77](https://gexijin.github.io/learnR/the-state-dataset.html) dataset in R for our example.
# It shows state-wise information of an assortment of variables that affect quality of life.

# ggplot has informative [tutorial website](https://ggplot2.tidyverse.org/)
# with [cheatsheets](https://github.com/rstudio/cheatsheets/blob/main/data-visualization.pdf)

# create a scatterplot from two continuous variables
require(ggplot2)
head(state.x77)
dim(state.x77)
dat <- as.data.frame(state.x77)
p <- ggplot(dat, aes(Illiteracy, Income)) + geom_point()
p

# add a confidence interval
p + geom_smooth()

# test a model
mod <- lm(Income~Illiteracy,dat)
summary(mod)

# for categorical data let's look at the fuel economy dataset
p <- ggplot(diamonds, aes(clarity, price)) + geom_boxplot()
p

# now let's look at bar plots nad grouping
mpg <- ggplot2::mpg
p <- ggplot(mpg, aes(class))
p + geom_bar()

# in each category show the split of a second variable
# here, for each class, show the split between different cylinder vehicles
# notice how the fill color 
p + geom_bar(aes(fill = factor(cyl)))

# increase font size for presentations and publications
# use `theme`
p <- p + geom_bar(aes(fill = factor(cyl)))
p + theme(text = element_text(size=20)
)

# looking at dichotomous variable

# linear models assume noise distribution is Gaussian (bell curve).
# A common class of models is those with binary outcomes. These are modelled using logistic regression. They are also used to model Poisson distributed functions. 

# https://vitalflux.com/generalized-linear-models-explained-with-examples/

# Instead of fitting a line, logistic regression fits a sigmoid.


# Examples: binary variables, such as survival or treatment response
# source: https://github.com/jeffprosise/Machine-Learning/blob/master/Data/titanic.csv
dat  <- read.delim("titanic.csv",sep=",")

p <- ggplot(dat, aes(Survived)) 
p + geom_bar()
p + geom_bar(aes(fill = factor(Sex)))

mod <- glm(
    Survived ~ Sex + Pclass + Embarked,
    data = dat,
    family = "binomial")

summary(mod)
plotMissMat(dat)

# add up the number of NAs in each column of data.
x <- colSums(is.na(dat))
# plot it.
barplot(x,main="Number of missing values, Titanic data")

# apply this to the Pima Indians Diabetes dataset
#install.packages(mlbench)
dat <- PimaIndiansDiabetes2
# plot missingness
x <- colSums(is.na(dat))
barplot(x, las=3)
p <- ggplot(dat, aes(glucose)) + geom_boxplot(aes(fill=diabetes))
plotMissMat(dat)
mod <- glm(diabetes ~ glucose + pregnant + age, dat, family="binomial")
summary(mod)
