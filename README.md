
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PCLassoReg

<!-- badges: start -->
<!-- badges: end -->

Protein complex-based group regression models for risk protein complex
identification.

## Installation

To install the released version from [CRAN](https://CRAN.R-project.org)
with:

``` r
install.packages("PCLassoReg")
```

To install the latest development version from GitHub:

``` r
devtools::install_github("weiliu123/PCLassoReg")
```

## Details

The package implements protein complex-based group regression models
(PCLasso and PCLasso2) for risk protein complex identification. PCLasso
is a prognostic model that identifies risk protein complexes associated
with survival. PCLasso2 is a classification model that identifies risk
protein complexes associated with classes.

## Example

``` r
library(PCLassoReg)

##### PCLasso #####
# load data
data(survivalData)
data(PCGroups)

x <- survivalData$Exp
y <- survivalData$survData

# get human protein complexes
PC.Human <- getPCGroups(Groups = PCGroups, Organism = "Human",
                        Type = "EntrezID")

set.seed(20150122)
idx.train <- sample(nrow(x), round(nrow(x)*2/3))
x.train <- x[idx.train,]
y.train <- y[idx.train,]
x.test <- x[-idx.train,]
y.test <- y[-idx.train,]

# fit cv.PCLasso model
cv.fit1 <- cv.PCLasso(x = x.train, y = y.train, group = PC.Human, nfolds = 5)

# predict risk scores of samples in x.test
s <- predict(object = cv.fit1, x = x.test, type="link",
             lambda=cv.fit1$cv.fit$lambda.min)

# Nonzero coefficients
sel.groups <- predict(object = cv.fit1, type="groups",
                      lambda = cv.fit1$cv.fit$lambda.min)
sel.vars.unique <- predict(object = cv.fit1, type="vars.unique",
                           lambda = cv.fit1$cv.fit$lambda.min)

##### PCLasso2 #####
# load data
data(classData)
data(PCGroups)

x = classData$Exp
y = classData$Label

# get human protein complexes
PC.Human <- getPCGroups(Groups = PCGroups, Organism = "Human",
                        Type = "GeneSymbol")

set.seed(20150122)
idx.train <- sample(nrow(x), round(nrow(x)*2/3))
x.train <- x[idx.train,]
y.train <- y[idx.train]
x.test <- x[-idx.train,]
y.test <- y[-idx.train]

# fit model
cv.fit1 <- cv.PCLasso2(x = x.train, y = y.train, group = PC.Human,
                       penalty = "grLasso", family = "binomial", nfolds = 5)

# predict risk scores of samples in x.test
s <- predict(object = cv.fit1, x = x.test, type="class",
             lambda=cv.fit1$cv.fit$lambda.min)

# Nonzero coefficients
sel.groups <- predict(object = cv.fit1, type="groups",
                      lambda = cv.fit1$cv.fit$lambda.min)
sel.vars.unique <- predict(object = cv.fit1, type="vars.unique",
                           lambda = cv.fit1$cv.fit$lambda.min)
```

# Reference

PCLasso2: a protein complex-based, group Lasso-logistic model for risk
protein complex discovery. To be published.

PCLasso: a protein complex-based, group lasso-Cox model for accurate
prognosis and risk protein complex discovery. Brief Bioinform, 2021.

Park, H., Niida, A., Miyano, S. and Imoto, S. (2015) Sparse overlapping
group lasso for integrative multi-omics analysis. Journal of
computational biology: a journal of computational molecular cell
biology, 22, 73-84.
