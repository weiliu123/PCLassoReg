% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.PCLasso2.R
\name{cv.PCLasso2}
\alias{cv.PCLasso2}
\title{Cross-validation for \code{PCLasso2}}
\usage{
cv.PCLasso2(
  x,
  y,
  group,
  penalty = c("grLasso", "grMCP", "grSCAD"),
  family = c("binomial", "gaussian", "poisson"),
  nfolds = 5,
  gamma = 8,
  standardize = TRUE,
  ...
)
}
\arguments{
\item{x}{A n x p design matrix of gene/protein expression measurements with n
samples and p genes/proteins, as in \code{PCLasso2}.}

\item{y}{The response vector.}

\item{group}{A list of groups as in \code{PCLasso}. The feature
(gene/protein) names in \code{group} should be consistent with the feature
(gene/protein) names in \code{x}.}

\item{penalty}{The penalty to be applied to the model. For group selection,
one of grLasso, grMCP, or grSCAD. See \code{grpreg} in the R package
\code{grpreg} for details.}

\item{family}{Either "binomial" or "gaussian", depending on the response.}

\item{nfolds}{The number of cross-validation folds. Default is 5.}

\item{gamma}{Tuning parameter of the \code{grSCAD}/\code{grMCP} penalty.
Default is 8.}

\item{standardize}{Logical flag for \code{x} standardization, prior to
fitting the model. Default is \code{TRUE}.}

\item{...}{Arguments to be passed to \code{cv.grpreg} in the R package
\code{grpreg}.}
}
\value{
An object with S3 class "cv.PCLasso2" containing: \item{cv.fit}{ An
object of class "cv.grpreg".} \item{complexes.dt}{ Complexes with  features
(genes/proteins) not included in \code{x} being filtered out. }
}
\description{
Perform k-fold cross validations for the PCLasso2 model with grouped
covariates over a grid of values for the regularization parameter
\code{lambda}.
}
\details{
The function calls \code{PCLasso2} \code{nfolds} times, each time
leaving out 1/\code{nfolds} of the data. The cross-validation error is based
on the deviance. The numbers for each class are balanced across the folds;
i.e., the number of outcomes in which y is equal to 1 is the same for each
fold, or possibly off by 1 if the numbers do not divide evenly. See
\code{cv.grpreg} in the R package \code{grpreg} for details.
}
\examples{
# load data
data(classData)
data(PCGroups)

x = classData$Exp
y = classData$Label

PC.Human <- getPCGroups(Groups = PCGroups, Organism = "Human",
Type = "GeneSymbol")

# fit model
cv.fit1 <- cv.PCLasso2(x, y, group = PC.Human, penalty = "grLasso",
family = "binomial", nfolds = 5)

\donttest{
cv.fit1 <- cv.PCLasso2(x, y, group = PC.Human, penalty = "grSCAD",
family = "binomial", nfolds = 5, gamma = 10)
cv.fit1 <- cv.PCLasso2(x, y, group = PC.Human, penalty = "grMCP",
family = "binomial", nfolds = 5, gamma = 15)}

}
\references{
PCLasso2: a protein complex-based, group Lasso-logistic model for risk
protein complex discovery. To be published.

Park, H., Niida, A., Miyano, S. and Imoto, S. (2015) Sparse overlapping group
lasso for integrative multi-omics analysis. Journal of computational biology:
a journal of computational molecular cell biology, 22, 73-84.
}
\seealso{
\code{\link{predict.cv.PCLasso2}}
}
\author{
Wei Liu
}
