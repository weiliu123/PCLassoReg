% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCLasso2.R
\name{PCLasso2}
\alias{PCLasso2}
\title{Protein complex-based group Lasso-logistic model}
\usage{
PCLasso2(
  x,
  y,
  group,
  penalty = c("grLasso", "grMCP", "grSCAD"),
  family = c("binomial", "gaussian", "poisson"),
  gamma = 8,
  standardize = TRUE,
  ...
)
}
\arguments{
\item{x}{A n x p matrix of gene/protein expression measurements with n samples
and p genes/proteins.}

\item{y}{The response vector.}

\item{group}{A list of groups. The feature (gene/protein) names in \code{group}
should be consistent with the feature (gene/protein) names in \code{x}.}

\item{penalty}{The penalty to be applied to the model. For group selection,
one of grLasso, grMCP, or grSCAD. See \code{grpreg} in the R package
\code{grpreg} for details.}

\item{family}{Either "binomial" or "gaussian", depending on the response.}

\item{gamma}{Tuning parameter of the \code{grSCAD}/\code{grMCP} penalty.
Default is 8.}

\item{standardize}{Logical flag for \code{x} standardization, prior to
fitting the model. Default is \code{TRUE}.}

\item{...}{Arguments to be passed to \code{grpreg} in the R package
\code{grpreg}.}
}
\value{
An object with S3 class \code{PCLasso2} containing:
\item{fit}{An object of class \code{grpreg}}
\item{Complexes.dt}{Complexes with  features (genes/proteins) not included
in \code{x} being filtered out. }
}
\description{
Protein complex-based group Lasso-logistic model
}
\details{
The PCLasso2 model is a classification model that selects important
predictors at the protein complex level to achieve accurate classification
and identify risk protein complexes. The PCLasso2 model has three inputs: a
protein expression matrix, a vector of binary response variables, and a
number of known protein complexes. It estimates the correlation between
protein expression and response variable at the level of protein complexes.
Similar to traditional Lasso-logistic model, PCLasso2 is based on the
logistic regression model and estimates the logistic regression coefficients
by maximizing likelihood function with regularization penalty. The
difference is that PCLasso2 selects features at the level of protein
complexes rather than individual proteins. Considering that proteins usually
function by forming protein complexes, PCLasso2 regards proteins belonging
to the same protein complex as a group and constructs a group Lasso penalty
(l1/l2 penalty) based on the sum (i.e. l1 norm) of the l2 norms of the
regression coefficients of the group members to perform the selection of
features at the group level. With the group Lasso penalty, PCLasso2 trains
the logistic regression model and obtains a sparse solution at the protein
complex level, that is, the proteins belonging to a protein complex are
either wholly included or wholly excluded from the model. PCLasso2 outputs a
prediction model and a small set of protein complexes included in the model,
which are referred to as risk protein complexes. The PCSCAD and PCMCP are
performed by setting the penalty parameter \code{penalty} as \code{grSCAD}
and \code{grMCP}, respectively.
}
\examples{
# load data
data(classData)
data(PCGroups)

x = classData$Exp
y = classData$Label

PC.Human <- getPCGroups(Groups = PCGroups, Organism = "Human",
Type = "GeneSymbol")

# fit PCLasso2 model
fit.PCLasso2 <- PCLasso2(x, y, group = PC.Human, penalty = "grLasso",
family = "binomial")

# fit PCSCAD model
fit.PCSCAD <- PCLasso2(x, y, group = PC.Human, penalty = "grSCAD",
family = "binomial", gamma = 10)

# fit PCMCP model
fit.PCMCP <- PCLasso2(x, y, group = PC.Human, penalty = "grMCP",
family = "binomial", gamma = 9)
}
\references{
PCLasso2: a protein complex-based, group Lasso-logistic model for risk
protein complex discovery. To be published.

PCLasso: a protein complex-based, group lasso-Cox model for accurate
prognosis and risk protein complex discovery. Brief Bioinform, 2021.

Park, H., Niida, A., Miyano, S. and Imoto, S. (2015) Sparse overlapping group
lasso for integrative multi-omics analysis. Journal of computational biology:
a journal of computational molecular cell biology, 22, 73-84.
}
\seealso{
\code{\link{predict.PCLasso2}}, \code{\link{cv.PCLasso2}}
}
