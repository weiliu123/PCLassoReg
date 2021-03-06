% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict.PCLasso}
\alias{predict.PCLasso}
\title{Make predictions from a PCLasso model}
\usage{
\method{predict}{PCLasso}(
  object,
  x = NULL,
  type = c("link", "response", "survival", "median", "norm", "coefficients", "vars",
    "nvars", "vars.unique", "nvars.unique", "groups", "ngroups"),
  lambda,
  ...
)
}
\arguments{
\item{object}{Fitted \code{PCLasso} model object.}

\item{x}{Matrix of values at which predictions are to be made. The features
(genes/proteins) contained in \code{x} should be consistent with those
contained in \code{x} in the \code{PCLasso} function.  Not used for
type="coefficients" or for some of the type settings in \code{predict}.}

\item{type}{Type of prediction: "link" returns the linear predictors;
"response" gives the risk (i.e., exp(link)); "vars" returns the indices for
the nonzero coefficients; "vars.unique" returns unique features
(genes/proteins) with nonzero coefficients (If a feature belongs to
multiple groups and multiple groups are selected, the feature will be
repeatedly selected. Compared with "var", "var.unique" will filter out
repeated features.); "groups" returns the groups with at least one nonzero
coefficient; "nvars" returns the number of nonzero coefficients;
"nvars.unique" returns the number of unique features (genes/proteins) with
nonzero coefficients; "ngroups" returns the number of groups with at least
one nonzero coefficient; "norm" returns the L2 norm of the coefficients in
each group."survival" returns the estimated survival function; "median"
estimates median survival times.}

\item{lambda}{Values of the regularization parameter \code{lambda} at which
predictions are requested. For values of \code{lambda} not in the sequence
of fitted models, linear interpolation is used.}

\item{...}{Arguments to be passed to \code{predict.grpsurv} in the R package
\code{grpreg}.}
}
\value{
The object returned depends on \code{type}.
}
\description{
Similar to other predict methods, this function returns
predictions from a fitted \code{PCLasso} object.
}
\details{
See \code{predict.grpsurv} in the R package \code{grpreg} for
details.
}
\examples{
# load data
data(survivalData)
data(PCGroups)

x <- survivalData$Exp
y <- survivalData$survData
PC.Human <- getPCGroups(Groups = PCGroups, Organism = "Human",
Type = "EntrezID")

set.seed(20150122)
idx.train <- sample(nrow(x), round(nrow(x)*2/3))
x.train <- x[idx.train,]
y.train <- y[idx.train,]
x.test <- x[-idx.train,]
y.test <- y[-idx.train,]

# fit PCLasso model
fit.PCLasso <- PCLasso(x = x.train, y = y.train, group = PC.Human,
                  penalty = "grLasso")

# predict risk scores of samples in x.test
s <- predict(object = fit.PCLasso, x = x.test, type="link",
lambda=fit.PCLasso$fit$lambda)

s <- predict(object = fit.PCLasso, x = x.test, type="link",
lambda=fit.PCLasso$fit$lambda[10])

# Nonzero coefficients
sel.groups <- predict(object = fit.PCLasso, type="groups",
                      lambda = fit.PCLasso$fit$lambda)
sel.ngroups <- predict(object = fit.PCLasso, type="ngroups",
                       lambda = fit.PCLasso$fit$lambda)
sel.vars.unique <- predict(object = fit.PCLasso, type="vars.unique",
                          lambda = fit.PCLasso$fit$lambda)
sel.nvars.unique <- predict(object = fit.PCLasso, type="nvars.unique",
                            lambda = fit.PCLasso$fit$lambda)
sel.vars <- predict(object = fit.PCLasso, type="vars",
                    lambda=fit.PCLasso$fit$lambda)
sel.nvars <- predict(object = fit.PCLasso, type="nvars",
                     lambda=fit.PCLasso$fit$lambda)

# For values of lambda not in the sequence of fitted models,
# linear interpolation is used.
sel.groups <- predict(object = fit.PCLasso, type="groups",
                      lambda = c(0.1, 0.05))
sel.ngroups <- predict(object = fit.PCLasso, type="ngroups",
                       lambda = c(0.1, 0.05))
sel.vars.unique <- predict(object = fit.PCLasso, type="vars.unique",
                           lambda = c(0.1, 0.05))
sel.nvars.unique <- predict(object = fit.PCLasso, type="nvars.unique",
                            lambda = c(0.1, 0.05))
sel.vars <- predict(object = fit.PCLasso, type="vars",
                    lambda=c(0.1, 0.05))
sel.nvars <- predict(object = fit.PCLasso, type="nvars",
                     lambda=c(0.1, 0.05))

}
\seealso{
\code{\link{PCLasso}}
}
