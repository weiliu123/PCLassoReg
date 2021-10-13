#' Cross-validation for \code{PCLasso2}
#'
#' Perform k-fold cross validations for the PCLasso2 model with grouped
#' covariates over a grid of values for the regularization parameter
#' \code{lambda}.
#'
#' @param x A n x p design matrix of gene expression measurements with n samples
#'   and p genes, as in \code{PCLasso2}.
#' @param y The response vector.
#' @param group A list of groups as in \code{PCLasso}. The feature (gene) names
#'   in \code{group} should be consistent with the feature (gene) names in
#'   \code{x}.
#' @param penalty The penalty to be applied to the model. For group selection,
#' one of grLasso, grMCP, or grSCAD. See \code{grpreg} in the R package
#' \code{grpreg} for details.
#' @param family Either "binomial" or "gaussian", depending on the response.
#' @param nfolds The number of cross-validation folds. Default is 5.
#' @param standardize Logical flag for \code{x} standardization, prior to
#'   fitting the model. Default is \code{TRUE}.
#' @param ... Arguments to be passed to \code{cv.grpreg} in the R package
#'   \code{grpreg}.
#'
#' @details
#' The function calls \code{PCLasso2} \code{nfolds} times, each time leaving out
#' 1/\code{nfolds} of the data. The cross-validation error is based on the
#' deviance. The numbers for each class are balanced across the folds;
#' i.e., the number of outcomes in which y is equal to 1 is the same for each
#' fold, or possibly off by 1 if the numbers do not divide evenly. See
#' \code{cv.grpreg} in the R package \code{grpreg} for details.
#'
#' @return An object with S3 class "cv.PCLasso2" containing:
#' \item{cv.fit}{
#' An object of class "cv.grpreg".}
#' \item{complexes.dt}{
#'   Complexes with  features (proteins) not included in \code{x} being filtered
#'   out. }
#' @import grpreg
#' @export
#'
#' @examples
#' # load data
#' data(classData)
#' data(PCGroups)
#'
#' x = classData$Exp
#' y = classData$Label
#'
#' PC.Human <- getPCGroups(Groups = PCGroups, Organism = "Human",
#' Type = "GeneSymbol")
#'
#' # fit model
#' cv.fit1 <- cv.PCLasso2(x, y, group = PC.Human, penalty = "grLasso",
#' family = "binomial", nfolds = 10)
#' cv.fit1 <- cv.PCLasso2(x, y, group = PC.Human, penalty = "grSCAD",
#' family = "binomial", nfolds = 5, gamma = 10)
#' cv.fit1 <- cv.PCLasso2(x, y, group = PC.Human, penalty = "grMCP",
#' family = "binomial", nfolds = 10, gamma = 15)
#' @references
#' PCLasso2: a protein complex-based, group Lasso-logistic model for risk
#' protein complex discovery. To be published.
#'
#' Park, H., Niida, A., Miyano, S. and Imoto, S. (2015) Sparse overlapping group
#' lasso for integrative multi-omics analysis. Journal of computational biology:
#' a journal of computational molecular cell biology, 22, 73-84.
#' @author Wei Liu
#' @seealso \code{\link{predict.cv.PCLasso2}}

cv.PCLasso2 <-
    function(x, y, group,
             penalty = c("grLasso", "grMCP", "grSCAD"),
             family=c("binomial", "gaussian", "poisson"),
             nfolds = 5,
             standardize = TRUE,...){
        penalty <- match.arg(penalty)
        family = match.arg(family)

        if(standardize){
            x <- scale(x, center = TRUE, scale = TRUE)
        }

        # feature set in all groups
        featureSet <- unique(unlist(group))

        # common features in groups and expression matrix x
        commonFeat <- intersect(colnames(x), featureSet)

        # filter undetected genes in expression matrix x
        x <- x[,commonFeat]

        # filter undetected genes in groups
        # Construct groups whose expressions are detected
        group.dt <- vector(mode = "list", length = 0)
        idx <- 0
        for(i in 1:length(group)){
            group.i <- intersect(group[[i]], commonFeat)
            if(length(group.i) > 1){
                idx <- idx + 1
                group.dt[[idx]] <- group.i
                names(group.dt)[idx] <- names(group)[i]
            }
        }

        # Filter duplicate groups (generated due to undetected genes)
        group.dt <- group.dt[!duplicated(group.dt)]

        # extended genes
        commonFeat.ext <- unlist(group.dt)

        # New names of extended genes
        # The new name consists of "group+_+gene name"
        commonFeat.extName <- c()
        for(i in 1:length(group.dt)){
            names.i <- paste0(names(group.dt)[i], "_", group.dt[[i]])
            commonFeat.extName <- c(commonFeat.extName, names.i)
        }

        # group of extended genes
        groupOfFeats <- c()
        for(i in 1:length(group.dt)){
            group.i <- rep(names(group.dt)[i], length = length(group.dt[[i]]))
            groupOfFeats <- c(groupOfFeats, group.i)
        }

        # extended dataset
        x.ext <- x[, commonFeat.ext]
        colnames(x.ext) <- commonFeat.extName

        # grpsurv
        cv.fit <- grpreg::cv.grpreg(X=x.ext,
                                    y=y,
                                    group = groupOfFeats,
                                    penalty = penalty,
                                    nfolds = nfolds,
                                    family = family,...)

        res <- list(cv.fit = cv.fit, complexes.dt = group.dt)
        class(res) <- "cv.PCLasso2"

        return(res)
    }
