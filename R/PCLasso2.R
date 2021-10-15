#' Protein complex-based group Lasso-logistic model
#'
#' @param x A n x p matrix of gene/protein expression measurements with n samples
#' and p genes/proteins.
#' @param y The response vector.
#' @param group A list of groups. The feature (gene/protein) names in \code{group}
#' should be consistent with the feature (gene/protein) names in \code{x}.
#' @param penalty The penalty to be applied to the model. For group selection,
#' one of grLasso, grMCP, or grSCAD. See \code{grpreg} in the R package
#' \code{grpreg} for details.
#' @param family Either "binomial" or "gaussian", depending on the response.
#' @param gamma Tuning parameter of the \code{grSCAD}/\code{grMCP} penalty.
#' Default is 8.
#' @param standardize Logical flag for \code{x} standardization, prior to
#' fitting the model. Default is \code{TRUE}.
#' @param ... Arguments to be passed to \code{grpreg} in the R package
#' \code{grpreg}.
#'
#'@details The PCLasso2 model is a classification model that selects important
#'  predictors at the protein complex level to achieve accurate classification
#'  and identify risk protein complexes. The PCLasso2 model has three inputs: a
#'  protein expression matrix, a vector of binary response variables, and a
#'  number of known protein complexes. It estimates the correlation between
#'  protein expression and response variable at the level of protein complexes.
#'  Similar to traditional Lasso-logistic model, PCLasso2 is based on the
#'  logistic regression model and estimates the logistic regression coefficients
#'  by maximizing likelihood function with regularization penalty. The
#'  difference is that PCLasso2 selects features at the level of protein
#'  complexes rather than individual proteins. Considering that proteins usually
#'  function by forming protein complexes, PCLasso2 regards proteins belonging
#'  to the same protein complex as a group and constructs a group Lasso penalty
#'  (l1/l2 penalty) based on the sum (i.e. l1 norm) of the l2 norms of the
#'  regression coefficients of the group members to perform the selection of
#'  features at the group level. With the group Lasso penalty, PCLasso2 trains
#'  the logistic regression model and obtains a sparse solution at the protein
#'  complex level, that is, the proteins belonging to a protein complex are
#'  either wholly included or wholly excluded from the model. PCLasso2 outputs a
#'  prediction model and a small set of protein complexes included in the model,
#'  which are referred to as risk protein complexes. The PCSCAD and PCMCP are
#'  performed by setting the penalty parameter \code{penalty} as \code{grSCAD}
#'  and \code{grMCP}, respectively.
#' @return An object with S3 class \code{PCLasso2} containing:
#' \item{fit}{An object of class \code{grpreg}}
#' \item{Complexes.dt}{Complexes with  features (genes/proteins) not included
#'     in \code{x} being filtered out. }
#' @import grpreg
#' @export
#' @references
#' PCLasso2: a protein complex-based, group Lasso-logistic model for risk
#' protein complex discovery. To be published.
#'
#' PCLasso: a protein complex-based, group lasso-Cox model for accurate
#' prognosis and risk protein complex discovery. Brief Bioinform, 2021.
#'
#' Park, H., Niida, A., Miyano, S. and Imoto, S. (2015) Sparse overlapping group
#' lasso for integrative multi-omics analysis. Journal of computational biology:
#'     a journal of computational molecular cell biology, 22, 73-84.
#' @seealso \code{\link{predict.PCLasso2}}, \code{\link{cv.PCLasso2}}
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
#' # fit PCLasso2 model
#' fit.PCLasso2 <- PCLasso2(x, y, group = PC.Human, penalty = "grLasso",
#' family = "binomial")
#'
#' # fit PCSCAD model
#' fit.PCSCAD <- PCLasso2(x, y, group = PC.Human, penalty = "grSCAD",
#' family = "binomial", gamma = 10)
#'
#' # fit PCMCP model
#' fit.PCMCP <- PCLasso2(x, y, group = PC.Human, penalty = "grMCP",
#' family = "binomial", gamma = 9)
PCLasso2 <- function(x, y, group,
                    penalty = c("grLasso", "grMCP", "grSCAD"),
                    family=c("binomial", "gaussian", "poisson"),
                    gamma = 8,
                    standardize = TRUE,...){

    penalty = match.arg(penalty)
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

    Complexes.dt <- vector(mode = "list", length = 0)
    idx <- 0
    for(i in 1:length(group)){
        Complexes.i <- intersect(group[[i]], commonFeat)
        if(length(Complexes.i) > 1){
            idx <- idx + 1
            Complexes.dt[[idx]] <- Complexes.i
            names(Complexes.dt)[idx] <- names(group)[i]
        }
    }

    # 过滤重复的Complex（由于有蛋白未检测到而产生）
    Complexes.dt <- Complexes.dt[!duplicated(Complexes.dt)]

    # extended expression profiles -----------
    # extended genes
    commonFeat.ext <- unlist(Complexes.dt)

    # New names of extended genes
    commonFeat.extName <- c()
    for(i in 1:length(Complexes.dt)){
        names.i <- paste0(names(Complexes.dt)[i], "_", Complexes.dt[[i]])
        commonFeat.extName <- c(commonFeat.extName, names.i)
    }

    # group of extended genes
    groupOfFeats <- c()
    for(i in 1:length(Complexes.dt)){
        group.i <- rep(names(Complexes.dt)[i], length = length(Complexes.dt[[i]]))
        groupOfFeats <- c(groupOfFeats, group.i)
    }

    # groupOfFeats <- factor(groupOfFeats)

    # extended dataset
    x.ext <- x[, commonFeat.ext]
    colnames(x.ext) <- commonFeat.extName

    # grpsurv
    fit <- grpreg::grpreg(X=x.ext,
                   y=y,
                   group=groupOfFeats,
                   penalty = penalty,
                   gamma = gamma,
                   family = family,...)


    res <- list(fit = fit, complexes.dt = Complexes.dt)

    class(res) <- c("PCLasso2")

    return(res)
}
