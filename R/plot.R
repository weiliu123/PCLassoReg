#' Plot coefficients from a PCLasso object
#'
#' @description Produces a plot of the coefficient paths for a fitted
#' \code{PCLasso} object.
#' @param x Fitted \code{PCLasso} model.
#' @param norm If TRUE, plot the norm of each group, rather than the individual
#' coefficients.
#' @param ... Other graphical parameters to \code{plot}.
#'
#' @seealso \code{\link{PCLasso}}
#' @export
#'
#' @examples
#' # load data
#' data(survivalData)
#' data(PCGroups)
#'
#' x = survivalData$Exp
#' y = survivalData$survData
#'
#' PC.Human <- getPCGroups(Groups = PCGroups, Organism = "Human",
#' Type = "EntrezID")
#'
#' # fit PCLasso model
#' fit.PCLasso <- PCLasso(x, y, group = PC.Human, penalty = "grLasso")
#'
#' # plot the norm of each group
#' plot(fit.PCLasso, norm = TRUE)
#'
#' # plot the individual coefficients
#' plot(fit.PCLasso, norm = FALSE)
plot.PCLasso <-
function(x, norm = TRUE, ...){
    plot(x$fit, norm = norm, ...)
}

#' Plot the cross-validation curve from a \code{cv.PCLasso} object
#'
#' @description Plot the cross-validation curve from a \code{cv.PCLasso} object,
#' along with standard error bars.
#' @param x Fitted \code{cv.PCLasso} model.
#' @param type What to plot on the vertical axis. "cve" plots the
#' cross-validation error (deviance); "rsq" plots an estimate of the fraction of
#' the deviance explained by the model (R-squared); "snr" plots an estimate of
#' the signal-to-noise ratio; "all" produces all of the above.
#' @param norm If TRUE, plot the norm of each group, rather than the individual
#' coefficients.
#' @param ... Other graphical parameters to \code{plot}
#' @details Error bars representing approximate +/- 1 SE (68\% confidence
#' intervals) are plotted along with the estimates at value of lambda. See
#' \code{plot.cv.grpreg} in the R package \code{grpreg} for details.
#'
#' @export
#' @seealso \code{\link{cv.PCLasso}}
#' @examples
#' # load data
#' data(survivalData)
#' data(PCGroups)
#'
#' x = survivalData$Exp
#' y = survivalData$survData
#'
#' PC.Human <- getPCGroups(Groups = PCGroups, Organism = "Human",
#' Type = "EntrezID")
#'
#' # fit model
#' cv.fit1 <- cv.PCLasso(x, y, group = PC.Human, penalty = "grLasso",
#' nfolds = 10)
#'
#' # plot the norm of each group
#' plot(cv.fit1, norm = TRUE)
#'
#' # plot the individual coefficients
#' plot(cv.fit1, norm = FALSE)
#'
#' # plot the cross-validation error (deviance)
#' plot(cv.fit1, type = "cve")
plot.cv.PCLasso <-
    function(x, type = c("cve", "rsq", "snr", "all"),
             norm = NULL, ...){
        if(is.null(norm)){
            type <- match.arg(type)
            plot(x$cv.fit, type = type, ...)
        }else{
            plot(x$cv.fit$fit, norm = norm, ...)
        }
    }


#' Plot coefficients from a PCLasso2 object
#'
#' @description Produces a plot of the coefficient paths for a fitted
#' \code{PCLasso2} object.
#' @param x Fitted \code{PCLasso2} model.
#' @param norm If TRUE, plot the norm of each group, rather than the individual
#' coefficients.
#' @param ... Other graphical parameters to \code{plot}.
#'
#' @seealso \code{\link{PCLasso2}}
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
#' # fit PCLasso2 model
#' fit.PCLasso2 <- PCLasso2(x, y, group = PC.Human, penalty = "grLasso")
#'
#' # plot the norm of each group
#' plot(fit.PCLasso2, norm = TRUE)
#'
#' # plot the individual coefficients
#' plot(fit.PCLasso2, norm = FALSE)
plot.PCLasso2 <-
    function(x, norm = TRUE, ...){
        plot(x$fit, norm = norm, ...)
    }

#' Plot the cross-validation curve from a \code{cv.PCLasso2} object
#'
#' @description Plot the cross-validation curve from a \code{cv.PCLasso2}
#' object, along with standard error bars.
#' @param x Fitted \code{cv.PCLasso2} model.
#' @param type What to plot on the vertical axis. "cve" plots the
#' cross-validation error (deviance); "rsq" plots an estimate of the fraction of
#' the deviance explained by the model (R-squared); "snr" plots an estimate of
#' the signal-to-noise ratio; "all" produces all of the above.
#' @param norm If TRUE, plot the norm of each group, rather than the individual
#' coefficients.
#' @param ... Other graphical parameters to \code{plot}
#' @details Error bars representing approximate +/- 1 SE (68\% confidence
#' intervals) are plotted along with the estimates at value of lambda. See
#' \code{plot.cv.grpreg} in the R package \code{grpreg} for details.
#'
#' @export
#' @seealso \code{\link{cv.PCLasso2}}
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
#'
#' # plot the norm of each group
#' plot(cv.fit1, norm = TRUE)
#'
#' # plot the individual coefficients
#' plot(cv.fit1, norm = FALSE)
#'
#' # plot the cross-validation error (deviance)
#' plot(cv.fit1, type = "cve")
plot.cv.PCLasso2 <-
    function(x, type = c("cve", "rsq", "snr", "all"),
             norm = NULL, ...){
        if(is.null(norm)){
            type <- match.arg(type)
            plot(x$cv.fit, type = type, ...)
        }else{
            plot(x$cv.fit$fit, norm = norm, ...)
        }
    }

