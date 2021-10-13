#' Make predictions from a PCLasso model
#'
#' @description Similar to other predict methods, this function returns
#' predictions from a fitted \code{PCLasso} object.
#'
#' @param object Fitted \code{PCLasso} model object.
#' @param x Matrix of values at which predictions are to be made. The features
#' (genes) contained in \code{x} should be consistent with those contained in
#' \code{x} in the \code{PCLasso} function.  Not used for type="coefficients"
#' or for some of the type settings in \code{predict}.
#' @param type Type of prediction: "link" returns the linear predictors;
#'   "response" gives the risk (i.e., exp(link)); "vars" returns the indices for
#'   the nonzero coefficients; "vars.unique" returns unique features
#'   (genes) with nonzero coefficients (If a feature belongs to multiple groups
#'   and multiple groups are selected, the feature will be repeatedly selected.
#'   Compared with "var", "var.unique" will filter out repeated features.);
#'   "groups" returns the groups with at least one nonzero coefficient; "nvars"
#'   returns the number of nonzero coefficients; "nvars.unique" returns the
#'   number of unique features (genes) with nonzero coefficients; "ngroups"
#'   returns the number of groups with at least one nonzero coefficient; "norm"
#'   returns the L2 norm of the coefficients in each group."survival" returns
#'   the estimated survival function; "median" estimates median survival times.
#' @param lambda Values of the regularization parameter \code{lambda} at which
#'   predictions are requested. For values of \code{lambda} not in the sequence
#'   of fitted models, linear interpolation is used.
#' @param ... Arguments to be passed to \code{predict.grpsurv} in the R package
#' \code{grpreg}.
#' @details
#' See \code{predict.grpsurv} in the R package \code{grpreg} for details.
#' @return The object returned depends on \code{type}.
#' @seealso \code{\link{PCLasso}}
#' @importFrom stats predict
#' @export
#'
#' @examples
#' # load data
#' data(survivalData)
#' data(PCGroups)
#'
#' x <- survivalData$Exp
#' y <- survivalData$survData
#' PC.Human <- getPCGroups(Groups = PCGroups, Organism = "Human",
#' Type = "EntrezID")
#'
#' set.seed(20150122)
#' idx.train <- sample(nrow(x), round(nrow(x)*2/3))
#' x.train <- x[idx.train,]
#' y.train <- y[idx.train,]
#' x.test <- x[-idx.train,]
#' y.test <- y[-idx.train,]
#'
#' # fit PCLasso model
#' fit.PCLasso <- PCLasso(x = x.train, y = y.train, group = PC.Human,
#'                   penalty = "grLasso")
#'
#' # predict risk scores of samples in x.test
#' s <- predict(object = fit.PCLasso, x = x.test, type="link",
#' lambda=fit.PCLasso$fit$lambda)
#'
#' s <- predict(object = fit.PCLasso, x = x.test, type="link",
#' lambda=fit.PCLasso$fit$lambda[10])
#'
#' # Nonzero coefficients
#' sel.groups <- predict(object = fit.PCLasso, type="groups",
#'                       lambda = fit.PCLasso$fit$lambda)
#' sel.ngroups <- predict(object = fit.PCLasso, type="ngroups",
#'                        lambda = fit.PCLasso$fit$lambda)
#' sel.vars.unique <- predict(object = fit.PCLasso, type="vars.unique",
#'                           lambda = fit.PCLasso$fit$lambda)
#' sel.nvars.unique <- predict(object = fit.PCLasso, type="nvars.unique",
#'                             lambda = fit.PCLasso$fit$lambda)
#' sel.vars <- predict(object = fit.PCLasso, type="vars",
#'                     lambda=fit.PCLasso$fit$lambda)
#' sel.nvars <- predict(object = fit.PCLasso, type="nvars",
#'                      lambda=fit.PCLasso$fit$lambda)
#'
#' # For values of lambda not in the sequence of fitted models,
#' # linear interpolation is used.
#' sel.groups <- predict(object = fit.PCLasso, type="groups",
#'                       lambda = c(0.1, 0.05))
#' sel.ngroups <- predict(object = fit.PCLasso, type="ngroups",
#'                        lambda = c(0.1, 0.05))
#' sel.vars.unique <- predict(object = fit.PCLasso, type="vars.unique",
#'                            lambda = c(0.1, 0.05))
#' sel.nvars.unique <- predict(object = fit.PCLasso, type="nvars.unique",
#'                             lambda = c(0.1, 0.05))
#' sel.vars <- predict(object = fit.PCLasso, type="vars",
#'                     lambda=c(0.1, 0.05))
#' sel.nvars <- predict(object = fit.PCLasso, type="nvars",
#'                      lambda=c(0.1, 0.05))
#'
predict.PCLasso <-
    function(object, x = NULL,
             type = c("link", "response", "survival", "median", "norm", "coefficients",
                      "vars", "nvars","vars.unique", "nvars.unique", "groups", "ngroups"),
             lambda, ...){

        type <- match.arg(type)
        if(type == "vars.unique"){
            vars.tmp <- predict(object = object$fit,
                                type = "vars", lambda = lambda, ...)

            if(is.list(vars.tmp)){
                vars.list <- vector(mode = "list", length = length(vars.tmp))
                names(vars.list) <- names(vars.tmp)
                for(vars.list.i in 1:length(vars.tmp)){
                    if(length(vars.tmp[[vars.list.i]]) > 0){
                        vars.list[[vars.list.i]] <-
                            unique(ext2GeneID(rownames(object$fit$beta)[vars.tmp[[vars.list.i]]]))
                    }else{
                        vars.list[[vars.list.i]] <- vars.tmp[[vars.list.i]]
                    }

                }
                vars.list
            }else{if(length(lambda) > 1){
                vars.vector <- rep(NA, length = length(vars.tmp))
                names(vars.vector) <- names(vars.tmp)
                for(ii in 1:length(vars.tmp)){
                    vars.vector[ii] <-
                        ext2GeneID(rownames(object$fit$beta)[vars.tmp[ii]])
                }
                vars.vector
            }else{
                unique(ext2GeneID(rownames(object$fit$beta)[vars.tmp]))
            }
            }

        }else if(type == "nvars.unique"){
            vars.tmp <- predict(object = object$fit,
                                type = "vars", lambda = lambda, ...)

            if(is.list(vars.tmp)){
                vars.list <- vector(mode = "list", length = length(vars.tmp))
                names(vars.list) <- names(vars.tmp)
                nvars.vector <- rep(0, length = length(vars.tmp))
                names(nvars.vector) <- names(vars.tmp)
                for(vars.list.i in 1:length(vars.tmp)){
                    if(length(vars.tmp[[vars.list.i]]) > 0){
                        vars.list[[vars.list.i]] <-
                            unique(ext2GeneID(rownames(object$fit$beta)[vars.tmp[[vars.list.i]]]))
                        nvars.vector[vars.list.i] <-
                            length(vars.list[[vars.list.i]])
                    }
                }
                nvars.vector
            }else{
                if(length(lambda) > 1){
                    nvars.vector <- rep(NA, length = length(vars.tmp))
                    names(nvars.vector) <- names(vars.tmp)
                    for(ii in 1:length(vars.tmp)){
                        nvars.vector[ii] <-
                            length(ext2GeneID(rownames(object$fit$beta)[vars.tmp[ii]]))
                    }
                    nvars.vector
                }else{
                    length(unique(ext2GeneID(rownames(object$fit$beta)[vars.tmp])))
                }

            }


        }else{
            if(is.null(x)){
                predict(object = object$fit, type = type,
                        lambda = lambda, ...)
            }else{
                # extended genes
                commonFeat.ext <- unlist(object$complexes.dt)

                # New names of extended genes
                # The new name consists of "complexes+_+gene name"
                commonFeat.extName <- c()
                for(i in 1:length(object$complexes.dt)){
                    names.i <- paste0(names(object$complexes.dt)[i], "_",
                                      object$complexes.dt[[i]])
                    commonFeat.extName <- c(commonFeat.extName, names.i)
                }

                # extended dataset
                x.ext <- x[, commonFeat.ext]
                colnames(x.ext) <- commonFeat.extName

                predict(object = object$fit, X = x.ext,
                        type = type, lambda = lambda, ...)
            }
        }
    }



#' Make predictions from a cross-validated PCLasso model
#'
#' @description
#' Similar to other predict methods, this function returns predictions from a
#' fitted "cv.PCLasso" object, using the optimal value chosen for \code{lambda}.
#'
#' @param object Fitted \code{cv.PCLasso} model object.
#' @param x
#' Matrix of values at which predictions are to be made. The features (genes)
#' contained in \code{x} should be consistent with those contained in \code{x}
#' in the \code{cv.PCLasso} function.  Not used for type="coefficients" or for
#' some of the type settings in \code{predict}.
#' @param type
#' Type of prediction: "link" returns the linear predictors; "response" gives
#' the risk (i.e., exp(link)); "vars" returns the indices for the nonzero
#' coefficients; "vars.unique" returns unique features (genes) with nonzero
#' coefficients (If a feature belongs to multiple groups and multiple groups are
#' selected, the feature will be repeatedly selected. Compared with "var",
#' "var.unique" will filter out repeated features.); "groups" returns the groups
#' with at least one nonzero coefficient; "nvars" returns the number of nonzero
#' coefficients; "nvars.unique" returens the number of unique features (genes)
#' with nonzero coefficients; "ngroups" returns the number of groups with at
#' least one nonzero coefficient; "norm" returns the L2 norm of the coefficients
#' in each group."survival" returns the estimated survival function; "median"
#' estimates median survival times.
#' @param lambda
#' Values of the regularization parameter \code{lambda} at which predictions are
#' requested. For values of  \code{lambda} not in the sequence of fitted models,
#' linear interpolation is used.
#' @param ...
#' Arguments to be passed to \code{predict.cv.grpsurv} in the R package
#' \code{grpreg}.
#'
#' @return
#' The object returned depends on \code{type}.
#' @method predict cv.PCLasso
#' @importFrom stats predict
#' @export
#'
#' @seealso \code{\link{cv.PCLasso}}
#'
#' @examples
#' # load data
#' data(survivalData)
#' data(PCGroups)
#'
#' x <- survivalData$Exp
#' y <- survivalData$survData
#' PC.Human <- getPCGroups(Groups = PCGroups, Organism = "Human",
#' Type = "EntrezID")
#'
#' set.seed(20150122)
#' idx.train <- sample(nrow(x), round(nrow(x)*2/3))
#' x.train <- x[idx.train,]
#' y.train <- y[idx.train,]
#' x.test <- x[-idx.train,]
#' y.test <- y[-idx.train,]
#'
#' # fit cv.PCLasso model
#' cv.fit1 <- cv.PCLasso(x = x.train,
#'                       y = y.train,
#'                       group = PC.Human,
#'                       nfolds = 5)
#'
#' # predict risk scores of samples in x.test
#' s <- predict(object = cv.fit1, x = x.test, type="link",
#'              lambda=cv.fit1$cv.fit$lambda.min)
#'
#' # Nonzero coefficients
#' sel.groups <- predict(object = cv.fit1, type="groups",
#'                       lambda = cv.fit1$cv.fit$lambda.min)
#' sel.ngroups <- predict(object = cv.fit1, type="ngroups",
#'                        lambda = cv.fit1$cv.fit$lambda.min)
#' sel.vars.unique <- predict(object = cv.fit1, type="vars.unique",
#'                            lambda = cv.fit1$cv.fit$lambda.min)
#' sel.nvars.unique <- predict(object = cv.fit1, type="nvars.unique",
#'                             lambda = cv.fit1$cv.fit$lambda.min)
#' sel.vars <- predict(object = cv.fit1, type="vars",
#'                     lambda=cv.fit1$cv.fit$lambda.min)
#' sel.nvars <- predict(object = cv.fit1, type="nvars",
#'                      lambda=cv.fit1$cv.fit$lambda.min)

predict.cv.PCLasso <-
    function(object, x = NULL,
             type = c("link", "response", "survival", "median", "norm", "coefficients",
                      "vars", "nvars","vars.unique", "nvars.unique", "groups", "ngroups"),
             lambda, ...){

        type <- match.arg(type)
        if(type == "vars.unique"){

            vars.tmp <- predict(object = object$cv.fit$fit,
                                type = "vars", lambda = lambda, ...)

            if(is.list(vars.tmp)){
                vars.list <- vector(mode = "list", length = length(vars.tmp))
                names(vars.list) <- names(vars.tmp)
                for(vars.list.i in 1:length(vars.tmp)){
                    if(length(vars.tmp[[vars.list.i]]) > 0){
                        vars.list[[vars.list.i]] <-
                            unique(ext2GeneID(rownames(object$cv.fit$fit$beta)[vars.tmp[[vars.list.i]]]))
                    }else{
                        vars.list[[vars.list.i]] <- vars.tmp[[vars.list.i]]
                    }

                }
                vars.list
            }else{
                if(length(lambda) > 1){
                    vars.vector <- rep(NA, length = length(vars.tmp))
                    names(vars.vector) <- names(vars.tmp)
                    for(ii in 1:length(vars.tmp)){
                        vars.vector[ii] <-
                            ext2GeneID(rownames(object$cv.fit$fit$beta)[vars.tmp[ii]])
                    }
                    vars.vector
                }else{
                    unique(ext2GeneID(rownames(object$cv.fit$fit$beta)[vars.tmp]))
                }
            }
        }else if(type == "nvars.unique"){
            vars.tmp <- predict(object = object$cv.fit$fit,
                                type = "vars", lambda = lambda, ...)

            if(is.list(vars.tmp)){
                vars.list <- vector(mode = "list", length = length(vars.tmp))
                names(vars.list) <- names(vars.tmp)
                nvars.vector <- rep(0, length = length(vars.tmp))
                names(nvars.vector) <- names(vars.tmp)
                for(vars.list.i in 1:length(vars.tmp)){
                    if(length(vars.tmp[[vars.list.i]]) > 0){
                        vars.list[[vars.list.i]] <-
                            unique(ext2GeneID(rownames(object$cv.fit$fit$beta)[vars.tmp[[vars.list.i]]]))
                        nvars.vector[vars.list.i] <-
                            length(vars.list[[vars.list.i]])
                    }
                }
                nvars.vector
            }else{
                if(length(lambda) > 1){
                    nvars.vector <- rep(NA, length = length(vars.tmp))
                    names(nvars.vector) <- names(vars.tmp)
                    for(ii in 1:length(vars.tmp)){
                        nvars.vector[ii] <-
                            length(ext2GeneID(rownames(object$cv.fit$fit$beta)[vars.tmp[ii]]))
                    }
                    nvars.vector
                }else{
                    length(unique(ext2GeneID(rownames(object$cv.fit$fit$beta)[vars.tmp])))
                }
            }

        }else{
            if(is.null(x)){
                predict(object = object$cv.fit$fit, type = type,
                        lambda = lambda, ...)
            }else{
                # extended genes
                commonFeat.ext <- unlist(object$complexes.dt)

                # New names of extended genes
                # The new name consists of "group+_+gene name"
                commonFeat.extName <- c()
                for(i in 1:length(object$complexes.dt)){
                    names.i <- paste0(names(object$complexes.dt)[i], "_",
                                      object$complexes.dt[[i]])
                    commonFeat.extName <- c(commonFeat.extName, names.i)
                }

                # extended dataset
                x.ext <- x[, commonFeat.ext]
                colnames(x.ext) <- commonFeat.extName

                predict(object = object$cv.fit$fit, X = x.ext,
                        type = type, lambda = lambda, ...)
            }
        }
    }


#' Make predictions from a PCLasso2 model
#'
#' @description Similar to other predict methods, this function returns
#' predictions from a fitted \code{PCLasso2} object.
#'
#' @param object Fitted \code{PCLasso2} model object.
#' @param x Matrix of values at which predictions are to be made. The features
#' (genes) contained in \code{x} should be consistent with those contained in
#' \code{x} in the \code{PCLasso2} function.  Not used for type="coefficients"
#' or for some of the type settings in \code{predict}.
#' @param type Type of prediction: "link" returns the linear predictors;
#'   "response" gives the risk (i.e., exp(link)); "class" returns the binomial
#'   outcome with the highest probability; "vars" returns the indices for the
#'   nonzero coefficients; "vars.unique" returns unique features (genes)
#'   with nonzero coefficients (If a feature belongs to multiple groups and
#'   multiple groups are selected, the feature will be repeatedly selected.
#'   Compared with "var", "var.unique" will filter out repeated features.);
#'   "groups" returns the groups with at least one nonzero coefficient; "nvars"
#'   returns the number of nonzero coefficients; "nvars.unique" returns the
#'   number of unique features (genes) with nonzero coefficients; "ngroups"
#'   returns the number of groups with at least one nonzero coefficient; "norm"
#'   returns the L2 norm of the coefficients in each group.
#' @param lambda Values of the regularization parameter \code{lambda} at which
#'   predictions are requested. For values of \code{lambda} not in the sequence
#'   of fitted models, linear interpolation is used.
#' @param ... Arguments to be passed to \code{predict.grpreg} in the R package
#' \code{grpreg}.
#' @details
#' See \code{predict.grpreg} in the R package \code{grpreg} for details.
#' @return The object returned depends on \code{type}.
#' @seealso \code{\link{PCLasso2}}
#' @importFrom stats predict
#' @export
#'
#' @examples
#' # load data
#' data(classData)
#' data(PCGroups)
#'
#' x <- classData$Exp
#' y <- classData$Label
#' PC.Human <- getPCGroups(Groups = PCGroups, Organism = "Human",
#' Type = "GeneSymbol")
#'
#' set.seed(20150122)
#' idx.train <- sample(nrow(x), round(nrow(x)*2/3))
#' x.train <- x[idx.train,]
#' y.train <- y[idx.train]
#' x.test <- x[-idx.train,]
#' y.test <- y[-idx.train]
#'
#' # fit PCLasso2 model
#' fit.PCLasso2 <- PCLasso2(x = x.train, y = y.train, group = PC.Human,
#'                   penalty = "grLasso", family = "binomial")
#'
#' # predict risk scores of samples in x.test
#' s <- predict(object = fit.PCLasso2, x = x.test, type="link",
#' lambda=fit.PCLasso2$fit$lambda)
#'
#' s <- predict(object = fit.PCLasso2, x = x.test, type="class",
#' lambda=fit.PCLasso2$fit$lambda[10])
#'
#' # Nonzero coefficients
#' sel.groups <- predict(object = fit.PCLasso2, type="groups",
#'                       lambda = fit.PCLasso2$fit$lambda)
#' sel.ngroups <- predict(object = fit.PCLasso2, type="ngroups",
#'                        lambda = fit.PCLasso2$fit$lambda)
#' sel.vars.unique <- predict(object = fit.PCLasso2, type="vars.unique",
#'                           lambda = fit.PCLasso2$fit$lambda)
#' sel.nvars.unique <- predict(object = fit.PCLasso2, type="nvars.unique",
#'                             lambda = fit.PCLasso2$fit$lambda)
#' sel.vars <- predict(object = fit.PCLasso2, type="vars",
#'                     lambda=fit.PCLasso2$fit$lambda)
#' sel.nvars <- predict(object = fit.PCLasso2, type="nvars",
#'                      lambda=fit.PCLasso2$fit$lambda)
#'
#' # For values of lambda not in the sequence of fitted models,
#' # linear interpolation is used.
#' sel.groups <- predict(object = fit.PCLasso2, type="groups",
#'                       lambda = c(0.1, 0.05))
#' sel.ngroups <- predict(object = fit.PCLasso2, type="ngroups",
#'                        lambda = c(0.1, 0.05))
#' sel.vars.unique <- predict(object = fit.PCLasso2, type="vars.unique",
#'                            lambda = c(0.1, 0.05))
#' sel.nvars.unique <- predict(object = fit.PCLasso2, type="nvars.unique",
#'                             lambda = c(0.1, 0.05))
#' sel.vars <- predict(object = fit.PCLasso2, type="vars",
#'                     lambda=c(0.1, 0.05))
#' sel.nvars <- predict(object = fit.PCLasso2, type="nvars",
#'                      lambda=c(0.1, 0.05))
predict.PCLasso2 <-
    function(object, x = NULL,
             type = c("link", "response", "class", "norm", "coefficients",
                      "vars", "nvars","vars.unique", "nvars.unique", "groups",
                      "ngroups"),
             lambda, ...){

        type <- match.arg(type)

        if(type == "vars.unique"){
            vars.tmp <- predict(object = object$fit,
                                type = "vars", lambda = lambda, ...)

            if(is.list(vars.tmp)){
                vars.list <- vector(mode = "list", length = length(vars.tmp))
                names(vars.list) <- names(vars.tmp)
                for(vars.list.i in 1:length(vars.tmp)){
                    if(length(vars.tmp[[vars.list.i]]) > 0){
                        vars.list[[vars.list.i]] <-
                            unique(ext2GeneID(rownames(object$fit$beta)[vars.tmp[[vars.list.i]]]))
                    }else{
                        vars.list[[vars.list.i]] <- vars.tmp[[vars.list.i]]
                    }

                }
                vars.list
            }else{if(length(lambda) > 1){
                vars.vector <- rep(NA, length = length(vars.tmp))
                names(vars.vector) <- names(vars.tmp)
                for(ii in 1:length(vars.tmp)){
                    vars.vector[ii] <-
                        ext2GeneID(rownames(object$fit$beta)[vars.tmp[ii]])
                }
                vars.vector
            }else{
                unique(ext2GeneID(rownames(object$fit$beta)[vars.tmp]))
            }
            }

        }else if(type == "nvars.unique"){
            vars.tmp <- predict(object = object$fit,
                                type = "vars", lambda = lambda, ...)

            if(is.list(vars.tmp)){
                vars.list <- vector(mode = "list", length = length(vars.tmp))
                names(vars.list) <- names(vars.tmp)
                nvars.vector <- rep(0, length = length(vars.tmp))
                names(nvars.vector) <- names(vars.tmp)
                for(vars.list.i in 1:length(vars.tmp)){
                    if(length(vars.tmp[[vars.list.i]]) > 0){
                        vars.list[[vars.list.i]] <-
                            unique(ext2GeneID(rownames(object$fit$beta)[vars.tmp[[vars.list.i]]]))
                        nvars.vector[vars.list.i] <-
                            length(vars.list[[vars.list.i]])
                    }
                }
                nvars.vector
            }else{
                if(length(lambda) > 1){
                    nvars.vector <- rep(NA, length = length(vars.tmp))
                    names(nvars.vector) <- names(vars.tmp)
                    for(ii in 1:length(vars.tmp)){
                        nvars.vector[ii] <-
                            length(ext2GeneID(rownames(object$fit$beta)[vars.tmp[ii]]))
                    }
                    nvars.vector
                }else{
                    length(unique(ext2GeneID(rownames(object$fit$beta)[vars.tmp])))
                }

            }


        }else{
            if(is.null(x)){
                predict(object = object$fit, type = type,
                        lambda = lambda, ...)
            }else{
                # extended genes
                commonFeat.ext <- unlist(object$complexes.dt)

                # New names of extended genes
                # The new name consists of "complexes+_+gene name"
                commonFeat.extName <- c()
                for(i in 1:length(object$complexes.dt)){
                    names.i <- paste0(names(object$complexes.dt)[i], "_",
                                      object$complexes.dt[[i]])
                    commonFeat.extName <- c(commonFeat.extName, names.i)
                }

                # extended dataset
                x.ext <- x[, commonFeat.ext]
                colnames(x.ext) <- commonFeat.extName

                predict(object = object$fit, X = x.ext,
                        type = type, lambda = lambda, ...)
            }
        }
    }


#' Make predictions from a cross-validated PCLasso2 model
#'
#' @description
#' Similar to other predict methods, this function returns predictions from a
#' fitted "cv.PCLasso2" object, using the optimal value chosen for \code{lambda}.
#'
#' @param object Fitted \code{cv.PCLasso2} model object.
#' @param x
#' Matrix of values at which predictions are to be made. The features (genes)
#' contained in \code{x} should be consistent with those contained in \code{x}
#' in the \code{cv.PCLasso2} function.  Not used for type="coefficients" or for
#' some of the type settings in \code{predict}.
#' @param type Type of prediction: "link" returns the linear predictors;
#'   "response" gives the risk (i.e., exp(link)); "class" returns the binomial
#'   outcome with the highest probability; "vars" returns the indices for the
#'   nonzero coefficients; "vars.unique" returns unique features (genes)
#'   with nonzero coefficients (If a feature belongs to multiple groups and
#'   multiple groups are selected, the feature will be repeatedly selected.
#'   Compared with "var", "var.unique" will filter out repeated features.);
#'   "groups" returns the groups with at least one nonzero coefficient; "nvars"
#'   returns the number of nonzero coefficients; "nvars.unique" returns the
#'   number of unique features (genes) with nonzero coefficients; "ngroups"
#'   returns the number of groups with at least one nonzero coefficient; "norm"
#'   returns the L2 norm of the coefficients in each group.
#' @param lambda
#' Values of the regularization parameter \code{lambda} at which predictions are
#' requested. For values of  \code{lambda} not in the sequence of fitted models,
#' linear interpolation is used.
#' @param ...
#' Arguments to be passed to \code{predict.cv.grpreg} in the R package
#' \code{grpreg}.
#'
#' @return
#' The object returned depends on \code{type}.
#' @method predict cv.PCLasso2
#' @importFrom stats predict
#' @export
#'
#' @seealso \code{\link{cv.PCLasso2}}
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
#' #' set.seed(20150122)
#' idx.train <- sample(nrow(x), round(nrow(x)*2/3))
#' x.train <- x[idx.train,]
#' y.train <- y[idx.train]
#' x.test <- x[-idx.train,]
#' y.test <- y[-idx.train]
#'
#' # fit model
#' cv.fit1 <- cv.PCLasso2(x = x.train, y = y.train, group = PC.Human,
#' penalty = "grLasso", family = "binomial", nfolds = 10)
#'
#' # predict risk scores of samples in x.test
#' s <- predict(object = cv.fit1, x = x.test, type="class",
#'              lambda=cv.fit1$cv.fit$lambda.min)
#'
#' # Nonzero coefficients
#' sel.groups <- predict(object = cv.fit1, type="groups",
#'                       lambda = cv.fit1$cv.fit$lambda.min)
#' sel.ngroups <- predict(object = cv.fit1, type="ngroups",
#'                        lambda = cv.fit1$cv.fit$lambda.min)
#' sel.vars.unique <- predict(object = cv.fit1, type="vars.unique",
#'                            lambda = cv.fit1$cv.fit$lambda.min)
#' sel.nvars.unique <- predict(object = cv.fit1, type="nvars.unique",
#'                             lambda = cv.fit1$cv.fit$lambda.min)
#' sel.vars <- predict(object = cv.fit1, type="vars",
#'                     lambda=cv.fit1$cv.fit$lambda.min)
#' sel.nvars <- predict(object = cv.fit1, type="nvars",
#'                      lambda=cv.fit1$cv.fit$lambda.min)
predict.cv.PCLasso2 <-
    function(object, x = NULL,
             type = c("link", "response", "class", "norm", "coefficients",
                      "vars", "nvars","vars.unique", "nvars.unique", "groups",
                      "ngroups"),
             lambda, ...){

        type <- match.arg(type)
        if(type == "vars.unique"){

            vars.tmp <- predict(object = object$cv.fit$fit,
                                type = "vars", lambda = lambda, ...)

            if(is.list(vars.tmp)){
                vars.list <- vector(mode = "list", length = length(vars.tmp))
                names(vars.list) <- names(vars.tmp)
                for(vars.list.i in 1:length(vars.tmp)){
                    if(length(vars.tmp[[vars.list.i]]) > 0){
                        vars.list[[vars.list.i]] <-
                            unique(ext2GeneID(rownames(object$cv.fit$fit$beta)[vars.tmp[[vars.list.i]]]))
                    }else{
                        vars.list[[vars.list.i]] <- vars.tmp[[vars.list.i]]
                    }

                }
                vars.list
            }else{
                if(length(lambda) > 1){
                    vars.vector <- rep(NA, length = length(vars.tmp))
                    names(vars.vector) <- names(vars.tmp)
                    for(ii in 1:length(vars.tmp)){
                        vars.vector[ii] <-
                            ext2GeneID(rownames(object$cv.fit$fit$beta)[vars.tmp[ii]])
                    }
                    vars.vector
                }else{
                    unique(ext2GeneID(rownames(object$cv.fit$fit$beta)[vars.tmp]))
                }
            }
        }else if(type == "nvars.unique"){
            vars.tmp <- predict(object = object$cv.fit$fit,
                                type = "vars", lambda = lambda, ...)

            if(is.list(vars.tmp)){
                vars.list <- vector(mode = "list", length = length(vars.tmp))
                names(vars.list) <- names(vars.tmp)
                nvars.vector <- rep(0, length = length(vars.tmp))
                names(nvars.vector) <- names(vars.tmp)
                for(vars.list.i in 1:length(vars.tmp)){
                    if(length(vars.tmp[[vars.list.i]]) > 0){
                        vars.list[[vars.list.i]] <-
                            unique(ext2GeneID(rownames(object$cv.fit$fit$beta)[vars.tmp[[vars.list.i]]]))
                        nvars.vector[vars.list.i] <-
                            length(vars.list[[vars.list.i]])
                    }
                }
                nvars.vector
            }else{
                if(length(lambda) > 1){
                    nvars.vector <- rep(NA, length = length(vars.tmp))
                    names(nvars.vector) <- names(vars.tmp)
                    for(ii in 1:length(vars.tmp)){
                        nvars.vector[ii] <-
                            length(ext2GeneID(rownames(object$cv.fit$fit$beta)[vars.tmp[ii]]))
                    }
                    nvars.vector
                }else{
                    length(unique(ext2GeneID(rownames(object$cv.fit$fit$beta)[vars.tmp])))
                }
            }

        }else{
            if(is.null(x)){
                predict(object = object$cv.fit$fit, type = type,
                        lambda = lambda, ...)
            }else{
                # extended genes
                commonFeat.ext <- unlist(object$complexes.dt)

                # New names of extended genes
                # The new name consists of "group+_+gene name"
                commonFeat.extName <- c()
                for(i in 1:length(object$complexes.dt)){
                    names.i <- paste0(names(object$complexes.dt)[i], "_",
                                      object$complexes.dt[[i]])
                    commonFeat.extName <- c(commonFeat.extName, names.i)
                }

                # extended dataset
                x.ext <- x[, commonFeat.ext]
                colnames(x.ext) <- commonFeat.extName

                predict(object = object$cv.fit$fit, X = x.ext,
                        type = type, lambda = lambda, ...)
            }
        }
    }


