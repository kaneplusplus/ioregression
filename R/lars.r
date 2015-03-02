#' Perform a lasso or stepwise regression
#'
#' This function takes a pre-computed iolm object,
#' and applies the Least Angle Regression (LARs)
#' algorithm to it. By using the iolm object, this
#' function can run without needing to make an additional
#' pass over the data.
#'
#' @importFrom       lars updateR nnls.lars delcol
#' @param object     an iolm object with the desired formula
#' @param type       One of "lasso", "lar", "forward.stagewise" or "stepwise". The
#'                   names can be abbreviated to any unique substring. Default is
#'                   "lasso".
#' @param normalize  If TRUE, each variable is standardized to have unit L2 norm,
#'                   otherwise it is left alone. Default is TRUE.
#' @param intercept  if TRUE, an intercept is included in the model (and not
#'                   penalized), otherwise no intercept is included. If missing,
#'                   will be infered from the formular in \code{object}.
#' @param eps        An effective zero
#' @param max.steps  Limit the number of steps taken
#' @param out.type   The desired class of the output. Either "iolars" or "lars".
#'                   The latter requires installing the lars package in order to
#'                   print, predict, ect.
#' @export
iolars = function(object, type = c("lasso", "lar", "forward.stagewise","stepwise"),
                  normalize=TRUE, intercept, eps = .Machine$double.eps, max.steps = NULL,
                  out.type = c("iolars", "lars")) {
  if (!inherits(object, "iolm"))
    stop("input to object must be an iolm object! See ?ioregression::iolm for more info.")

  out.type = match.arg(out.type)
  if (attr(object$terms, "intercept")) {
    XtX = as.matrix(object$xtx)[-1,-1]
    XtY = as.matrix(object$xty)[-1]
    mean_x = as.numeric(object$mean_x)[-1]
    noms = names(object$coefficients)[-1]
    if (missing(intercept)) intercept = TRUE
  } else {
    XtX = as.matrix(object$xtx)
    XtY = as.matrix(object$xty)
    mean_x = as.numeric(object$mean_x)
    noms = names(object$coefficients)
    if (missing(intercept)) intercept = FALSE
  }

  obj = lars_par(XtX, XtY, as.numeric(object$yty),
                  n = object$n, mean_x=mean_x, mean_y=object$sum_y / object$n,
                  type = type, normalize=normalize,
                  intercept=intercept, eps = eps, max.steps = max.steps)
  obj$call = match.call()

  if (out.type == "lars") return(obj)

  obj$info = cbind(lambda=c(obj$lambda,0), action=c(0L, obj$action),
                    df=obj$df, RSS=obj$RSS, R2=obj$R2, obj=obj$RSS,
                    Cp=obj$Cp)
  rownames(obj$info) = rep("",nrow(obj$info))
  attr(obj$beta,"scaled:scale") = NULL

  if (intercept) {
    obj$beta = cbind(obj$mu - obj$beta %*% obj$meanx, obj$beta)
    colnames(obj$beta) = c("(Intercept)", noms)
  } else {
    names(obj$beta) = noms
  }

  class(obj) = "iolars"
  obj
}

#' Get coefficients of an iolars object
#'
#' @param object  an iolars object to get the coefficients matrix for
#' @param ...     other inputs; currently unused
#' @export
coef.iolars = function (object, ...) object$beta

#' Print an iolars object
#'
#' @method print iolars
#' @param x   an iolars object to print the summary of
#' @param ... other inputs; currently unused
#' @export
print.iolars = function (x, ...) print(x$info)

#' This is a modified version of the lars function from the lars package,
#' which allows for pre-computing XtY.
#'
#' @param XtX        precomputed XtX matrix
#' @param XtY        precomputed XtY vector
#' @param YtY        precomputed YtY value
#' @param n          number of data points in the original X matrix
#' @param mean_x     vector of means for each column of the X matrix
#' @param mean_y     mean of the response vector Y
#' @param type       One of "lasso", "lar", "forward.stagewise" or "stepwise". The
#'                   names can be abbreviated to any unique substring. Default is
#'                   "lasso".
#' @param normalize  If TRUE, each variable is standardized to have unit L2 norm,
#'                   otherwise it is left alone. Default is TRUE.
#' @param intercept  if TRUE, an intercept is included in the model (and not
#'                   penalized), otherwise no intercept is included. Default is
#'                   TRUE.
#' @param eps        An effective zero
#' @param max.steps Limit the number of steps taken
lars_par <-
function(XtX, XtY, YtY, n, mean_x, mean_y,
          type = c("lasso", "lar", "forward.stagewise","stepwise"),
          normalize=TRUE, intercept=FALSE, eps = .Machine$double.eps,
          max.steps = NULL)
{
### Lazy hack for now, rather than fiddling with the rest of the code
  Gram = XtX
  Cvec = XtY
  trace = FALSE
  x = matrix(0)
  y = 0L
  use.Gram = TRUE

### program automatically centers and standardizes predictors by default.
###
### Original program by Brad Efron September 2001
### Recoded by Trevor Hastie November 2001
### Computational efficiency December 22, 2001
### Bug fixes and singularities February 2003
### Conversion to R April 2003
### stepwise and non-standardize options added May 2007
### Copyright Brad Efron and Trevor Hastie
###
### Edited to use precomputed XtX and XtY by Taylor Arnold 2015
  call <- match.call()
  type <- match.arg(type)
  TYPE <- switch(type,
                 lasso = "LASSO",
                 lar = "LAR",
                 forward.stagewise = "Forward Stagewise",
                 stepwise = "Forward Stepwise")
  if(trace)
    cat(paste(TYPE, "sequence\n"))

  # Taylor Edit:
  # OLD: nm <- dim(x)
  nm = c(n,ncol(Gram))
  n <- nm[1]
  m <- nm[2]
  im <- inactive <- seq(m)
  one <- rep(1, n)
  vn <- dimnames(x)[[2]]
### Center x and y, and scale x, and save the means and sds
  if(intercept){
    mu <- mean_y
    Cvec = Cvec - n*mean_x*mean_y
    Gram = Gram - outer(mean_x,mean_x)*n
  } else {
    mean_x = rep(0,m)
    mu <- 0
  }
  if(normalize){
    normx <- sqrt(diag(Gram))
    names(normx) <- NULL
    Gram = Gram / outer(normx,normx)
    Cvec = Cvec / normx
  } else {
    normx <- rep(1,m)
  }
  ignores <- NULL
  # Taylor Edits:
  # Cvec <- drop(t(y) %*% x)
  # ssy <- sum(y^2) ### Some initializations
  # residuals <- y
  ssy = 0
  residuals = 0

  if(is.null(max.steps))
    max.steps <- 8*min(m, n-intercept)
  beta <- matrix(0, max.steps + 1, m) # beta starts at 0
  lambda=double(max.steps)
  Gamrat <- NULL
  arc.length <- NULL
  R2 <- 1
  RSS <- ssy
  first.in <- integer(m)
  active <- NULL  # maintains active set
  actions <- as.list(seq(max.steps))
                                        # a signed index list to show what comes in and out
  drops <- FALSE  # to do with type=="lasso" or "forward.stagewise"
  Sign <- NULL  # Keeps the sign of the terms in the model
  R <- NULL ###
### Now the main loop over moves
###
  k <- 0
  while((k < max.steps) & (length(active) < min(m - length(ignores),n-intercept)) )
    {
      action <- NULL
      C <- Cvec[inactive] #
### identify the largest nonactive gradient
      Cmax <- max(abs(C))
      if(Cmax<eps*100){ # the 100 is there as a safety net
        if(trace)cat("Max |corr| = 0; exiting...\n")
        break
      }
      k <- k + 1
      lambda[k]=Cmax
### Check if we are in a DROP situation
      if(!any(drops)) {
        new <- abs(C) >= Cmax - eps
        C <- C[!new]  # for later
        new <- inactive[new]  # Get index numbers
### We keep the choleski R  of X[,active] (in the order they enter)
        for(inew in new) {
          if(use.Gram) {
            R <- lars::updateR(Gram[inew, inew], R, drop(Gram[
                                                        inew, active]), Gram = TRUE,eps=eps)
          }
          else {
            R <- lars::updateR(x[, inew], R, x[, active], Gram
                         = FALSE,eps=eps)
          }
          if(attr(R, "rank") == length(active)) {
            ##singularity; back out
            nR <- seq(length(active))
            R <- R[nR, nR, drop = FALSE]
            attr(R, "rank") <- length(active)
            ignores <- c(ignores, inew)
            action <- c(action,  - inew)
            if(trace)
              cat("LARS Step", k, ":\t Variable", inew,
                  "\tcollinear; dropped for good\n")  #
          }
          else {
            if(first.in[inew] == 0)
              first.in[inew] <- k
            active <- c(active, inew)
            Sign <- c(Sign, sign(Cvec[inew]))
            action <- c(action, inew)
            if(trace)
              cat("LARS Step", k, ":\t Variable", inew,
                  "\tadded\n")  #
          }
        }
      }
      else action <-  - dropid
      Gi1 <- backsolve(R, backsolve(R, Sign, ncol(R), transpose = TRUE))
### Now we have to do the forward.stagewise dance
### This is equivalent to NNLS
      dropouts<-NULL
      if(type == "forward.stagewise") {
        directions <- Gi1 * Sign
        if(!all(directions > 0)) {
          if(use.Gram) {
            nnls.object <- lars::nnls.lars(active, Sign, R,
                                     directions, Gram[active, active], trace =
                                     trace, use.Gram = TRUE,eps=eps)
          }
          else {
            nnls.object <- lars::nnls.lars(active, Sign, R,
                                     directions, x[, active], trace = trace,
                                     use.Gram = FALSE,eps=eps)
          }
          positive <- nnls.object$positive
          dropouts <-active[-positive]
          action <- c(action, -dropouts)
          active <- nnls.object$active
          Sign <- Sign[positive]
          Gi1 <- nnls.object$beta[positive] * Sign
          R <- nnls.object$R
          C <- Cvec[ - c(active, ignores)]
        }
      }
      A <- 1/sqrt(sum(Gi1 * Sign))
      w <- A * Gi1  # note that w has the right signs
      if(!use.Gram) u <- drop(x[, active, drop = FALSE] %*% w)  ###
### Now we see how far we go along this direction before the
### next competitor arrives. There are several cases
###
### If the active set is all of x, go all the way
      if( (length(active) >=  min(n-intercept, m - length(ignores) ) )|type=="stepwise") {
        gamhat <- Cmax/A
      }
      else {
        if(use.Gram) {
          a <- drop(w %*% Gram[active,  - c(active,ignores), drop = FALSE])
        }
        else {
          a <- drop(u %*% x[,  - c(active, ignores), drop=FALSE])
        }
        gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))
### Any dropouts will have gam=0, which are ignored here
        gamhat <- min(gam[gam > eps], Cmax/A)
      }
      if(type == "lasso") {
        dropid <- NULL
        b1 <- beta[k, active] # beta starts at 0
        z1 <-  - b1/w
        zmin <- min(z1[z1 > eps], gamhat)
        if(zmin < gamhat) {
          gamhat <- zmin
          drops <- z1 == zmin
        }
        else drops <- FALSE
      }
      beta[k + 1,  ] <- beta[k,  ]
      beta[k + 1, active] <- beta[k + 1, active] + gamhat * w
      if(use.Gram) {
        Cvec <- Cvec - gamhat * Gram[, active, drop = FALSE] %*% w
      }
      else {
        residuals <- residuals - gamhat * u
        Cvec <- drop(t(residuals) %*% x)
      }
      Gamrat <- c(Gamrat, gamhat/(Cmax/A))
      arc.length <- c(arc.length, gamhat)
### Check if we have to drop any guys
      if(type == "lasso" && any(drops)) {
        dropid <- seq(drops)[drops]
                                        #turns the TRUE, FALSE vector into numbers
        for(id in rev(dropid)) {
          if(trace)
            cat("Lasso Step", k+1, ":\t Variable", active[
                                                        id], "\tdropped\n")
          if ((p <- dim(R)[1]) == 1)
           R = NULL
          else
            R <- lars::delcol(R, rep(1, p), id)[[1]][-p, , drop = FALSE]
          attr(R, "rank") = p - 1
        }
        dropid <- active[drops] # indices from 1:m
        beta[k+1,dropid]<-0  # added to make sure dropped coef is zero
        active <- active[!drops]
        Sign <- Sign[!drops]
      }
      if(!is.null(vn))
        names(action) <- vn[abs(action)]
      actions[[k]] <- action
      inactive <- im[ - c(active, ignores)]
      if(type=="stepwise")Sign=Sign*0
    }
  beta <- beta[seq(k + 1), ,drop=FALSE ]  #
  lambda=lambda[seq(k)]
  dimnames(beta) <- list(paste(0:k), vn)  ### Now compute RSS and R2
  if(trace)
    cat("Computing residuals, RSS etc .....\n")

  # Taylor Edit:
  beta <- scale(beta, FALSE, normx)
  MtM = XtX - outer(mean_x,mean_x)*n
  RSS = drop(YtY + diag(beta %*% MtM %*% t(beta)) - 2 * beta %*% XtY
             - mu^2*n + 2 * n * mu * beta %*% mean_x)
  R2 <- drop(1 - RSS/RSS[1])
  actions=actions[seq(k)]
  netdf=sapply(actions,function(x)sum(sign(x)))
  df=cumsum(netdf)### This takes into account drops
  if(intercept) df=c(Intercept=1,df+1)
    else df=c(Null=0,df)
  rss.big=rev(RSS)[1]
  df.big=n-rev(df)[1]
  if(rss.big<eps|df.big<eps) sigma2=NaN
   else sigma2=rss.big/df.big
  Cp <- drop(RSS/sigma2 - n + 2 * df)
  attr(Cp,"sigma2")=sigma2
  attr(Cp,"n")=n
  attributes(Gram)$dimnames <- NULL
  object <- list(call = call, type = TYPE, df=df, lambda=lambda, R2 = R2, RSS = RSS, Cp = Cp,
                 actions = actions[seq(k)], entry = first.in, Gamrat = Gamrat,
                 arc.length = arc.length, Gram = if(use.Gram) Gram else NULL,
                 beta = beta, mu = mu, normx = normx, meanx = mean_x)
  class(object) <- "lars"
  object
}
