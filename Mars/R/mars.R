#' Multivariate Adaptive Regression Splines (MARS)
#'TEST
#' Fit Friedman's Multivariate Adaptive Regression Splines (MARS) model.
#'
#' @param formula an R formula
#' @param data a data frame containing the data
# ....
# .....

mars <- function(formula,data,control=mars.control()) {
  cc <- match.call() # save the call
  mf <- model.frame(formula,data)
  y <- model.response(mf)
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf)[,-1,drop=FALSE]
  x_names <- colnames(x)
  control <- validate_mars.control(control)
  fwd <- fwd_stepwise(y,x,control)
  bwd <- bwd_stepwise(fwd,control)
  fit <- lm(y~.-1,data=data.frame(y=y,bwd$B)) # notice -1 added
  out <- c(list(call=cc,formula=formula,y=y,B=bwd$B,Bfuncs=bwd$Bfuncs,
                x_names=x_names),fit)
  class(out) <- c("mars",class(fit))
  out
}


fwd_stepwise <- function(y,x,control=mars.control()){
  Mmax = control$Mmax

  # Initialize:
  N <- length(y) # sample size
  n <- ncol(x) # number of predictors = number of X
  B <- init_B(N,Mmax) # Exercise: write init_B()
  Bfuncs <- vector("list", length = Mmax+1)

  # Looping for forward selection:
  for(i in 1:(Mmax/2)) { # contrast to indexing 2...Mmax in Friedman
    lof_best <- Inf
    M <- 2*i-1
    for(m in 1:M) { # choose a basis function to split
      for (v in setdiff(1:n, Bfuncs[[m]][, "v"])) {
        tt <- split_points(x[,v],B[,m]) # Exercise: write split_points()
        for(t in tt) {
          Bnew <- data.frame(B[,(1:M)], # drop m-th col: B[,-m]
                             # replace parent B[,m] with Btem1,Btem2
                             Btem1=B[,m]*h(x[,v],+1, t),
                             Btem2=B[,m]*h(x[,v],-1, t))
          gdat <- data.frame(y=y,Bnew)
          lof <- LOF(y~.,gdat,mc) #  Use your LOF() from week 4
          if(lof < lof_best) {
            lof_best <- lof
            best_split <- c(m=m,v=v,t=t)} } } }
    mstar <- best_split["m"]; vstar <- best_split["v"]; tstar <- best_split["t"]
    cat("[Info] best (m,v,t,lof): (",mstar,vstar,tstar,lof_best,")\n")
    B[,M+1] <- B[,mstar]*h(x[,vstar],-1,tstar)
    B[,M+2] <- B[,mstar]*h(x[,vstar],+1,tstar)
    Bfuncs[[M+1]] = rbind(Bfuncs[[mstar]], c(s = -1, vstar, tstar))
    Bfuncs[[M+2]] = rbind(Bfuncs[[mstar]], c(s = +1, vstar, tstar)) # Update parent basis with the right cihld basis
  } # end loop over M
  colnames(B) <- paste0("B",(0:(ncol(B)-1)))
  return(list(y=y,B=B, Bfuncs=Bfuncs))
}


init_B <- function(N,Mmax) {
  B <- data.frame(matrix(NA,nrow=N,ncol=(Mmax+1)))
  B[,1] <- 1
  names(B) <- c("B0",paste0("B",1:Mmax))
  return(B)
}

bwd_stepwise <- function(fwd,control) {
  y <- fwd$y
  B <- fwd$B
  Mmax <- ncol(fwd$B) - 1
  Jstar <- 2:(Mmax+1)
  Kstar <- Jstar
  dat <- data.frame(y = y, B[,-1])
  LOFstar <- LOF(y~.-1, dat, control)
  for(M in (Mmax+1):2) {
    b <- Inf
    L <- Kstar
    for(m in L) {
      K <- setdiff(L,m)
      dat <- data.frame(y = fwd$y, fwd$B[,K])
      lof <- LOF(y~.-1, dat, control)
      if(lof < b) {
        b <- lof
        Kstar <- K
      }
      if(lof < LOFstar) {
        LOFstar <- lof
        Jstar <- K
      }
    }
  }
  Jstar <- c(1,Jstar)
  return(list(y = fwd$y, B = fwd$B[,Jstar], Bfuncs = fwd$Bfuncs[Jstar]))
}

LOF <- function(form,data,control) {
  fit <- lm(formula, data)
  rss <- sum(resid(fit)^2)

  # Calculate number of rows and columns of basis matrix
  N <- nrow(data)
  M <- length(coefficients(fit)) - 1
  Ctilde <- sum(hatvalues(fit))+(mars_control$d)*M

  # Calculate GCV criterion
  GCV <- rss*N / (N - Ctilde)^2

  return(GCV)
}

h <- function(x,s,t) {
  return(pmax(0,s*(x-t)))
}

split_points <- function(xv,Bm) {
  out <- sort(unique(xv[Bm>0]))
  return(out[-length(out)])
}


#------------------------------------------------------------------------
# constructor, validator and helper for class mars.control
#------------------------------------------------------------------------
#
new_mars.control <- function(control) {
  structure(control,class="mars.control")
}

validate_mars.control <- function(control) {
  stopifnot(is.integer(control$Mmax),is.numeric(control$d),
            is.logical(control$trace))
  if(control$Mmax < 2){
    warning("Mmax must be >= 2; Reset it to 2")
    control$Mmax <- 2}
  if(control$Mmax %% 2 > 0){
    control$Mmax <- 2*ceiling(control$Mmax/2)
    warning("Mmax should be an even integer. Reset it to ", control$Mmax) }
  control
}


#' Constructor for `mars.control` objects
#'
#' This function constructs a `mars.control` object that specifies
#' parameters used in the model fitting procedure.
#'
#' @param Mmax Maximum number of basis functions. Should be an even integer. Default value is 2.
# .....
# ...

mars.control <- function(Mmax=2,d=3,trace=FALSE) {
  Mmax <- as.integer(Mmax)
  control <- list(Mmax=Mmax,d=d,trace=trace)
  control <- validate_mars.control(control)
  new_mars.control(control)
}



