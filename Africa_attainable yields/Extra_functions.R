#- 1. Percent NA one input-==================================================

Nalength <-function(x){
  
  length(which(is.na(x)))/length(x)
}

#- 2. Percent NA multiple input==============================================

MultiNalength <- function(data, mar=c(13,5,4,2)){
  
  out <- lapply(data, function(x){
    
    length(which(is.na(x)))/length(x)
    
  })
  
  names(out) <- names(data)
  
  # Plot--------------------
  par(mar=mar)
  k <- data.frame(x = names(out), y=unname(unlist(out)))
  # Calculate means
  means <- tapply(k$y, k$x, mean, na.rm = TRUE)
  # Barplot
  par(mar = mar)
  bar_colors <- ifelse(means < 1, "green", "red")
  barplot(means, ylab = "Proportion, NA", xlab = "", col = bar_colors,las = 2)
  
}

# 3. Data availability-======================================================

AvailableData <- function(data, mar=c(13,5,4,2)){
  
  out <- lapply(data, function(x){
    
    1-length(which(is.na(x)))/length(x)
    
  })
  
  names(out) <- names(data)
  
  # Plot--------------------
  par(mar = mar, bty = "7")
  k <- data.frame(x = names(out), y=unname(unlist(out)))
  # Calculate means
  means <- tapply(k$y, k$x, mean, na.rm = TRUE)
  # Barplot
  vals <- pmin(pmax(means, 0), 1)          # clamp to [0,1]
  bar_cols <- hsv(h = (1/3)*vals, s = 0.95, v = 0.95)  # 0=red … 1/3=green
  bar_cols[is.na(means)] <- "grey80"
  barplot(means, ylab = "Proportion, Available data", xlab = "", 
          col = bar_cols,las = 2, ylim=c(0,max(means)+0.1))
  
}

# 4. Data proportion availability ==============================================

Props <- function(bas, check, data, mar=c(13,5,4,2)){
  
  out1 <- vector()
  
  for(i in 1: length(check)){
    
    out1[i] <- length(which(is.na(data[,bas])==F & is.na(data[,check[i]])==F))/length(which(is.na(data[,bas])==F))
  }
  
  y=c(1, out1)
  
  out2<- data.frame(Factor=c(bas,check), Response=y, n=y*length(which(is.na(data[,bas])==F)))
  
  
  # Plot --------------
  par(mar=mar)
  means <- tapply(out2$Response, out2$Factor, mean, na.rm = TRUE)
  
  vals <- pmin(pmax(means, 0), 1)          # clamp to [0,1]
  bar_cols <- hsv(h = (1/3)*vals, s = 0.95, v = 0.95)  # 0=red … 1/3=green
  bar_cols[is.na(means)] <- "grey80"
  barplot(means, ylab = "Proportion, Available data", xlab = "", 
          col = bar_cols,las = 2, ylim=c(0,max(means)))
  
  return(out2) 
}


# 5. ble_profile extension =====================================================


ble_extention <- function(data, start, sigh, model, plot=T, optim.method="BFGS", equation=NULL){
  
  profile <- lapply(start, function(start) { 
    lapply(sigh, function(sigh) { 
      ble_profile(data = data, start = start, sigh = sigh, model = model, plot = F, 
                  optim.method=optim.method, equation=equation)
    }) 
  }) 
  
  
  log_likelihood <- unlist(lapply(profile, function(sublist) { 
    sapply(sublist, function(x) x[["log-likelihood"]]) 
  })) 
  
  Merror<- unlist(lapply(profile, function(sublist) { 
    sapply(sublist, function(x) x[["Merror"]]) 
  })) 
  
  maximum_likelihood <- Merror[which(log_likelihood==max(log_likelihood[is.finite(log_likelihood)], na.rm = T))]# profile maximized at 6 t/ha
  
  ##-----Plot the max values-----
  if(plot){
    
    me_profile<-data.frame(x=Merror,y=log_likelihood)
    me_profile<-me_profile[order(Merror),] # ordering data from smallest to largest Merror
    
    vec<-vector()
    
    suppressWarnings(
      
      for(i in unique(me_profile$x)){
        maxy<-max(me_profile[which(me_profile$x==i),]$y, na.rm = T)
        vec<-c(vec,maxy)
      }
    )
    
    
    # Ploting the largest value at each Merror
    
    par(mar=c(5,5,4,4))
    
    plot(unique(me_profile$x),vec,pch=16,
         xlab=expression(bold(sigma[e])),
         ylab=expression(bold(log-Likelihood)),
         cex.lab=1.8, cex.axis=1.8)
    lines(unique(me_profile$x),vec, lty=5, lwd=1.5)
    abline(v=Merror[which(log_likelihood==max(log_likelihood, na.rm = T))],
           lty=5, col="red", lwd=1.3) # sigh with largest log-likelihood
    
  }
  
  
  #--------Return result -------
  
  paste("log-likelihood of Measurement error is maximised at", maximum_likelihood)
  
}



# 6. cbvn extension ==============================================================================================


cbvn_extention <- function(data,model, equation=NULL, start, sigh, UpLo="U", optim.method="BFGS",
                           Hessian=FALSE, plot=TRUE, line_smooth=1000, lwd=2, l_col="red",...){
  
  models <- lapply(start, function(start){
    
    tryCatch(
      
      cbvn(data=data, model=model, equation=equation, start=start, sigh=sigh, UpLo=UpLo, optim.method=optim.method,
           Hessian=Hessian, plot=plot, line_smooth=line_smooth, lwd=2, l_col=l_col,...),
      error = function(e) NA)
  })
  
  
  
  model <- models[[which.min(unlist(lapply(X=models, FUN=function(a){
    
    b <- tryCatch(a$AIC[2,1], error = function(e) NA)
    
    return(b)
    
  })))]]
  
  return(model)
  
}

# 7. Akaike weights---------------------------------------------------------------------------------------------------

akweight<-function(aic1,aic2){
  lam_v<-min(aic1,aic2)
  
  del1<-aic1-lam_v
  del2<-aic2-lam_v
  
  wt1<-exp(-del1/2)/(exp(-del1/2)+exp(-del2/2))
  wt2<-exp(-del2/2)/(exp(-del1/2)+exp(-del2/2))
  
  op<-matrix(0,nrow=2,ncol=2)
  rownames(op)<-c("AIC","Akaike weight")
  op[1,]<-c(aic1,aic2)
  op[2,]<-c(wt1,wt2)
  return(op)
}

# 8. Categorical BLA ===========================================================


ble_Cat <- function(x, y, sigh, theta=NULL, optim_method = "L-BFGS-B", 
                    plot = TRUE, upper=0.99) {
  
  # Convert x to factor if not already
  k <- as.factor(x)
  nlev <- length(unique(k))  # Number of levels in the factor
  
  
  # Use different sighs
  
  likelihood <- vector()
  
  for(i in 1:length(sigh)){
    
    # Define the negative log-likelihood for the bounded model
    nll_factor_bound <- function(theta, y, k, nlev, sigh=sigh[i]) {
      mu <- theta[1:nlev]
      bound <- theta[(nlev+1):(2*nlev)]
      sd <- theta[(2*nlev+1)]
      
      n <- length(y)
      llik <- 0
      
      # Calculate max response per level. CM addition
      max_y_per_level <- tapply(y, k, max, na.rm = TRUE)
      
      for (i in 1:n) {
        
        i.lev <- as.numeric(k[i])
        
        # # Constraint: Ensure bound does not exceed max(y) for that level. CM addition
        # if (bound[i.lev] > max_y_per_level[i.lev]) {
        #   return(Inf)  # Penalize invalid solutions
        # }
        
        
        dens <- coffcturb(y[i], mu[i.lev], sd, -Inf, bound[i.lev], sigh)
        llik <- llik + log(dens)
      }
      
      return(-1 * llik)
    }
    
    # Approximate partial derivatives for scaling boundary model for factor
    
    nll_factor_bound_pd<-function(theta,y,k,nlev,sigh=sigh[i]){
      
      eps=1e-4
      nr<-length(theta)
      part<-vector("numeric",nr)
      
      for (i in 1:nr){
        del<-rep(0,nr)
        del[i]<-eps
        part[i]<-(nll_factor_bound((theta+del),y,k,nlev,sigh)-
                    nll_factor_bound(theta,y,k,nlev,sigh))/eps
      }
      
      return(part)
    }
    
    # Define Turban's function for censored normal distribution with measurement error
    coffcturb <- function(x, mu, sig, a, c, sigh=sigh[i]) {
      k <- ((mu - c) / sig)
      d <- ((mu - a) / sig)
      alpha <- ((sigh * sigh) * (x - mu)) / ((sigh * sigh) + (sig * sig))
      beta <- sqrt((sigh * sigh * sig * sig) / ((sigh * sigh) + (sig * sig)))
      gamma <- (beta * sqrt(2 * pi)) / (2 * pi * sig * sigh * (pnorm(d) - pnorm(k)))
      
      com1 <- -((x - mu)^2) / (2 * ((sigh * sigh) + (sig * sig)))
      com2 <- gamma * exp(com1)
      f <- com2 * (pnorm((x - a - alpha) / beta) - pnorm((x - c - alpha) / beta))
      
      # Rescale for censored values
      f <- f * pnorm(c, mu, sig)
      
      # Add contribution at x from mass at c
      f <- f + (dnorm((x - c), 0, sigh) * (1 - pnorm(c, mu, sig)))
      
      return(f)
    }
    
    ## The null-model ----------------------------------------------------------------------
    
    #  Function for negative log likelihood for factor, unbounded
    nll_factor<-function(theta,y,k,nlev){
      
      nlev<-length(theta)-1
      mu<-theta[1:nlev]
      sd<-theta[(nlev+1)]
      n<-length(y)
      llik<-0
      
      for(i in 1:n){
        i.lev<-as.numeric(k[i])
        dens<-dnorm(y[i],mu[i.lev],sd,log=F)
        llik<-llik+log(dens)
      }
      return(-1*llik)
    }
    
    # Function to approximate partial derivatives for scaling boundary model for factor
    nll_factor_pd<-function(theta,y,k,nlev){
      eps=1e-4
      nr<-length(theta)
      part<-vector("numeric",nr)
      
      for (i in 1:nr){
        del<-rep(0,nr)
        del[i]<-eps
        part[i]<-(nll_factor((theta+del),y,k,nlev)-
                    nll_factor(theta,y,k,nlev))/eps
      }
      return(part)
    }
    
    ## Optimization ------------------------------------------------------------------------
    
    # Step 1: Compute scaling for optimization----------------------------------------------
    
    
    ## theta to use-------------------------------------------------------------------
    
    if (is.null(theta)) {
      
      means <- tapply(y,x,"mean")
      upper <- tapply(y, x, function(z) quantile(z, probs=upper))
      lmod<-lm(y~x)
      res <- sd(lmod$residuals)
      
      theta <- unname(c(means,upper,res))
      
    }else{
      
      theta <- theta
    }
    
    ## Calculate approximate partial derivatives at the guess for scaling
    
    likelihood[i] <- tryCatch({
      
      if(optim_method=="L-BFGS-B"){
        
        # Define lower and upper bounds for parameters. CM edits
        lower_bounds <- c(rep(-Inf, nlev), rep(min(y), nlev), 1e-100)  # Boundaries cannot be below min(y)
        upper_bounds <- c(rep(Inf, nlev), tapply(y, k, max, na.rm = TRUE), Inf)  # Ensure upper bound ≤ max(y)
        
        
        pd<-nll_factor_bound_pd(theta,y,k,nlev,sigh = sigh[i])
        scale<-1/abs(pd)
        
        # first estimate
        suppressWarnings(mlest <- optim(theta, nll_factor_bound, method = optim_method,
                                        y = y, k = k, nlev = nlev, sigh = sigh[i],
                                        lower = lower_bounds, upper = upper_bounds))
        
        pd<-nll_factor_bound_pd(mlest$par,y,k,nlev,sigh= sigh[i])
        scale<-1/abs(pd)
        
        # Step 2: Refined estimation
        mlest2 <- optim(mlest$par, nll_factor_bound, method = optim_method,
                        y = y, k = k, nlev = nlev, sigh = sigh[i],,
                        lower = lower_bounds, upper = upper_bounds, hessian = TRUE)
        
      }else{
        pd<-nll_factor_bound_pd(theta,y,k,nlev,sigh= sigh[i])
        scale<-1/abs(pd)
        
        # first estimate
        suppressWarnings(mlest <- optim(theta, nll_factor_bound, method = optim_method,
                                        y = y, k = k, nlev = nlev, sigh = sigh[i]))
        
        pd<-nll_factor_bound_pd(mlest$par,y,k,nlev,sigh = sigh[i])
        scale<-1/abs(pd)
        
        # Step 2: Refined estimation
        mlest2 <- optim(mlest$par, nll_factor_bound, method = optim_method,
                        y = y, k = k, nlev = nlev, sigh = sigh[i],
                        hessian = TRUE)
      }
      
      
      # Step 3: Compute AIC for boundary model
      
      mlest2$value*-1}, error = function(e) NA)
    
  }
  
  
  # Step 4: Visualization (if plot = TRUE)
  if (plot) {
    
    par(mar=c(5,5,4,2))
    
    plot(sigh,likelihood,xlab="Measurement Error", ylab="Likelihood", pch=16)
    
  }
  
  return(list(Measurement_error=sigh, Likelihood = likelihood))
  
  
}



# 9. Outlier Removal============================================================

outlier_remove <- function(data, y, d=1.5){
  
  Q1 <- quantile(data$y, 0.25, na.rm = TRUE)
  Q3 <- quantile(data$y, 0.75, na.rm = TRUE)
  IQR <- Q3-Q1
  lower_bound <- Q1 - d*IQR
  upper_bound <- Q3 + d*IQR
  
  data[data$y >= lower_bound &  data$y <= upper_bound ,]
}



# 10. Standard error from Hessian ==============================================

seHessian<-function(a, hessian = FALSE, silent = FALSE){
  namesp <- colnames(a)
  mathessian <- a
  mathessian <- ifelse(mathessian == -Inf, -1e+09, mathessian)
  mathessian <- ifelse(mathessian == +Inf, 1e+09, mathessian)
  sigma <- try(solve(mathessian), silent = TRUE)
  if (inherits(sigma, "try-error")) {
    if (!silent)
      warning("Error in Hessian matrix inversion")
    mathessianx <- try(as.matrix(getFromNamespace("nearPD",
                                                  ns = "Matrix")(mathessian)$mat), silent = TRUE)
    if (inherits(mathessianx, "try-error")) {
      if (!silent)
        warning("Error in estimation of the Nearest Positive Definite Matrix. Calculates the Moore-Penrose generalized inverse. Use result with caution.")
      sigma <- try(ginv(mathessian), silent = TRUE)
      if (is.null(colnames(sigma)) | is.null(rownames(sigma))) {
        colnames(sigma) <- rownames(sigma) <- colnames(mathessian)
      }
    }
    else {
      if (!silent)
        warning("Calculates the Nearest Positive Definite Matrix. Use result with caution.")
      sigma <- try(solve(mathessianx), silent = TRUE)
    }
  }
  if (!inherits(sigma, "try-error")) {
    if (all(diag(sigma) >= 0)) {
      res <- sqrt(diag(sigma))
    }
    else {
      s. <- svd(sigma)
      R <- t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
      res <- structure((matrix(rep(1, nrow(R)), nrow = 1,
                               byrow = TRUE) %*% R)[1, ], .Names = colnames(mathessian))
      if (any(res < 0)) {
        d <- diag(as.matrix(getFromNamespace("nearPD",
                                             ns = "Matrix")(sigma)$mat))
        names(d) <- colnames(mathessian)
        res <- ifelse(d < 0, NA, sqrt(d))
      }
      if (any(is.na(res))) {
        a <- sigma
        n = dim(a)[1]
        root = matrix(0, n, n)
        for (i in 1:n) {
          sum = 0
          if (i > 1) {
            sum = sum(root[i, 1:(i - 1)]^2)
          }
          x = a[i, i] - sum
          if (x < 0) {
            x = 0
          }
          root[i, i] = sqrt(x)
          if (i < n) {
            for (j in (i + 1):n) {
              if (root[i, i] == 0) {
                x = 0
              }
              else {
                sum = 0
                if (i > 1) {
                  sum = root[i, 1:(i - 1)] %*% t(t(root[j,
                                                        1:(i - 1)]))
                }
                x = (a[i, j] - sum)/root[i, i]
              }
              root[j, i] = x
            }
          }
        }
        colnames(root) <- rownames(root) <- colnames(mathessian)
        pseudoV <- root %*% t(root)
        d <- diag(pseudoV)
        if (any(d != 0) & all(d >= 0)) {
          res <- sqrt(d)
          if (!silent)
            warning("Estimates using pseudo-variance based on Gill & King (2004)")
        }
        else {
          if (!silent)
            warning("Approximation of Cholesky matrix based on Rebonato and Jackel (2000)")
          if (!silent)
            warning("Estimates using pseudo-variance based on Gill & King (2004)")
          newMat <- a
          cholError <- TRUE
          iter <- 0
          while (cholError) {
            iter <- iter + 1
            newEig <- eigen(newMat)
            newEig2 <- ifelse(newEig$values < 0, 0, newEig$values)
            newMat <- newEig$vectors %*% diag(newEig2) %*%
              t(newEig$vectors)
            newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))
            cholStatus <- try(u <- chol(newMat), silent = TRUE)
            cholError <- ifelse(inherits(cholStatus,
                                         "try-error"), TRUE, FALSE)
          }
          root <- cholStatus
          colnames(root) <- rownames(root) <- colnames(mathessian)
          pseudoV <- root %*% t(root)
          res <- sqrt(diag(pseudoV))
        }
      }
    }
  }
  SEInf <- namesp[!namesp %in% names(res)]
  res <- c(res, structure(rep(+Inf, length(SEInf)), .Names = SEInf))
  if (hessian) {
    return(list(SE = res, hessian = mathessian))
  }
  else {
    return(res)
  }
}


# 11. Categorical boundary

ccnm <- function(x, y, sigh, theta=NULL, optim_method = "L-BFGS-B", plot = FALSE,nullplot=F,
                       fit_BL=FALSE, equation=NULL, start=NULL, show_Model=TRUE,
                       xlab="x", ylab="y",legend_pos="topright",...) {
  
  # Convert x to factor if not already
  Model<-"categorical"
  k <- as.factor(x)
  nlev <- length(unique(k))  # Number of levels in the factor
  
  # Define the negative log-likelihood for the bounded model
  nll_factor_bound <- function(theta, y, k, nlev, sigh) {
    mu <- theta[1:nlev]
    bound <- theta[(nlev+1):(2*nlev)]
    sd <- theta[(2*nlev+1)]
    
    n <- length(y)
    llik <- 0
    
    # Calculate max response per level. CM addition
    max_y_per_level <- tapply(y, k, max, na.rm = TRUE)
    
    for (i in 1:n) {
      
      i.lev <- as.numeric(k[i])
      
      # # Constraint: Ensure bound does not exceed max(y) for that level. CM addition
      # if (bound[i.lev] > max_y_per_level[i.lev]) {
      #   return(Inf)  # Penalize invalid solutions
      # }
      
      
      dens <- coffcturb(y[i], mu[i.lev], sd, -Inf, bound[i.lev], sigh)
      llik <- llik + log(dens)
    }
    
    return(-1 * llik)
  }
  
  # Approximate partial derivatives for scaling boundary model for factor
  nll_factor_bound_pd<-function(theta,y,k,nlev,sigh){
    
    eps=1e-4
    nr<-length(theta)
    part<-vector("numeric",nr)
    
    for (i in 1:nr){
      del<-rep(0,nr)
      del[i]<-eps
      part[i]<-(nll_factor_bound((theta+del),y,k,nlev,sigh)-
                  nll_factor_bound(theta,y,k,nlev,sigh))/eps
    }
    
    return(part)
  }
  
  # Define Turban's function for censored normal distribution with measurement error
  coffcturb <- function(x, mu, sig, a, c, sigh) {
    k <- ((mu - c) / sig)
    d <- ((mu - a) / sig)
    alpha <- ((sigh * sigh) * (x - mu)) / ((sigh * sigh) + (sig * sig))
    beta <- sqrt((sigh * sigh * sig * sig) / ((sigh * sigh) + (sig * sig)))
    gamma <- (beta * sqrt(2 * pi)) / (2 * pi * sig * sigh * (pnorm(d) - pnorm(k)))
    
    com1 <- -((x - mu)^2) / (2 * ((sigh * sigh) + (sig * sig)))
    com2 <- gamma * exp(com1)
    f <- com2 * (pnorm((x - a - alpha) / beta) - pnorm((x - c - alpha) / beta))
    
    # Rescale for censored values
    f <- f * pnorm(c, mu, sig)
    
    # Add contribution at x from mass at c
    f <- f + (dnorm((x - c), 0, sigh) * (1 - pnorm(c, mu, sig)))
    
    return(f)
  }
  
  ## The null-model ----------------------------------------------------------------------
  
  #  Function for negative log likelihood for factor, unbounded
  nll_factor<-function(theta,y,k,nlev){
    
    nlev<-length(theta)-1
    mu<-theta[1:nlev]
    sd<-theta[(nlev+1)]
    n<-length(y)
    llik<-0
    
    for(i in 1:n){
      i.lev<-as.numeric(k[i])
      dens<-dnorm(y[i],mu[i.lev],sd,log=F)
      llik<-llik+log(dens)
    }
    return(-1*llik)
  }
  
  # Function to approximate partial derivatives for scaling boundary model for factor
  nll_factor_pd<-function(theta,y,k,nlev){
    eps=1e-4
    nr<-length(theta)
    part<-vector("numeric",nr)
    
    for (i in 1:nr){
      del<-rep(0,nr)
      del[i]<-eps
      part[i]<-(nll_factor((theta+del),y,k,nlev)-
                  nll_factor(theta,y,k,nlev))/eps
    }
    return(part)
  }
  
  ## Optimization ------------------------------------------------------------------------
  
  # Step 1: Compute scaling for optimization----------------------------------------------
  
  
  ## theta to use-------------------------------------------------------------------
  
  if (is.null(theta)) {
    
    means <- tapply(y,x,"mean")
    upper <- tapply(y, x, function(z) quantile(z, 0.99))
    lmod<-lm(y~x)
    res <- sd(lmod$residuals)
    
    theta <- unname(c(means,upper,res))
    
  }else{
    
    theta <- theta
  }
  
  ## Calculate approximate partial derivatives at the guess for scaling
  
  if(optim_method=="L-BFGS-B"){
    
    # Define lower and upper bounds for parameters. CM edits
    lower_bounds <- c(rep(-Inf, nlev), rep(min(y), nlev), 1e-100)  # Boundaries cannot be below min(y)
    upper_bounds <- c(rep(Inf, nlev), tapply(y, k, max, na.rm = TRUE), Inf)  # Ensure upper bound ≤ max(y)
    
    
    pd<-nll_factor_bound_pd(theta,y,k,nlev,sigh)
    scale<-1/abs(pd)
    
    # first estimate
    suppressWarnings(mlest <- optim(theta, nll_factor_bound, method = optim_method,
                                    y = y, k = k, nlev = nlev, sigh = sigh,
                                    lower = lower_bounds, upper = upper_bounds))
    
    pd<-nll_factor_bound_pd(mlest$par,y,k,nlev,sigh)
    scale<-1/abs(pd)
    
    # Step 2: Refined estimation
    mlest2 <- optim(mlest$par, nll_factor_bound, method = optim_method,
                    y = y, k = k, nlev = nlev, sigh = sigh,,
                    lower = lower_bounds, upper = upper_bounds, hessian = TRUE)
    
  }else{
    pd<-nll_factor_bound_pd(theta,y,k,nlev,sigh)
    scale<-1/abs(pd)
    
    # first estimate
    suppressWarnings(mlest <- optim(theta, nll_factor_bound, method = optim_method,
                                    y = y, k = k, nlev = nlev, sigh = sigh))
    
    pd<-nll_factor_bound_pd(mlest$par,y,k,nlev,sigh)
    scale<-1/abs(pd)
    
    # Step 2: Refined estimation
    mlest2 <- optim(mlest$par, nll_factor_bound, method = optim_method,
                    y = y, k = k, nlev = nlev, sigh = sigh,
                    hessian = TRUE)
  }
  
  
  # Extract bounded model estimates
  
  bounded_model <- mlest2$par
  
  # Determine standard error of bounded values
  
  stderr <- seHessian(mlest2$hessian, hessian = FALSE, silent = FALSE)[(nlev+1) : (nlev*2)]
  
  # Add results in a matrix for parameter and null model
  
  lables <- levels(k)
  l<-2*nlev
  Parameterz <- matrix(bounded_model[1:l], nrow = nlev, ncol = 2)
  Parameterz <-cbind(Parameterz,stderr)
  rownames(Parameterz) <- lables
  colnames(Parameterz) <- c("Mean","Censor","std error")
  
  ###null parameters
  
  
  # Step 3: Compute AIC for boundary model
  AIC_b <- 2 * (mlest2$value + 2 * nlev + 2)
  
  # Step 4: Visualization (if plot = TRUE)
  if (plot) {
    par(mfrow = c(1, 1))
    
    # Strip chart with jittered points
    stripchart(y ~ k, method = "jitter", jitter = 0.1, pch = 16, 
               col = "black", cex = 0.6, vertical = TRUE, 
               ylab = ylab, xlab = xlab,...)
    
    # Add estimated upper bounds as red horizontal lines and solid points
    for (i in 1:nlev) {
      
      #lines(c(i - 0.15, i + 0.15), rep(bounded_model[nlev + i], 2), col = "red", lwd = 2)
      
      points(i,bounded_model[nlev + i], col = "red", pch = 16, cex=1.3)
      
      #lines(c(i, i), c(Parameterz[,2][i]-Parameterz[,3][i], Parameterz[,2][i]+Parameterz[,3][i]),col = "red", lwd = 2)
      
      arrows(i, Parameterz[,2][i]-Parameterz[,3][i], i, Parameterz[,2][i]+Parameterz[,3][i], 
             angle=90, code=3, length=0.1, col = "red", lwd = 2)  # Vertical line with horizontal ends
      
    }
    
    
    if(fit_BL==TRUE & show_Model==TRUE){
      
      legend(legend_pos, legend = c("Data", "Estimated bounds", "Standrd error","Fitted Model"), 
             col = c("black", "red", "red","blue"), pch = c(16,16, NA, NA), lty = c(NA,NA, 1, 2), lwd = c(NA,NA, 2, 2))
      
    } else if(nullplot==TRUE){
      
      legend(legend_pos, legend = c("Data", "Estimated bounds", "Null model mean", "Standrd error", "Standrd error null"), 
             col = c("black", "red", "blue", "red", "blue"), pch = c(16, 16, 16,NA,NA), 
             lty = c(NA,NA,NA,1, 1), lwd = c(NA,NA,NA,2, 2))
      
    }else{
      
      legend(legend_pos, legend = c("Data", "Estimated bounds", "Standrd error"), 
             col = c("black", "red", "red"), pch = c(16, 16, NA), lty = c(NA,NA, 1), lwd = c(NA,NA, 2))
    }
    
  }
  
  
  ## The NULL model---------------------------------------------------------------------
  
  theta_nb_guess<-c(mlest$par[1:nlev],mlest$par[(2*nlev+1)])
  
  # Calculate approximate partial derivatives at the guess for scaling
  
  pd<-nll_factor_pd(theta_nb_guess,y,k,nlev)
  scale<-1/abs(pd)
  
  # first estimate
  
  suppressWarnings(null.est<-optim(theta_nb_guess,nll_factor,method="BFGS",
                                   y=y,k=k,nlev=nlev,
                                   control = list(parscale = scale)))
  
  #refined estimate
  
  null.est2<-optim(null.est$par,nll_factor,method="BFGS",
                   y=y,k=k,nlev=nlev,
                   control = list(parscale = scale),
                   hessian="T")
  
  AIC_null<-2*(null.est2$value+nlev+1)
  
  unbounded_model <- null.est2$par[1:nlev]
  
  stderr_mean <- seHessian(null.est2$hessian, hessian = FALSE, silent = FALSE)[1 : nlev]
  
  lables <- levels(k)
  null_parz <- matrix(unbounded_model, nrow = nlev)
  null_parz <-cbind(null_parz,stderr_mean)
  rownames(null_parz) <- lables
  colnames(null_parz) <- c("Mean","std error")
  
  if(plot==TRUE & nullplot==TRUE){
    
    for (i in 1:nlev) {
      
      points(i,unbounded_model[i], col = "blue", pch = 16, cex=1.3)
      
      arrows(i, null_parz[,1][i]-null_parz[,2][i], i, null_parz[,1][i]+null_parz[,2][i], 
             angle=90, code=3, length=0.1, col = "blue", lwd = 2)  # Vertical line with horizontal ends
      
    }
    
  }
  
  #--------------------------------------------------------------------------------------
  ## AIC values 
  
  AikakeIC <- rbind(AIC_b,AIC_null)
  rownames(AikakeIC)<-c("Bounded","Unbounded")
  colnames(AikakeIC)<-c("")
  
  ### Adding Boundary Model -------------------------------------------------------------
  
  if(fit_BL){
    
    #### Names in start and rearranging them  --------------------------------------------
    
    are_entries_named <- function(vec) {
      # Check if names attribute is not NULL
      if (is.null(names(vec))) {
        return(FALSE)
      }
      
      # Check if all entries have non-NA and non-empty names
      has_valid_names <- all(!is.na(names(vec))) && all(names(vec) != "")
      return(has_valid_names)
    }
    
    if(are_entries_named(start)==TRUE){
      start<-start[order(names(start))]
    } else{
      start<-start
    }
    
    start_names <- names(start) ## keep names for output
    start<-unname(start) # removes names from start
    
    ########## Sum-sq function -----------------------------------------------------------
    
    SUMSQ <- function(start, x, y, equation){
      
      start_list <- as.list(start[1:length(start)])
      names(start_list) <- names(start)[1:(length(start))]
      
      y_pred <- mapply(equation, x=x, MoreArgs = start_list)
      
      return(sum((y-y_pred)^2))
      
      
    }
    
    ########## optimization--------------------------------------------------------------
    
    
    
    params1 <- suppressWarnings(optim(start, SUMSQ, equation = equation, 
                                      x = 1:length(Parameterz[,2]), y=unname(Parameterz[,2]),
                                      method = optim_method, hessian = TRUE))
    
    model_parameters <- matrix(params1$par,length(params1$par))
    rownames(model_parameters) <- start_names
    colnames(model_parameters) <- "Estimate"
    
    ## Boundary line ---------------------------------------------------------------------
    
    if (plot) {
      
      if (show_Model==TRUE){
        
        para <- as.list(params1$par)
        names(para) <- names(start)
        
        bound_pred <- mapply(equation, x=1:length(Parameterz[,2]), MoreArgs = para)
        
        lines(1:length(Parameterz[,2]), bound_pred, col = "blue", lwd = 2, lty = 2)
      }
      
      
    }
    
    # Return All Results------------------------------------------------------------------
    
    return(list(
      Model=Model,
      `Data Properties` = Parameterz,
      `Null Model` =null_parz,
      AIC= AikakeIC,
      `Model parameters` = model_parameters,
      `Residue Sum square`= params1$value
    ))
    
    
  }else{
    
    return(list(
      Model=Model,
      `Data Properties` = Parameterz,
      `Null Model` =null_parz,
      AIC= AikakeIC
    ))
    
  }
  
}


# 12. predict categorical

predict_cat <- function(object,x) {
  
  mat <- object[[1]]
  
  if (is.null(rownames(mat))) { # Ensure the matrix has row names
    stop("The matrix must have row names.")
  }
  
  idx <- match(x, rownames(mat))# Match input names to matrix row names
  response <- mat[idx, 2]# Extract the first column (yield) using matched indices
  return(unname(response))
}


# 12 b. predictBL addition


predictBL2 <-function(object,x){
  
  if(object$Model=="blm"){
    y<-tryCatch(lapply(x,
                       function(a,b) b$Parameters[1,1] + b$Parameters[2,1]*a,
                       b=object),error=function(e) NA)
    return(unlist(y))
  }
  
  if(object$Model=="lp"){
    y<-tryCatch(lapply(x,
                       function(a,b) min(b$Parameters[1,1] + b$Parameters[2,1]*a,b$Parameters[3,1],na.rm = F),
                       b=object),error=function(e) NA)
    return(unlist(y))
  }
  
  if(object$Model=="logistic"|object[[1]]=="logistic"){
    
    y<-tryCatch(lapply(x,
                       function(x,b) b$Parameters[3,1]/(1+exp(b$Parameters[2,1]*(b$Parameters[1,1]-x))),
                       b=object),error=function(e) NA)
    return(unlist(y))
  }
  
  if(object$Model=="inv-logistic"|object[[1]]=="inv-logistic"){
    
    y<-tryCatch(lapply(x,
                       function(x,b) b$Parameters[3,1]- (b$Parameters[3,1]/(1+exp(b$Parameters[2,1]*(b$Parameters[1,1]-x)))),
                       b=object),error=function(e) NA)
    return(unlist(y))
  }
  
  if(object$Model=="logisticND"|object[[1]]=="logisticND"){
    
    y<-tryCatch(lapply(x,
                       function(x,b) b$Parameters[3,1]/(1+(b$Parameters[1,1]*exp(-b$Parameters[2,1]*x))),
                       b=object),error=function(e) NA)
    return(unlist(y))
  }
  
  if(object$Model=="double-logistic"|object[[1]]=="double-logistic"){
    
    y<-tryCatch(lapply(x,
                       function(x,b) {
                         (b$Parameters[3,1]/(1+exp(b$Parameters[2,1]*(b$Parameters[1,1]-x))))-(b$Parameters[4,1]/(1+exp(b$Parameters[6,1]*(b$Parameters[5,1]-x))))
                       },
                       b=object),error=function(e) NA)
    return(unlist(y))
  }
  
  if(object$Model=="qd"){
    
    y<-tryCatch(lapply(x,
                       function(x,b) {b$Parameters[1,1] + b$Parameters[2,1]*x + b$Parameters[3,1]*x^2},
                       b=object),error=function(e) NA)
    return(unlist(y))
  }
  
  if(object$Model=="trapezium"){
    
    b<-object
    yr<-b$Parameters[1,1]+b$Parameters[2,1]*x
    yf<-b$Parameters[4,1]+b$Parameters[5,1]*x
    ym<-rep(b$Parameters[3,1],length(x))
    
    dat<-data.frame(yr,yf,ym)
    y<-apply(dat, 1, min)
    return(y)
    
  }
  
  if(object$Model=="schmidt"){
    
    y<-tryCatch(lapply(x,
                       function(x,b) b$Parameters[3,1]-b$Parameters[1,1]*(x-b$Parameters[2,1])^2,
                       b=object),error=function(e) NA)
    return(unlist(y))
  }
  
  if(object$Model=="mit"){
    
    y<-tryCatch(lapply(x,
                       function(x,b) b$Parameters[3,1]-b$Parameters[1,1]*b$Parameters[2,1]^x,
                       b=object),error=function(e) NA)
    return(unlist(y))
  }
  
  if(object$Model=="other"){
    
    predict_y<-function(x){
      do.call(object$Equation, c(list(x=x),as.list(c(object$Parameters[,1]))))
    }
    
    y <- tryCatch(lapply(x, predict_y),error=function(e) NA)
    return(unlist(y))
  }
  
  if(object$Model=="categorical"){
    
    predict_cat <- function(object,x) {
      mat <- object[[2]]
      if (is.null(rownames(mat))) { # Ensure the matrix has row names
        stop("The matrix must have row names.")
      }
      idx <- match(x, rownames(mat))# Match input names to matrix row names
      response <- mat[idx, 2]# Extract the first column (yield) using matched indices
      return(unname(response))
    }
    
    y = predict_cat(object,x)
    return(y)
  }
}

















































