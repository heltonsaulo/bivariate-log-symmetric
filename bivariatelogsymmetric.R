

####################################################################################
################################### Packages #######################################
####################################################################################

require(ssym)
require(resample)
require(normalp) #### Powerexp
require(maxLik)
require(pearson7) #### Pearson VII
require(VGAM)  ###LAplace
require(Bessel)

####################################################################################
################################## Random numbers ##################################
####################################################################################


rmlogsym.Choleski <- function(n, etas, sigmas, rho,  xi, family) {
  eta1   <- etas[1]
  eta2   <- etas[2]
  sigma1 <- sigmas[1]
  sigma2 <- sigmas[2]
  rho    <- rho
  mu     <- c(log(eta1),log(eta2))
  Sigma  <- matrix(c(sigma1^2,rho*sigma1*sigma2,rho*sigma1*sigma2,sigma2^2), 2,2,byrow = TRUE) 
  
    if(family == "Normal" | family == "Student" | family =="Contnormal" | family=="Powerexp"|
       family == "Hyperbolic" | family == "Slash"){
      # generate n random vectors from MVN(mu, Sigma)
      # dimension is inferred from mu and Sigma
      d <- length(mu)
      Q <- chol(Sigma) # Choleski factorization of Sigma
      Z <- matrix(rvgs(n*d,family,xi), nrow=n, ncol=d)
      X <- Z %*% Q + matrix(mu, n, d, byrow=TRUE)
    }
  
    if(family == "PearsonVII"){
      d <- length(mu)
      Q <- chol(Sigma) 
      Z <- matrix(rpearson7(n*d), nrow=n, ncol=d)
      X <- Z %*% Q + matrix(mu, n, d, byrow=TRUE)
    }
    
  if(family == "Laplace"){
    d <- length(mu)
    Q <- chol(Sigma) 
    Z <- matrix(rlaplace(n*d), nrow=n, ncol=d) ####pacote stats
    X <- Z %*% Q + matrix(mu, n, d, byrow=TRUE)
  }
  
  if(family == "Logistic"){
    d <- length(mu)
    Q <- chol(Sigma) 
    Z <- matrix(rlogis(n*d), nrow=n, ncol=d)
    X <- Z %*% Q + matrix(mu, n, d, byrow=TRUE)
  }

return(exp(X))
  
}

####################################################################################
#################################### BLS density ###################################
####################################################################################

dbivlogsym <- function(t1, 
                       t2, 
                       eta1, 
                       eta2, 
                       sigma1, 
                       sigma2, 
                       rho, 
                       xi, 
                       family,
                       log = FALSE){
  
  if(family!="Normal" & family!="Logistic" & family!="Student" & family!="Laplace" & family!="Powerexp"&
     family!="PearsonVII"& family!="Contnormal"& family!="Slash"&family!="Hyperbolic")
    stop("family of distributions specified by the user is not supported!!",call.=FALSE)
  
  t1_tilde <- log((t1/eta1)^(1/sigma1))
  t2_tilde <- log((t2/eta2)^(1/sigma2))
  pp_rho   <- 1-(rho^2)
  int_cte  <- (t1_tilde^2-2*rho*t1_tilde*t2_tilde+t2_tilde^2)/pp_rho
  
  if(family=="Normal"){
    xi  <- 0
    Zgc <- 2*pi
    pdf <- 1/(t1*t2*sigma1*sigma2*sqrt(pp_rho)*Zgc) * exp(-int_cte/2)
  }
  
  if(family=="Logistic"){ 
    xi  <- 0
    Zgc <- pi/2
    pdf <- 1/(t1*t2*sigma1*sigma2*sqrt(pp_rho)*Zgc) * exp(-int_cte)/((1+exp(-int_cte))^2)
  }  
  
  if(family=="Student"){
    if(xi[1]<=0) stop("nu must be positive!!",call.=FALSE)
    nu <- xi[1]
    Zgc <- (gamma(nu/2)*nu*pi)/gamma((nu+2)/2)
    pdf <- 1/(t1*t2*sigma1*sigma2*sqrt(pp_rho)*Zgc) * (1 + int_cte/nu)^(-(nu+2)/2)
  }
  
  if(family=="Laplace"){
    xi <- 0
    Zgc <- pi
    pdf <- 1/(t1*t2*sigma1*sigma2*sqrt(pp_rho)*Zgc)*besselK(sqrt(2*int_cte), 0)
  }
  
  if(family == "Powerexp"){
    if(xi[1]<=-1  | xi[1]>=1) stop("xi must be within the interval (-1, 1)!!",call.=FALSE)
    xii <- xi[1]
    Zgc <- (2^(xii+1))*(1+xii)*gamma(1+xii)*pi
    pdf <- 1/(t1*t2*sigma1*sigma2*sqrt(pp_rho)*Zgc)*exp((-0.5)*int_cte^(1/(1+xii)))
  }
  
  if(family == "PearsonVII"){
    if(xi[1] <= 1) stop("xi must be greater than 1!!",call.=FALSE)
    if(xi[2] <= 0) stop("theta must be greater than 0!!",call.=FALSE)
    xi1   <- xi[1]
    theta <- xi[2]
    Zgc <- (gamma(xi1-1)*theta*pi)/gamma(xi1)
    pdf <- 1/(t1*t2*sigma1*sigma2*sqrt(pp_rho)*Zgc)*(1+(int_cte/theta))^(-xi1)
    
  }  
  
  
  if(family=="Slash"){
    if(xi[1]<=1) stop("nu must be greater than 1!!",call.=FALSE)
    nu <- xi[1]
    Zgc <- pi/(nu-1) * 2^((3-nu)/2)
    pdf <- 1/(t1*t2*sigma1*sigma2*sqrt(pp_rho)*Zgc) * int_cte^(-(nu+1)/2)*pgamma(int_cte/2, (nu+1)/2) * gamma((nu+1)/2)
  }
  
  if(family == "Contnormal"){
    if(xi[1]<=0  | xi[1]>=1) stop("the extra parameters must be within the interval (0, 1)!!",call.=FALSE)
    if(xi[2]<=0  | xi[2]>=1) stop("the extra parameters must be within the interval (0, 1)!!",call.=FALSE)
    v1 <- xi[1]
    v2 <- xi[2]
    Zgc <- 2*pi*((1/sqrt(v2)) + (1/v1) -1)
    pdf <- 1/(t1*t2*sigma1*sigma2*sqrt(pp_rho)*Zgc)*(sqrt(v2)*exp(-0.5*v2*int_cte))+(((1-v1)/v1)*exp(-0.5*int_cte))
  }
  
  if(family == "Hyperbolic"){
    if(xi[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
    nu2 <- xi[1]
    Zgc <- (2*pi*(nu2+1)*exp(-nu2))/(nu2^2)
    pdf <- 1/(t1*t2*sigma1*sigma2*sqrt(pp_rho)*Zgc)*exp(-nu2*sqrt(1+int_cte))
  }
  
  if (log==TRUE){pdf <-log(pdf)}
  
  return(pdf)
}
 
####################################################################################
#################################### ML estimation #################################
####################################################################################

mle.bivlogsym <- function(data,family,xi,print=TRUE) {
  
  
  if(family!="Normal" & family!="Logistic" & family!="Student" & family!="Laplace" & family!="Powerexp"&
     family!="PearsonVII"& family!="Contnormal"& family!="Slash"&family!="Hyperbolic")
    stop("family of distributions specified by the user is not supported!!",call.=FALSE)
  
  
    
   ## log-likelihood
    loglik <- function(par,data,xi,family) {    
      
      t1     <- data[,1]
      t2     <- data[,2]
      eta1   <- par[1]
      eta2   <- par[2]
      sigma1 <- par[3]
      sigma2 <- par[4]
      rho    <- par[5]
      
      loglikvalue <- sum(log( dbivlogsym(t1=t1, 
                                t2=t2, 
                                eta1 = eta1, 
                                eta2 = eta2, 
                                sigma1 = sigma1, 
                                sigma2 = sigma2, 
                                rho = rho, 
                                xi = xi, 
                                family = family,
                                log = FALSE)  + 1e-300 ) )
      
      return(- loglikvalue)  
    }
    
    ## initial guess
    if(family == "Normal" | family == "Student" | family =="Contnormal" | family=="Powerexp"|
       family == "Hyperbolic" | family == "Slash"){
    
    fit1 <- ssym.l(log(data[,1]) ~ rep(1,length(data[,1])) ,family=family, xi=xi, link.mu = "identity", link.phi = "log")
    fit2 <- ssym.l(log(data[,2]) ~ rep(1,length(data[,2])) ,family=family, xi=xi, link.mu = "identity", link.phi = "log")
     }
  
    
    if(family == "PearsonVII"| family == "Laplace"| family == "Logistic"){
      
      fit1 <- ssym.l(log(data[,1]) ~ rep(1,length(data[,1])) ,family="Normal", xi=xi, link.mu = "identity", link.phi = "log")
      fit2 <- ssym.l(log(data[,2]) ~ rep(1,length(data[,2])) ,family="Normal", xi=xi, link.mu = "identity", link.phi = "log")
      
    } 
    
    
    eta1.guess  <- exp(fit1$theta.mu)
    sigma1.gues <- exp(fit1$theta.phi)
    eta2.guess  <- exp(fit2$theta.mu)
    sigma2.gues <- exp(fit2$theta.phi)  
    rho.guess   <- cor(data[,1],data[,2])
    
    startvalues = c(eta1.guess,eta2.guess,sigma1.gues,sigma2.gues,rho.guess)
    opt <- stats::optim(par = startvalues,
                              fn = loglik,
                              family=family,data=data,xi=xi,
                              method = "Nelder-Mead",
                              hessian = TRUE)

    convv <- 0
    
    if (opt$conv != 0) {
     
      convv <- 1
         
      warning("optimization failed to converge")
    
    }
      
    log.lik.est = -opt$value
    estimates   = opt$par
    n           = nrow(data)
      
    #Criterios para selección de modelos
    AIC   = - 2 * log.lik.est + 2 * (5)
    AICc  = AIC + (2 * (5) * ((5) + 1)) / (n - (5) - 1)
    BIC   = - 2 * log.lik.est + log(n) * (5)
    
    eta1   <- estimates[1]
    eta2   <- estimates[2]
    sigma1 <- estimates[3]
    sigma2 <- estimates[4]
    rho    <- estimates[5]
    

    
    fisher = opt$hessian
    se     = sqrt(diag(solve(fisher)))
    
    confidence.int <- cbind(estimates - qnorm(0.975)*se, estimates + qnorm(0.975)*se)
    colnames(confidence.int) <- c("2.5%","97.5%")
    
    hess = as.matrix(opt$hessian)
    
    zstatpars  = estimates / se
    pvalorpars = 2 * pnorm(abs(zstatpars), lower.tail = F)
    

    tb <- miscTools::coefTable(estimates,se, df=(n-(5)))
    
    ## mahalanobis distance
    t1     <- data[,1]
    t2     <- data[,2]
    t1_tilde_hat <-  log((t1/eta1)^(1/sigma1))
    t2_tilde_hat <-  log((t2/eta2)^(1/sigma2))
    d2_maha <- (t1_tilde_hat^2 - 2*rho*t1_tilde_hat*t2_tilde_hat + t2_tilde_hat^2)/(1-rho^2)
    
    if(print == TRUE){
      cat("\n")
      cat("--------------------------------------------------------------\n")
      cat("             Bivariate log-symmetric distribution             \n")
      cat("--------------------------------------------------------------\n")
      cat("--------------------------------------------------------------\n")
      cat("Maximum Likelihood estimation \n")
      cat("Log-Likelihood:", log.lik.est, "\n")
      cat("AIC:", AIC, "AICc:", AICc, "BIC:", BIC, "\n")
      cat("Number of observations:", n, "\n")
      cat("Family:", family, "\n")
      cat("Xi (fixed):", xi, "\n")
      cat("--------------------------------------------------------------\n")
      cat("Coefficients:\n")
      printCoefmat(tb, signif.stars = TRUE, signif.legend = TRUE, digits = 4)
      cat("--------------------------------------------------------------\n")
    }
    

    
    rval = list(coefficients = estimates, 
                se = se, 
                conf.int = confidence.int,
                pvalor =   pvalorpars,
                Initial.values = startvalues, 
                converged = convv, 
                information.criterions = list(aic = AIC,bic = BIC,aicc=AICc), 
                loglik = log.lik.est,
                n = n, 
                family = family,
                d2_maha = d2_maha,
                Hessian = hess,
                data    = data)
    
    
    return(rval)
    
    
}


