InverseLinfty <- function(sigma, n, resol=1.5, mu=NULL, maxiter=50, threshold=1e-2, verbose = TRUE) {
  isgiven <- 1;
  if (is.null(mu)){
    isgiven <- 0;
  }
  
  p <- nrow(sigma);
  M <- matrix(0, p, p);
  xperc = 0;
  xp = round(p/10);
  for (i in 1:p) {
    if ((i %% xp)==0){
      xperc = xperc+10;
      if (verbose) {
        print(paste(xperc,"% done",sep="")); }
    }
    if (isgiven==0){
      mu <- (1/sqrt(n)) * qnorm(1-(0.1/(p^2)));
    }
    mu.stop <- 0;
    try.no <- 1;
    incr <- 0;
    while ((mu.stop != 1)&&(try.no<10)){
      last.beta <- beta
      output <- InverseLinftyOneRow(sigma, i, mu, maxiter=maxiter, threshold=threshold)
      beta <- output$optsol
      iter <- output$iter
      if (isgiven==1){
        mu.stop <- 1
      }
      else{
        if (try.no==1){
          if (iter == (maxiter+1)){
            incr <- 1;
            mu <- mu*resol;
          } else {
            incr <- 0;
            mu <- mu/resol;
          }
        }
        if (try.no > 1){
          if ((incr == 1)&&(iter == (maxiter+1))){
            mu <- mu*resol;
          }
          if ((incr == 1)&&(iter < (maxiter+1))){
            mu.stop <- 1;
          }
          if ((incr == 0)&&(iter < (maxiter+1))){
            mu <- mu/resol;
          }
          if ((incr == 0)&&(iter == (maxiter+1))){
            mu <- mu*resol;
            beta <- last.beta;
            mu.stop <- 1;
          }                        
        }
      }
      try.no <- try.no+1
    }
    M[i,] <- beta;
  }
  return(M)
}

InverseLinftyOneRow <- function ( sigma, i, mu, maxiter=50, threshold=1e-2 ) {
  p <- nrow(sigma);
  rho <- max(abs(sigma[i,-i])) / sigma[i,i];
  mu0 <- rho/(1+rho);
  beta <- rep(0,p);
  
  if (mu >= mu0){
    beta[i] <- (1-mu0)/sigma[i,i];
    returnlist <- list("optsol" = beta, "iter" = 0);
    return(returnlist);
  }
  
  diff.norm2 <- 1;
  last.norm2 <- 1;
  iter <- 1;
  iter.old <- 1;
  beta[i] <- (1-mu0)/sigma[i,i];
  beta.old <- beta;
  sigma.tilde <- sigma;
  diag(sigma.tilde) <- 0;
  vs <- -sigma.tilde%*%beta;
  
  while ((iter <= maxiter) && (diff.norm2 >= threshold*last.norm2)){    
    
    for (j in 1:p){
      oldval <- beta[j];
      v <- vs[j];
      if (j==i)
        v <- v+1;    
      beta[j] <- SoftThreshold(v,mu)/sigma[j,j];
      if (oldval != beta[j]){
        vs <- vs + (oldval-beta[j])*sigma.tilde[,j];
      }
    }
    
    iter <- iter + 1;
    if (iter==2*iter.old){
      d <- beta - beta.old;
      diff.norm2 <- sqrt(sum(d*d));
      last.norm2 <-sqrt(sum(beta*beta));
      iter.old <- iter;
      beta.old <- beta;
      if (iter>10)
        vs <- -sigma.tilde%*%beta;
    }
  }
  
  returnlist <- list("optsol" = beta, "iter" = iter)
  return(returnlist)
}

SoftThreshold <- function( x, lambda ) {

  if (x>lambda){
    return (x-lambda);}
  else {
    if (x< (-lambda)){
      return (x+lambda);}
    else {
      return (0); }
  }
}


BBH <- function(X, y, alpha, safety = 0.99) {

  X = as.matrix(X)
  y = as.numeric(y)
  n = nrow(X)
  p = ncol(X)

  Xs = knockoff:::normc(X, center = FALSE)
  colnames(Xs) = colnames(X)

  X.svd = knockoff:::decompose(Xs, randomize = FALSE)

  lambda_min = min(X.svd$d)^2
  s_boundary = min(2*lambda_min, 1)
  s = safety*s_boundary

  s_diff = pmax(0, 2*s - (s/X.svd$d)^2)

  Xk = (
    sweep(X.svd$u, 2, X.svd$d - s/X.svd$d, "*") +
      sweep(X.svd$u_perp, 2, sqrt(s_diff), "*")
  ) %*% t(X.svd$v)
  colnames(Xk) = colnames(X)

  Sigma = crossprod(Xs)
  Sigma = (Sigma + t(Sigma))/2

  D = Sigma - crossprod(Xs, Xk)
  D = (D + t(D))/2

  M = 2*Sigma - D
  M = (M + t(M))/2

  beta1 = as.numeric(solve(M, crossprod(Xs + Xk, y)))
  beta2 = as.numeric(solve(D, crossprod(Xs - Xk, y)))


  fit_sigma = lm.fit(cbind(Xs, Xk), y)
  nu = n - 2*p
  hsigma = sqrt(sum(fit_sigma$residuals^2)/nu)

  invM = chol2inv(chol(M))
  T1 = beta1/(hsigma*sqrt(2*diag(invM)))
  T2 = beta2/(hsigma*sqrt(2/diag(D)))

  p1 = 2*pt(-abs(T1), df = nu)
  p2 = 2*pt(-abs(T2), df = nu)


  ptilde = p2
  ptilde[p1 > sqrt(alpha)] = 1
  ptilde_BH = p.adjust(ptilde, method = "BH")
  selected = which(ptilde_BH <= sqrt(alpha))

  colnames(Xs)[selected]
}



gm  = function(y, x, q=0.1, family="gaussian", ncores = 1,rep.boot=100, var.est=FALSE){

    n = length(y)
    p = dim(x)[2]
    x= scale(x)
    lambda_max = 1*max(abs(t(x) %*% y))/(n^(3/2))
    lambda_min = lambda_max / 1e2
    nlambda=100
    k = (0:(nlambda-1)) / nlambda
    lambda = lambda_max * (lambda_min/lambda_max)^k
    fit0.lasso = cv.glmnet(x, y, family=family, lambda=lambda, nfold=10,alpha=1)
    lambda_ast = fit0.lasso$lambda.min

    coef_ast = coef(fit0.lasso, s='lambda.min')[-1]
    if(sum(coef_ast!=0)>0.9*n){
        ord_idx = order(abs(coef_ast[which(coef_ast!=0)]))
        coef_ast[which(coef_ast!=0)[1:(sum(coef_ast!=0)-0.9*n)]]=0
    }
    sigmaz_vec = numeric(p)
    z = matrix(0, nrow=n, ncol=p)

#################################################################
cat("Generate Gaussian mirrors... ...")

    for(j in 1:p){
        jidx =  intersect(which(coef_ast!=0), c(1:p)[-j])
        xj =  x[,j]
        xnj = x[,jidx]
        xnj.svd = svd(xnj)
        xnj.idx=1:(min(0.95*length(jidx),n))
        pxnj = xnj.svd$u[, xnj.idx] 
        pp= (1 - t(xj) %*% pxnj  %*%  solve(t(pxnj) %*% pxnj) %*% t(pxnj) %*% xj/n)

        if(pp<0.05){
            xnj.idx=1:10
            pxnj = xnj.svd$u[, xnj.idx] 
            sigmaz_vec[j] = (1 - t(xj) %*% pxnj  %*%  solve(t(pxnj) %*% pxnj) %*% t(pxnj) %*% xj/n)/(1- sum(coef_ast!=0)/n)
        }else{
            sigmaz_vec[j] = pp/(1- sum(coef_ast!=0)/n)
        }
        z[,j] = rnorm(n,0,sqrt(sigmaz_vec[j]))
    }


#########################################################################
cat("Calculate mirror statistics ... ...")
    w_vec = numeric(p)
    #cl <- makeCluster(ncores)
    #registerDoParallel(cl)
    para_list = foreach(z_idx_i = 1:p,.packages='glmnet') %do% {

        z_idx = z_idx_i
        z_indicator = numeric(p)
        z_indicator[z_idx]=1
        z_indicator = as.logical(z_indicator)
        x = scale(x)
        xnew = cbind(x[,z_indicator] + z[,z_indicator], x[,z_indicator] - z[,z_indicator], x[,as.logical(1-z_indicator)])

        fit1.lasso = glmnet(xnew, y, family=family, lambda=lambda_ast)
        z_idx_len = length(z_idx)
        b11 = coef(fit1.lasso, s = 'lambda.min')[2:(z_idx_len+1)]    
        b12 = coef(fit1.lasso, s = 'lambda.min')[(z_idx_len+2):(2*z_idx_len+1)]
        w =  abs(b11+b12)- abs(b11-b12)


        pred.lasso = predict(fit1.lasso, newx=xnew )
        residual.lasso = y - pred.lasso 

        if(var.est==TRUE){
            w.boot   =    rep( 0, rep.boot )
            b11.boot = rep(0, rep.boot )
            b12.boot = rep(0, rep.boot )

            for(ind.boot in 1:rep.boot)
            {
              z.boot = matrix(rnorm(n*sum(z_indicator),0,sigmaz_vec[z_idx_i]), nrow=n, ncol= sum(z_indicator) )
              xnew.boot = cbind(x[,z_indicator] + z.boot, x[,z_indicator] - z.boot, x[,as.logical(1-z_indicator)])
              y.boot = sample(residual.lasso, n, replace=T) + pred.lasso
              fit.lasso.boot = glmnet( xnew.boot, y.boot, family="gaussian", lambda=lambda_ast, alpha=1 )
              
              b11.boot[ind.boot] = coef(fit.lasso.boot)[ 2:(z_idx_len+1) ]
              b12.boot[ind.boot] = coef(fit.lasso.boot)[(z_idx_len+2):(2*z_idx_len+1)]
            } ## end the bootstrap with w.boot, a distribution of w_j|z_j.
            
            w_boot = abs( b11.boot + b12.boot) - abs(b11.boot - b12.boot)
        }
        if(var.est==TRUE){
            list(w=w, w_boot=w_boot)
        }else{
            w
        }
    }

    #stopCluster(cl)
    
    if(var.est==TRUE){
        w_vec =  unlist(lapply(para_list,function(x){x[[1]]}))
        w_boot =  array(unlist(lapply(para_list,function(x){x[[2]]})), c(rep.boot, p)) 
        gm.var = gm_var( w_vec, w_boot, q)$variance
    }else{
        w_vec =  unlist(para_list)
    }


    w_vec_nonzero = w_vec[which(w_vec!=0)]
    w_idx = order(w_vec_nonzero, decreasing =T)
    est_fdp = numeric(length(w_vec_nonzero))

    for( i in 1:length(which(w_vec>0))){
        cut = w_vec_nonzero[w_idx[i]]
        est_fdp[i] = sum(w_vec_nonzero < -cut)/(i+1)
    }

    if(length(which(est_fdp-q>0))==0){
        gm_idx= length(which(w_vec>0))
    }else{
        gm_idx = which(est_fdp-q>0)[1]-1
    }

gm_selected = which(w_vec!=0)[w_idx[1:gm_idx]]

if(var.est== TRUE){
    return(results= list(gm_selected=gm_selected, gm_statistics= w_vec, gm_var= gm.var))
}else{
    return(results= list(gm_selected=gm_selected, gm_statistics= w_vec))
}

}
