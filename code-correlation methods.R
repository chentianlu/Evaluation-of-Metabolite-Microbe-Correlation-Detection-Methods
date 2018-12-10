#Contents
#R code for SparCC, CCLasso and Cosine similarity
#SparCC
#CCLasso
#Cosine similarity



#SparCC:
  ###############################################################################
# File: SparCC.R
# Aim : SparCC 
#---------------------------------------------------------------------------------------------------------------------
# Author1 : Fang Huaying (Peking University)
# Email  : hyfang@pku.edu.cn
# Date   : 11/12/2014
# Author2 : Tianlu Chen (Shanghai Jiao Tong University)
# Email  : chentianlu@sjtu.edu.cn
# Date   : 02/22/2017
#---------------------------------------------------------------------------------------------------------------------

require(gtools);

#---------------------------------------------------------------------------------------------------------------------
# calculate both correlation r and pseudo p-values for SparCC by permutation
#   function: SparCC.both
#   input:
#         x ------ nxp count data matrix, row is sample, col is variable
#         n_boot ------ Bootstrap times, Default: 20
#   output: a list structure   
#         cor ------ correlation estimation
#         p ------ pseudo p-values, Default: 0.5
#
# call functions SparCC.count and SparCC.frac

SparCC.both <- function (x,n_boot = 20) {
  c <- ncol(x);
  n <- nrow(x);
  cor = p = matrix(0,c, c);
  result <- SparCC.count(x);
  cor <- result$cor.w;
  time_name<-paste("time",1:n_boot,sep="");
  
  for (i in 1:n_boot) {
    
    #   boot_x <- matrix(sample(x,replace = T),n,c,byrow = T);
    boot_x <- sample(x,replace = T);
    boot_result <- SparCC.count(boot_x);
    boot_cor <- boot_result$cor.w;
    #    assign(time_name[i],ifelse(abs(cor) <= abs(boot_cor),1,0));
    temp <-ifelse(abs(cor) <= abs(boot_cor),1,0);
    p<-p+temp
  }
  p <- p /n_boot;
  diag(p) <- 0;
  
  return(list(cor = cor, p = p));
}

#---------------------------------------------------------------------------------------------------------------------
# SparCC for counts known
# function: SparCC.count
# input:
# x ------ nxp count data matrix, row is sample, col is variable
# imax ------ resampling times from posterior distribution. default 20
# kmax ------ max iteration steps for SparCC. default is 10
# alpha ------ the threshold for strong correlation. default is 0.1
# Vmin ------ minimal variance if negative variance appears. default is 1e-4
# output: a list structure
# cov.w ------ covariance estimation
# cor.w ------ correlation estimation

SparCC.count <- function(x, imax = 20, kmax = 10, alpha = 0.1, Vmin = 1e-4) {
  # dimension for w (latent variables)
  p <- as.numeric(ncol(x));
  n <- as.numeric(nrow(x));
  # posterior distribution (alpha)
  x <- x + 1;
  # store generate data
  y <- matrix(0,n,p);
  # store covariance/correlation matrix
  cov.w <- cor.w <- matrix(0, p, p);
  indLow <- lower.tri(cov.w, diag = T);
  # store covariance/correlation for several posterior samples
  covs <- cors <- matrix(0, p * (p + 1) / 2, imax);
  for(i in 1:imax) {
    # generate fractions from posterior distribution
    y <- t(apply(x, 1, function(x) 
      gtools::rdirichlet(n = 1, alpha = x)));
    # estimate covariance/correlation
    cov_cor <- SparCC.frac(x = y, kmax = kmax, alpha = alpha, Vmin = Vmin);
    # store variance/correlation only low triangle 
    covs[, i] <- cov_cor$cov.w[indLow];
    cors[, i] <- cov_cor$cor.w[indLow];
  }
  # calculate median for several posterior samples
  cov.w[indLow] <- apply(covs, 1, median); 
  cor.w[indLow] <- apply(cors, 1, median);
  #
  cov.w <- cov.w + t(cov.w);
  diag(cov.w) <- diag(cov.w) / 2;
  cor.w <- cor.w + t(cor.w);
  diag(cor.w) <- 1;
  #
  return(list(cov.w = cov.w, cor.w = cor.w));
}

#---------------------------------------------------------------------------------------------------------------------
# SparCC for fractions known
# function: SparCC.frac
# input:
#  x ------ nxp fraction data matrix, row is sample, col is variable
#  kmax ------ max iteration steps for SparCC. default is 10
#  alpha ------ the threshold for strong correlation. default is 0.1
#  Vmin ------ minimal variance if negative variance appears. default is 1e-4
#  output: a list structure
#  cov.w ------ covariance estimation
#  cor.w ------ correlation estimation
SparCC.frac <- function(x, kmax = 10, alpha = 0.1, Vmin = 1e-4) {  
  # Log transformation
  x <- log(x);
  p <- ncol(x);
  # T0 = var(log(xi/xj)) variation matrix
  TT <- stats::var(x);
  T0 <- diag(TT) + rep(diag(TT), each = p) - 2 * TT;
  # Variance and correlation coefficients for Basic SparCC  
  rowT0 <- rowSums(T0);
  var.w <- (rowT0 - sum(rowT0) / (2 * p - 2))/(p - 2);
  var.w[var.w < Vmin] <- Vmin;
  #cor.w <- (outer(var.w, var.w, "+") - T0 ) / 
  #  sqrt(outer(var.w, var.w, "*")) / 2;
  Is <- sqrt(1/var.w);
  cor.w <- (var.w + rep(var.w, each = p) - T0) * Is * rep(Is, each = p) * 0.5;
  # Truncated correlation in [-1, 1]
  cor.w[cor.w <= - 1] <- - 1; 
  cor.w[cor.w >= 1] <- 1;
  # Left matrix of estimation equation
  Lmat <- diag(rep(p - 2, p)) + 1; 
  # Remove pairs
  rp <- NULL;
  # Left components
  cp <- rep(TRUE, p);
  # Do loops until max iteration or only 3 components left
  k <- 0;  
  while(k < kmax && sum(cp) > 3) {
    # Left T0 = var(log(xi/xj)) after removing pairs
    T02 <- T0;
    # Store current correlation to find the strongest pair
    curr_cor.w <- cor.w;
    # Remove diagonal
    diag(curr_cor.w) <- 0;
    # Remove removed pairs
    if(!is.null(rp)) {
      curr_cor.w[rp] <- 0;
    }
    # Find the strongest pair in vector form
    n_rp <- which.max(abs(curr_cor.w));
    # Remove the pair if geater than alpha
    if(abs(curr_cor.w[n_rp]) >= alpha) {
      # Which pair in matrix form
      t_id <- c(arrayInd(n_rp, .dim = c(p, p)));
      Lmat[t_id, t_id] <- Lmat[t_id, t_id] - 1;
      # Update remove pairs
      n_rp <- c(n_rp, (p + 1) * sum(t_id) - 2 * p - n_rp);
      rp <- c(rp, n_rp);
      # Update T02
      T02[rp] <- 0;
      # Which component left
      cp <- (diag(Lmat) > 0);
      # Update variance and truncated lower by Vmin
      var.w[cp] <- solve(Lmat[cp, cp], rowSums(T02[cp, cp]));
      var.w[var.w <= Vmin] <- Vmin;
      # Update correlation matrix and truncated by [-1, 1]
      #cor.w <- (outer(var.w, var.w, "+") - T0 ) / 
      #  sqrt(outer(var.w, var.w, "*")) / 2;    
      Is <- sqrt(1/var.w);
      cor.w <- (var.w + rep(var.w, each = p) - T0) * 
        Is * rep(Is, each = p) * 0.5;
      # Truncated correlation in [-1, 1]
      cor.w[cor.w <= - 1] <- - 1;
      cor.w[cor.w >= 1] <- 1;
    }
    else {
      break;
    }
    # 
    k <- k + 1;
  }
  # Covariance
  Is <- sqrt(var.w);
  cov.w <- cor.w * Is * rep(Is, each = p);
  #
  return(list(cov.w = cov.w, cor.w = cor.w));
}
#---------------------------------------------------------------------------------------------------------------------

































CCLasso:
  #function -CCLasso correlation test#
  ########################################################################
# File: cclasso.R
# Aim : Correlation inference for compositional data through lasso
#---------------------------------------------------------------------------------------------------------------------
# Author : Fang Huaying (Peking University)
# Email  : hyfang@pku.edu.cn
# Date   : 2016-01-08
# Version: 2.0
#---------------------------------------------------------------------------------------------------------------------
# Main function: cclasso(x, counts = FALSE, pseudo = 0.5, k_cv = 3, 
#                        lam_int = c(1e-4, 1), k_max = 20, n_boot = 20) 
#
#  Input:
#  x ------ n x p data matrix (row/column is sample/variable)
#        n samples & p compositional variables
#  counts ------ Is the compositional data matrix a count matrix? 
#            Default: FALSE
#  pseudo ------ pseudo count if counts = TRUE
#            Default: 0.5
#  k_cv ------ folds of cross validation
#            Default: 3     
#  lam_int ------ tuning parameter interval
#            Default: [1e-4, 1]
#  k_max ------ maximum iterations for golden section method
#            Default: 20
#  n_boot ------ Bootstrap times
#            Default: 20
#  Output: 
#  A list structure contains:
#  var_w ------ variance estimation
#  cor_w ------ correlation estimation
#  p_vals ------ p-values for elements of cor_w equal 0 or not
#  lambda ------ final tuning parameter
#  info_cv ------ information for cross validation
#---------------------------------------------------------------------------------------------------------------------
cclasso <- function(x, counts = FALSE, pseudo = 0.5, k_cv = 3, 
                    lam_int = c(1e-4, 1), k_max = 20, n_boot = 20) {
  n <- nrow(x);
  p <- ncol(x);
  
  if(counts) {
    x <- x + pseudo;
    x <- x / rowSums(x);
  }
  x <- log(x);
  vx2 <- stats::var(x);
  
  # Diagonal weight for loss function
  rmean_vx2 <- rowMeans(vx2);
  wd <- 1/diag(vx2 - rmean_vx2 - rep(rmean_vx2, each = p) + mean(rmean_vx2));
  wd2 <- sqrt(wd);
  
  # Some global parameters for optimization with single lambda
  rho <- 1;
  u_f <- eigen(diag(p) - 1/p)$vectors;
  wd_u <- (t(u_f) %*% (wd * u_f))[-p, -p];
  wd_u_eig <- eigen(wd_u);
  d0_wd <- 1 / ( (rep(wd_u_eig$values, each = p-1) + wd_u_eig$values) / 
                   (2 * rho) + 1 );
  u0_wd <- wd_u_eig$vectors;
  
  # Golden section method for the selection of lambda (log10 scale)
  sigma <- vx2;
  lam_int2 <- log10(range(lam_int));
  a1 <- lam_int2[1]; 
  b1 <- lam_int2[2];
  # Store lambda and corresponding cross validation's loss
  lams <- NULL; 
  fvals <- NULL;
  # Two trial points in first 
  a2 <- a1 + 0.382 * (b1 - a1); 
  b2 <- a1 + 0.618 * (b1 - a1);
  fb2 <- cv_loss_cclasso(lambda2 = 10^b2 / rho, x = x, k_cv = k_cv, 
                         sigma = sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
                         wd2 = wd2);
  lams <- c(lams, b2); 
  fvals <- c(fvals, fb2$cv_loss);
  fa2 <- cv_loss_cclasso(lambda2 = 10^a2 / rho, x = x, k_cv = k_cv, 
                         sigma = fb2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
                         wd2 = wd2);
  lams <- c(lams, a2); 
  fvals <- c(fvals, fa2$cv_loss);
  # Error tolerance for convergence
  err_lam2 <- 1e-1 * max(1, lam_int2);
  err_fval <- 1e-4;
  
  err <- b1 - a1;
  k <- 0;
  while(err > err_lam2 && k < k_max) {
    fval_max <- max(fa2$cv_loss, fb2$cv_loss);
    
    if(fa2$cv_loss > fb2$cv_loss) {
      a1 <- a2;      
      a2 <- b2;
      fa2 <- fb2;
      b2 <- a1 + 0.618 * (b1 - a1);
      fb2 <- cv_loss_cclasso(lambda2 = 10^b2 / rho, x = x, k_cv = k_cv, 
                             sigma = fa2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
                             wd2 = wd2);
      
      lams <- c(lams, b2); 
      fvals <- c(fvals, fb2$cv_loss);
    } else {
      b1 <- b2;
      b2 <- a2;
      fb2 <- fa2;
      a2 <- a1 + 0.382 * (b1 - a1);
      fa2 <- cv_loss_cclasso(lambda2 = 10^a2 / rho, x = x, k_cv = k_cv, 
                             sigma = fb2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
                             wd2 = wd2);
      
      lams <- c(lams, a2); 
      fvals <- c(fvals, fa2$cv_loss);
    }
    fval_min <- min(fa2$cv_loss, fb2$cv_loss);      
    
    k <- k + 1;
    err <- b1 - a1;
    if(abs(fval_max - fval_min) / (1 + fval_min) <= err_fval) {
      break;
    }
  }
  info_cv <- list(lams = lams, fvals = fvals, k = k + 2, 
                  lam_int = 10^c(a1, b1)); 
  if(a1 == lam_int2[1] || b1 == lam_int2[2]) {
    cat("WARNING:\n", "\tOptimal lambda is near boundary! ([", 10^a1, ",", 
        10^b1, "])\n", sep = "");
  }
  
  lambda <- 10^((a2 + b2)/2);
  # Bootstrap for cclasso
  lambda2 <- lambda / rho;
  info_boot <- boot_cclasso(x = x, sigma = fb2$sigma, lambda2 = lambda2, 
                            n_boot = n_boot, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
  
  return(list(var_w = info_boot$var_w, cor_w = info_boot$cor_w, 
              p_vals = info_boot$p_vals, lambda = lambda, info_cv = info_cv));
}
#---------------------------------------------------------------------------------------------------------------------
# Bootstrap for cclasso
boot_cclasso <- function(x, sigma, lambda2, n_boot = 20, 
                         wd, u_f, u0_wd, d0_wd) {
  n <- nrow(x);
  p <- ncol(x);
  
  # Store the result of bootstrap
  cors_boot <- matrix(0, nrow = p * (p - 1)/2, ncol = n_boot + 1);
  vars_boot <- matrix(0, nrow = p, ncol = n_boot + 1);
  cors_mat <- matrix(0, p, p);
  ind_low <- lower.tri(cors_mat);
  
  # Bootstrap procedure
  sam_boot <- matrix(sample(1:n, size = n * n_boot, replace = T), 
                     ncol = n_boot);
  for(k in 1:n_boot) {
    ind_samp <- sam_boot[, k];
    sigma2 <- cclasso_sub(sigma = sigma, vx = var(x[ind_samp, ]), 
                          lambda2 = lambda2, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
    
    vars_boot[, k] <- diag(sigma2);
    Is <- 1 / sqrt(vars_boot[, k]);
    cors_mat <- Is * sigma2 * rep(Is, each = p);
    cors_boot[, k] <- cors_mat[ind_low];
  }
  Is <- 1 / sqrt(diag(sigma));
  cors_mat <- sigma * Is * rep(Is, each = p);
  cors_boot[, n_boot + 1] <- cors_mat[ind_low];
  vars_boot[, n_boot + 1] <- diag(sigma);
  
  #------------------------------------------------------------------------------------------------------------------  
  # Variance estimation via bootstrap
  vars2 <- rowMeans(vars_boot);
  #------------------------------------------------------------------------------------------------------------------
  # Correlations' relationship for artificial null sample
  tol_cor <- 1e-3;
  sam_art0 <- matrix(rnorm(n * p), nrow = n) * rep(sqrt(vars2), each = n);
  cors_art0 <- cor(sam_art0)[ind_low];
  sam_art <- sam_art0 - log(rowSums(exp(sam_art0)));
  sigma_art <- cclasso_sub(sigma = sigma, vx = var(sam_art), 
                           lambda2 = lambda2, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
  Is <- 1 / sqrt(diag(sigma_art));
  cors_mat <- Is * sigma_art * rep(Is, each = p);
  cors_art2 <- cors_mat[ind_low];
  # Bias of estimation between absolute data and cclasso of compositional data
  cors0m2 <- log( ((1 + cors_art0) * (1 - cors_art2)) / ((1 + cors_art2) * 
                                                           (1 - cors_art0)) );
  tmp <- abs(cors_art2) >= tol_cor;
  bias02 <- ifelse(sum(tmp), median(abs(cors0m2)[tmp]), 0);
  # Modification of estimation for cclasso
  cors2 <- log( (1 + cors_boot) / (1 - cors_boot) );    
  cors2mod <- (cors_boot >= tol_cor) * (cors2 + bias02) + 
    (cors_boot <= -tol_cor) * (cors2 - bias02);
  cors2mod <- 1 - rowMeans(2 / (exp(cors2mod) + 1));
  cors2_mat <- diag(p);
  cors2_mat[ind_low] <- cors2mod;
  cors2_mat <- t(cors2_mat);
  cors2_mat[ind_low] <- cors2mod;
  # P-values with null distribution of correlation estimations of absolute data
  p_vals <- pt(cors2mod * sqrt((n - 2) / (1 - cors2mod^2)), df = n - 2);
  p_vals <- ifelse(p_vals <= 0.5, p_vals, 1 - p_vals);
  pval_mat <- diag(p);
  pval_mat[ind_low] <- p_vals;
  pval_mat <- t(pval_mat);
  pval_mat[ind_low] <- p_vals;
  #---------------------------------------------------------------------------
  
  return(list(var_w = vars2, cor_w = cors2_mat, p_vals = pval_mat));    
}
#-------------------------------------------------------------------------------
# cross validation's loss of cclasso for single lambda
cv_loss_cclasso <- function(lambda2, x, k_cv, sigma, 
                            wd, u_f, u0_wd, d0_wd, wd2) {
  n <- nrow(x);
  p <- ncol(x);
  
  n_b <- floor(n / k_cv);
  cv_loss <- 0;
  for(k in 1:k_cv) {
    itest <- (n_b * (k - 1) + 1):(n_b * k);
    vxk <- stats::var(x[itest, ]);
    vx2k <- stats::var(x[-itest, ]);
    
    sigma <- cclasso_sub(sigma = sigma, vx = vx2k, lambda2 = lambda2, 
                         wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
    
    dsig <- sigma - vxk;
    tmp <- rowMeans(dsig);
    dsig <- dsig - tmp - rep(tmp, each = p) + mean(tmp);
    cv_loss <- cv_loss + base::norm(wd2 * dsig, "F")^2;
  }
  
  return(list(cv_loss = cv_loss, sigma = sigma));
}
#-------------------------------------------------------------------------------
# cclasso for single lambda
cclasso_sub <- function(sigma, vx, lambda2, 
                        wd, u_f, u0_wd, d0_wd, 
                        k_max = 200, x_tol = 1e-4) {
  p <- ncol(sigma);
  sigma2 <- sigma;
  LAMBDA <- matrix(0, p, p);
  lambda2 <- matrix(lambda2, p, p);
  diag(lambda2) <- 0;
  
  k <- 0;
  err <- 1;
  while(err > x_tol && k < k_max) {
    # Update sigma
    x_sigma <- t(u_f) %*% ((sigma2 - vx) - LAMBDA) %*% u_f;
    x_sigma[-p,-p] <- u0_wd %*% ((t(u0_wd) %*% x_sigma[-p, -p] %*% u0_wd) * 
                                   d0_wd) %*% t(u0_wd);
    sigma_new <- vx + u_f %*% x_sigma %*% t(u_f);
    # Update sigma2
    A <- LAMBDA + sigma_new;
    sigma2_new <- (A > lambda2) * (A - lambda2) + (A < -lambda2) * 
      (A + lambda2);
    # Update Lambda
    LAMBDA <- LAMBDA + (sigma_new - sigma2_new);
    
    err <- max( abs(sigma_new - sigma)/(abs(sigma) + 1), 
                abs(sigma2_new - sigma2)/(abs(sigma2) + 1) ); 
    k <- k + 1;
    sigma <- sigma_new;
    sigma2 <- sigma2_new;
  }
  
  if(k >= k_max) {
    cat("WARNING of cclasso_sub:\n", "\tMaximum Iteration:", k_max, 
        "&& Relative error:", err, "!\n");
  }
  
  return(sigma);
}
#---------------------------------------------------------------------------------------------------------------------

























# cosine similarity#
cosine=function(x,boot_n=20)
{
  #normalization#
  x<-scale(x, center = TRUE, scale = TRUE)
  #cosine value count#
  x.row<-as.numeric(nrow(x))
  x.col<-as.numeric(ncol(x))
  bb <- matrix(rep(0,x.row^2),x.row,x.row)
  for(j in 1:x.row)
  {
    for(k in 1:x.row)
    {
      if(j<k)
        bb[j,k] = sum(t(x[j,])*x[k,])/sqrt((sum(x[j,]^2))*sum(x[k,]^2))
    }
  }  
  #set up the original p value#
  p<-matrix(rep(0,x.row^2),x.row,x.row)
  #for rountine test#
  for(l in 1:boot_n)
  {
    #bootstrap to generate shuffled datasets#
    boot_xC<-matrix(sample(x,replace = T),x.row,x.col)
    #count the cosine of shuffled datasets#
    bb1 <- matrix(rep(0,x.row^2),x.row,x.row)
    for(m in 1:x.row)
    {
      for(n in 1:x.row)
      {
        if(m<n)
          bb1[m,n]=sum(t(boot_xC[m,])*boot_xC[n,])/sqrt((sum(boot_xC[m,]^2))*sum(boot_xC[n,]^2))
      }
    } 
    #Pseudo p-value calculation#
    tempC<-ifelse(bb1<=bb,1,0)
    p<-tempC+p
  }
  p<-p/boot_n
  #get the right triangle data#
  p[lower.tri(p)]=0
  diag(p)=0
  return(list(r.cosine=bb,p=p))
}

