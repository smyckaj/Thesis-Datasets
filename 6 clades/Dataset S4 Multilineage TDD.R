#birth-death model
#################
fit_bdmulti=function (phylolist, tot_timelist, f.lamb, f.mu, lamb_par, mu_par, flist, 
          meth = "Nelder-Mead", cst.lamb = FALSE, cst.mu = FALSE, expo.lamb = FALSE, 
          expo.mu = FALSE, fix.mu = FALSE, dt = 0, cond = "crown") 
{
  #####################
  #a function for fitting a birth-death diversification model across multiple phylogenies simultaneously
  #input is lists of trees, total times and sampling fractions, the rest is identical to fit_bd
  #####################
  #the calculation of total n for aicc
  ntiplist=lapply(phylolist,Ntip)
  nobs=Reduce("+",ntiplist)
  
  if (fix.mu == FALSE) {
    init <- c(lamb_par, mu_par)
    p <- length(init)
    optimLH <- function(init) {
      lamb_par <- init[1:length(lamb_par)]
      mu_par <- init[(1 + length(lamb_par)):length(init)]
      f.lamb.par <- function(t) {
        abs(f.lamb(t, lamb_par))
      }
      f.mu.par <- function(t) {
        abs(f.mu(t, mu_par))
      }
      #the LH is now a sum of likelihoods for different trees, each tree can have a different total age (tot_time) 
      #and sampling proportion (f)
      lhsinpar=function(phylo, tot_time, f){likelihood_bd(phylo, tot_time, f.lamb.par, 
                                            f.mu.par, f, cst.lamb = cst.lamb, cst.mu = cst.mu, 
                                            expo.lamb = expo.lamb, expo.mu = expo.mu, dt = dt, 
                                            cond = cond)}
      likelihoodlist=mapply(lhsinpar, phylolist, tot_timelist, flist)
      LH=Reduce("+",likelihoodlist)
      
      return(-LH)
    }
    temp <- suppressWarnings(optim(init, optimLH, method = meth))
    lamb.par <- temp$par[1:length(lamb_par)]
    mu.par <- temp$par[(1 + length(lamb_par)):length(init)]
    f.lamb.par <- function(t) {
      f.lamb(t, lamb.par)
    }
    f.mu.par <- function(t) {
      f.mu(t, mu.par)
    }
    res <- list(model = "birth death", LH = -temp$value, 
                aicc = 2 * temp$value + 2 * p + (2 * p * (p + 1))/(nobs - 
                                                                     p - 1), lamb_par = lamb.par, mu_par = mu.par, 
                f.lamb = Vectorize(f.lamb.par), f.mu = Vectorize(f.mu.par))
  }
  else {
    init <- c(lamb_par)
    p <- length(init)
    optimLH <- function(init) {
      lamb_par <- init[1:length(lamb_par)]
      f.lamb.par <- function(t) {
        abs(f.lamb(t, lamb_par))
      }
      f.mu.par <- function(t) {
        abs(f.mu(t, mu_par))
      }
      #the LH is now a sum of likelihoods for different trees, each tree can have a different total age (tot_time) 
      #and sampling proportion (f)
      lhsinpar=function(phylo, tot_time, f){likelihood_bd(phylo, tot_time, f.lamb.par, 
                                                       f.mu.par, f, cst.lamb = cst.lamb, cst.mu = cst.mu, 
                                                       expo.lamb = expo.lamb, expo.mu = expo.mu, dt = dt, 
                                                       cond = cond)}
      likelihoodlist=mapply(lhsinpar, phylolist, tot_timelist, flist)
      LH=Reduce("+",likelihoodlist)
      
      return(-LH)
    }
    temp <- suppressWarnings(optim(init, optimLH, method = meth))
    lamb.par <- temp$par[1:length(lamb_par)]
    f.lamb.par <- function(t) {
      f.lamb(t, lamb.par)
    }
    f.mu.par <- function(t) {
      f.mu(t, mu_par)
    }
    res <- list(model = "birth.death", LH = -temp$value, 
                aicc = 2 * temp$value + 2 * p + (2 * p * (p + 1))/(nobs - 
                                                                     p - 1), lamb_par = lamb.par, f.lamb = Vectorize(f.lamb.par))
  }
  class(res) <- "fit.bd"
  return(res)
}


#environmental model
####################
fit_envmulti=function (phylolist, env_data, tot_timelist, f.lamb, f.mu, lamb_par, 
                       mu_par, df = NULL, flist, meth = "Nelder-Mead", cst.lamb = FALSE, 
                       cst.mu = FALSE, expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE, 
                       dt = 0, cond = "crown") 
{
  #####################
  #an environmental spline wrapper for handling fit_bdmulti instead fir_bd
  #####################
  if (is.null(df)) {
    df <- smooth.spline(x = env_data[, 1], env_data[, 2])$df
  }
  spline_result <- sm.spline(env_data[, 1], env_data[, 2], 
                             df = df)
  env_func <- function(t) {
    predict(spline_result, t)
  }
  lower_bound_control <- 0.1
  upper_bound_control <- 0.1
  lower_bound <- min(env_data[, 1])
  upper_bound <- max(env_data[, 1])
  time_tabulated <- seq(from = lower_bound * (1 - lower_bound_control), 
                        to = upper_bound * (1 + upper_bound_control), length.out = 1 + 
                          1e+06)
  env_tabulated <- env_func(time_tabulated)
  env_func_tab <- function(t) {
    b <- upper_bound * (1 + upper_bound_control)
    a <- lower_bound * (1 - lower_bound_control)
    n <- length(env_tabulated) - 1
    index <- 1 + as.integer((t - a) * n/(b - a))
    return(env_tabulated[index])
  }
  f.lamb.env <- function(t, y) {
    f.lamb(t, env_func_tab(t), y)
  }
  f.mu.env <- function(t, y) {
    f.mu(t, env_func_tab(t), y)
  }
  #here we use fit_bdmulti instead of fit_bd
  res <- fit_bdmulti(phylolist, tot_timelist, f.lamb.env, f.mu.env, lamb_par, 
                     mu_par, flist, meth=meth, cst.lamb, cst.mu, expo.lamb, expo.mu, 
                     fix.mu, dt, cond)
  res$model <- "environmental birth death"
  res$f.lamb <- function(t) {
    f.lamb(t, env_func_tab(t), res$lamb_par)
  }
  if (fix.mu == FALSE) {
    res$f.mu <- function(t) {
      f.mu(t, env_func_tab(t), res$mu_par)
    }
  }
  class(res) <- "fit.env"
  return(res)
}