library(magrittr)
library(reReg)
library(BB)

expit <- function(x) exp(x)/(1+exp(x))

inv <- function (t, z, exa, exb, fun) {
  mapply(t, FUN = function(u) {
    uf <- function(x) u - fun(x, z, exa, exb) ## / Lam.f(10, r, b, model)
    r1 <- dfsane(par = 1, function(y) uf(y), quiet = TRUE)
    if (!r1$convergence) {
      r2 <- exp(spg(par = 1, fn = function(y) uf(exp(y))^2, quiet = TRUE)$par)
      if (uf(r1$par) <= uf(r2)) return(r1$par)
      if (uf(r1$par) > uf(r2)) return(r2)
    } else return(r1$par)
    ## exp(optim(par = 1, fn = function(y) uf(exp(y))^2, 
    ##           control = list(warn.1d.NelderMead = FALSE))$par)
  })
}

# generate recurrent event data
simGSC <- function(n, summary = FALSE, para,
                   xmat, censoring, frailty, tau, origin,
                   Lam0) {
  call <- match.call()
  Z <- frailty
  Cen <- censoring
  X <- xmat
  p <- ncol(X)
  para0 <- list(alpha = rep(0, p), beta = rep(-1, p), eta = rep(0, p), theta = rep(1, p))
  namel <- names(para)
  para0[namel] <- para
  alpha <- para0$alpha
  beta <- para0$beta
  eta <- para0$eta
  theta <- para0$theta
  msg.mismatch <- function(x)
    paste("Parameter", substitute(x), "does not match with the number of covariates.")
  if (length(alpha) != p) stop(msg.mismatch(alpha))
  if (length(eta) != p) stop(msg.mismatch(eta))
  if (length(beta) != p) stop(msg.mismatch(beta))
  if (length(theta) != p) stop(msg.mismatch(theta))
  ## lapply(list(alpha, eta, beta, theta), mismatch, p = p)
  if (missing(Lam0)) {
    Lam <- function(t, z, exa, exb) z * exb * log(1 + t * exa) / exa / .5
    invLam <- function(t, z, exa, exb) (exp(.5 * t * exa / exb / z) - 1) / exa
  } else {
    Lam <- function(t, z, exa, exb) z * Lam0(t * exa) * exb / exa
    invLam <- function(t, z, exa, exb) inv(t, z, exa, exb, Lam)
  }
  if (n != length(origin) & length(origin) > 1)
    stop("Invalid length for 'origin'. See '?simGSC' for details.")
  simOne <- function(id, z, x, cen) {
    y <- min(cen, tau)
    tij <- NULL
    up <- Lam(y, z, c(exp(x %*% alpha)), c(exp(x %*% beta)))
    m <- -1
    while(sum(tij) < up) {
      tij <- c(tij, rexp(1))
      m <- m + 1
    }
    if (m > 0) {
      tij <- invLam(cumsum(tij[1:m]), z, c(exp(x %*% alpha)), c(exp(x %*% beta)))
      return(data.frame(id = id, Time = c(sort(tij), y),
                        event = c(rep(1, m), 0),
                        Z = z, m = m, x = t(x)))
    } else {
      return(data.frame(id = id, Time = y, event = 0,
                        Z = z, m = m, x = t(x)))
    }
  }
  dat <- data.frame(do.call(rbind, lapply(1:n, function(i) simOne(i, Z[i], X[i,], Cen[i]))))
  names(dat)[grep("x.", names(dat))] <- paste0("x", 1:p)
  if (length(origin) > 1) origin <- rep(origin, unlist(lapply(split(dat$id, dat$id), length)))
  dat$t.start <- do.call(c, lapply(split(dat$Time, dat$id), function(x)
    c(0, x[-length(x)]))) + origin
  dat$t.stop <- dat$Time + origin
  dat$Time <- NULL
  if (summary) {
    dg <- min(3, getOption("digits"))
    cat("Call: \n")
    print(call)
    cat("\n")
    cat("Summary:\n")
    cat("Sample size:                                   ", n, "\n")
    cat("Number of recurrent event observed:            ", sum(dat$event), "\n")
    cat("Average number of recurrent event per subject: ", round(sum(dat$event) / n, dg), "\n")
    cat("\n\n")
  }
  dat$m <- dat$Z <- NULL
  ord <- c("id", "t.start", "t.stop", "event")
  dat <- dat[,c(ord, setdiff(names(dat), ord))]
  attr(dat, "Call") <- call
  return(dat)
}

# M is continuous
data_generation <- function(n, a_par, m_par, m_sd, y_par, cum_baseline, tau, censoring = NULL) {
  C1 <- rnorm(n, 0, 1)
  C2 <- rbinom(n, 1, 0.6)
  A <- rbinom(n, 1, expit(cbind(rep(1, n), C1, C2) %*% a_par))
  M <- rnorm(n, mean = cbind(rep(1, n), A, C1, C2) %*% m_par, sd = m_sd)
  if(is.null(censoring)) censoring <- runif(n, 0, 2 * tau)
  dataset <- simGSC(n = n,
                    summary = T,
                    para = list(alpha = c(0,0,0,0), 
                                eta   = c(0,0,0,0), 
                                beta  = y_par,
                                theta = c(0,0,0,0)),
                    xmat = cbind(A, M, C1, C2),
                    censoring = censoring,
                    frailty = rep(1, n),
                    tau = tau,
                    origin = rep(0, n),
                    Lam0 = cum_baseline)
  colnames(dataset) <- c("id", "t.start", "t.stop", "event", "A", "M", "C1", "C2")
  return(dataset)
}

# M is bianry
data_generation2 <- function(n, a_par, m_par, y_par, cum_baseline, tau, censoring = NULL) {
  C1 <- rnorm(n, 0, 1)
  C2 <- rbinom(n, 1, 0.6)
  A <- rbinom(n, 1, expit(cbind(rep(1, n), C1, C2) %*% a_par))
  M <- rbinom(n, 1, expit(cbind(rep(1, n), A, C1, C2) %*% m_par))
  if(is.null(censoring)) censoring <- runif(n, 0, 2 * tau)
  dataset <- simGSC(n = n,
                    summary = T,
                    para = list(alpha = c(0,0,0,0), 
                                eta   = c(0,0,0,0), 
                                beta  = y_par,
                                theta = c(0,0,0,0)),
                    xmat = cbind(A, M, C1, C2),
                    censoring = censoring,
                    frailty = rep(1, n),
                    tau = tau,
                    origin = rep(0, n),
                    Lam0 = cum_baseline)
  colnames(dataset) <- c("id", "t.start", "t.stop", "event", "A", "M", "C1", "C2")
  return(dataset)
}

# data <- data_generation(n = 1e5,
#                         a_par = c(1, 1, -2),
#                         m_par = c(3, -1, -1, 1.5),
#                         m_sd = 2,
#                         y_par = c(1, 0.1, 0.2, -0.22),
#                         cum_baseline = function(x) 0.05 * x,
#                         tau = 24)

# data <- data_generation2(n = 1e3,
#                          a_par = c(1, 1, -2),
#                          m_par = c(3, -4, -3, 4),
#                          y_par = c(1, 0.1, 0.2, -0.22),
#                          cum_baseline = function(x) 0.05 * x,
#                          tau = 24)

