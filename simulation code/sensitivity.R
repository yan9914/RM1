source("fun_sim_misspecified.R")
source("fun_triply_robust.R")
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

data_generation <- function(n, a_par, m_par, m_sd, y_par, cum_baseline, tau) {
  C1 <- rnorm(n, 0, 1)
  C2 <- rbinom(n, 1, 0.6)
  A <- rbinom(n, 1, expit(cbind(rep(1, n), C1, C2) %*% a_par))
  M <- rnorm(n, mean = cbind(rep(1, n), A, C1, C2) %*% m_par, sd = m_sd)
  dataset <- simGSC(n = n,
                    summary = T,
                    para = list(alpha = c(0,0,0,0), 
                                eta   = c(0,0,0,0), 
                                beta  = y_par,
                                theta = c(0,0,0,0)),
                    xmat = cbind(A, M, C1, C2),
                    censoring = ifelse(C2 == 1,
                                       runif(n, 0, tau),
                                       runif(n, 0, 4 * tau)),
                    frailty = rep(1, n),
                    tau = tau,
                    origin = rep(0, n),
                    Lam0 = cum_baseline)
  colnames(dataset) <- c("id", "t.start", "t.stop", "event", "A", "M", "C1", "C2")
  return(dataset)
}

library(data.table)
library(magrittr)
library(survival)
library(parallel)
library(ggplot2)
library(latex2exp)

true_value <- readRDS("true_value_20221011.rds")
true_value_NDE <- stepfun(true_value$time, c(0, true_value$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8))
true_value_NIE <- stepfun(true_value$time, c(0, true_value$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8))

simf <- function(i) {
  set.seed(i)
  data <- data_generation(n = 3e3,
                          a_par = c(1, 1, -2),
                          m_par = c(3, -1, -1, 1.5),
                          m_sd = 2,
                          y_par = c(1, 0.1, 0.2, -0.22),
                          cum_baseline = function(x) 0.05 * x,
                          tau = 24)
  AM <- f_a_misspecified(data, "A", "M", c("C1", "C2"), g_computation = F)$estimate
  MM <- f_m_misspecified(data, "A", "M", c("C1", "C2"), g_computation = F)$estimate
  YM <- f_y_misspecified(data, "A", "M", c("C1", "C2"), g_computation = F)$estimate
  AC <- f(data, "A", "M", c("C1", "C2"), g_computation = F)$estimate
  AM_NDE <- stepfun(AM$time, c(0, AM$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NDE
  AM_NIE <- stepfun(AM$time, c(0, AM$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NIE
  MM_NDE <- stepfun(MM$time, c(0, MM$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NDE
  MM_NIE <- stepfun(MM$time, c(0, MM$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NIE
  YM_NDE <- stepfun(YM$time, c(0, YM$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NDE
  YM_NIE <- stepfun(YM$time, c(0, YM$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NIE
  AC_NDE <- stepfun(AC$time, c(0, AC$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NDE
  AC_NIE <- stepfun(AC$time, c(0, AC$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NIE
  cbind(c("AM", "AM", "MM", "MM", "YM", "YM", "AC", "AC"),
        c("NDE", "NIE", "NDE", "NIE", "NDE", "NIE", "NDE", "NIE"),
        rbind(AM_NDE, AM_NIE, MM_NDE, MM_NIE, YM_NDE, YM_NIE, AC_NDE, AC_NIE))
}

a <- proc.time()
cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores)
clusterExport(cl, c("true_value_NDE",
                    "true_value_NIE",
                    "simf",
                    "f_a_misspecified",
                    "f_m_misspecified",
                    "f_y_misspecified",
                    "f",
                    "data_generation",
                    "expit",
                    "inv",
                    "simGSC"))
clusterEvalQ(cl, c(library(data.table),
                   library(magrittr),
                   library(survival),
                   library(reReg),
                   library(BB)))
result <- parLapply(cl, 1:1000, function(x) simf(x))
stopCluster(cl)
proc.time() - a

result2 <- do.call("rbind", result) %>% 
  as.data.table %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))]
result2$V1 %<>% sapply(function(x) switch(x, "AM" = "A_miss", "MM" = "M_miss", "YM" = "Y_miss", "AC" = "correct"))

# saveRDS(result2, "sensitivity_censor_m_20221031.rds")

result_NIE <- result2[V2 == "NIE",][, .(mean(V3), mean(V4), mean(V5), mean(V6), mean(V7)), by = V1] %>%
  rbind(result2[V2 == "NIE",][, .(median(V3), median(V4), median(V5), median(V6), median(V7)), by = V1]) %>%
  rbind(result2[V2 == "NIE",][, .(sd(V3), sd(V4), sd(V5), sd(V6), sd(V7)), by = V1])
result_NDE <- result2[V2 == "NDE",][, .(mean(V3), mean(V4), mean(V5), mean(V6), mean(V7)), by = V1] %>%
  rbind(result2[V2 == "NDE",][, .(median(V3), median(V4), median(V5), median(V6), median(V7)), by = V1]) %>%
  rbind(result2[V2 == "NDE",][, .(sd(V3), sd(V4), sd(V5), sd(V6), sd(V7)), by = V1])

theme_set(theme(panel.background = element_rect(fill = NA, colour = "black"),
                panel.grid = element_blank() ,
                legend.position = c(0.05,0.05)  ,
                legend.justification = c(0,0),
                legend.background = element_rect(fill = NA, colour = "black", size = 0.25)))

result3 <- melt(result2, id.vars = c("V1", "V2"), measure.vars = c("V3", "V4", "V5", "V6", "V7"))

ggplot(result3) +
  geom_boxplot(aes(V1, value), outlier.size = 0.75) +
  scale_x_discrete(limits = c("correct", "A_miss", "M_miss", "Y_miss"),
                   labels = c("correct" = "all correct",
                              "A_miss" = "A misspecified",
                              "M_miss" = "M misspecified",
                              "Y_miss" = TeX("$\\tilde{N}$ misspecified"))) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust = 0.3)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25) +
  labs(x = "", y = "Bias") +
  ggh4x::facet_grid2(V2~variable, scales = "free",
                     labeller = labeller(.cols = as_labeller(c(V3 = "4.8 months (20%)",
                                                               V4 = "9.6 months (40%)",
                                                               V5 = "12 months (50%)",
                                                               V6 = "14.4 months (60%)",
                                                               V7 = "19.2 months (80%)"))))
