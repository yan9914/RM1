source("fun_data_generation.R")
library(data.table)
library(magrittr)
library(survival)

f_m_misspecified <- function(data, exposure, mediator, confounders, normal_mediator = TRUE, g_computation = FALSE, B_points = NULL) {
  setDT(data)
  subdata <- data[data[, .I[1], by = id]$V1, mget(c("id", exposure, mediator, confounders))]
  n <- nrow(subdata)
  A <- subdata[, get(exposure)] %>% as.matrix
  M <- subdata[, get(mediator)] %>% as.matrix
  C <- subdata[, mget(confounders)] %>% as.matrix
  Y_model <- formula(paste0("Surv(time = t.start, time2 = t.stop, event = event)~",
                            paste0(c(exposure, mediator, confounders, "cluster(id)"), collapse = "+"))) %>%
    coxph(data = data, control = coxph.control(timefix = FALSE))
  baseline <- basehaz(Y_model, centered = F)
  if(normal_mediator) {
    M_model <- lm(reformulate(paste0(c(exposure, confounders), collapse = "+"), mediator), data = subdata)
    M_model$coefficients[confounders] <- M_model$coefficients[confounders] + rnorm(2, 0.8, 0.1)
    M_model$coefficients[exposure] <- M_model$coefficients[exposure] + rnorm(1, -0.2, 0.1)
    if(g_computation) {
      eta_hat <- baseline$hazard * exp(M_model$coefficients["(Intercept)"] * Y_model$coefficients[mediator] + 1/2 * Y_model$coefficients[mediator]^2 * summary(M_model)$sigma^2) *
        mean(exp(C %*% M_model$coefficients[confounders] * Y_model$coefficients[mediator] + C %*% Y_model$coefficients[confounders]))
      if(is.null(B_points)) {
        return(
          list(
            estimate = data.frame(
              time = baseline$time,
              NDE = eta_hat * (exp(Y_model$coefficients[exposure])-1),
              NIE = eta_hat * (exp(Y_model$coefficients[exposure] + Y_model$coefficients[mediator] * M_model$coefficients[exposure]) - exp(Y_model$coefficients[exposure]))
            )
          )
        )
      } else {
        nde <- eta_hat * (exp(Y_model$coefficients[exposure])-1)
        nie <- eta_hat * (exp(Y_model$coefficients[exposure] + Y_model$coefficients[mediator] * M_model$coefficients[exposure]) - exp(Y_model$coefficients[exposure]))
        return(
          c(stepfun(baseline$time,
                    c(0, nde))(B_points),
            stepfun(baseline$time,
                    c(0, nie))(B_points))
        )
      }
      
    }
    g_com <- function(a, b) {
      eta_hat <- (diff(c(0, baseline$hazard)) * exp(M_model$coefficients["(Intercept)"] * Y_model$coefficients[mediator] + 1/2 * Y_model$coefficients[mediator]^2 * summary(M_model)$sigma^2)) %*%
        t(exp(C %*% M_model$coefficients[confounders] * Y_model$coefficients[mediator] + C %*% Y_model$coefficients[confounders]))
      eta_hat * exp(Y_model$coefficients[exposure] * a + Y_model$coefficients[mediator] * M_model$coefficients[exposure] * b)
    }
    A_model <- glm(formula(paste0(exposure, "~", paste(confounders, collapse = "+"))), data = subdata, family = "binomial")
    if(identical(data$t.stop %>% unique %>% length, baseline$time)) stop("time points are not all equal to package(coxph) output")
    at_risk_tb <- data[,(function(x) {tmp <- rep(1L, length(baseline$time)) ; tmp[baseline$time >= x[[1]][x[[2]] == 0]] <- 0L ; return(tmp)})(list(t.stop, event)), by = id] %>%
      .[, V1] %>% matrix(nrow = length(baseline$time))
    dn <- data[,(function(x) {tmp <- rep(0L, length(baseline$time)) ; tmp[baseline$time %in% x[[1]][x[[2]] == 1]] <- 1L ; return(tmp)})(list(t.stop, event)), by = id] %>%
      .[, V1] %>% matrix(nrow = length(baseline$time))
    dn_pred1 <- diff(c(0, baseline$hazard)) %*% exp(t(cbind(1, M, C) %*% Y_model$coefficients))
    dn_pred0 <- diff(c(0, baseline$hazard)) %*% exp(t(cbind(0, M, C) %*% Y_model$coefficients))
    m_ratio <- exp(-1/2/summary(M_model)$sigma^2 * ((M - cbind(1, 0, C) %*% head(M_model$coefficients, 4))^2 - (M - cbind(1, 1, C) %*% head(M_model$coefficients, 4))^2))
    ps <- cbind(1, C) %*% A_model$coefficients %>% expit
    N1 <- ((A == 1) / ps)
    N0 <- ((A == 0) / (1 - ps))
    D1 <- at_risk_tb %*% N1 %>% as.numeric / n
    D0 <- at_risk_tb %*% N0 %>% as.numeric / n
    nan_to_zero <- function(x) {x[is.nan(x)] <- 0 ; return(x)}
    if(is.null(B_points)) {
      return(
        list(
          estimate = data.frame(
            time = baseline$time,
            NDE = 1/n * (
              nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% (m_ratio * N1) / D1) -
                nan_to_zero(((dn - dn_pred0) * at_risk_tb) %*% N0 / D0) +
                nan_to_zero(((dn_pred1 - g_com(1, 0) - dn_pred0 + g_com(0, 0)) * at_risk_tb) %*% N0 / D0) +
                nan_to_zero(apply(g_com(1, 0) - g_com(0, 0), 1, sum))
            ) %>% as.numeric %>% cumsum,
            NIE = 1/n * (
              nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% ((1 - m_ratio) * N1) / D1) +
                nan_to_zero(((dn_pred1 - g_com(1, 1)) * at_risk_tb) %*% N1 / D1) -
                nan_to_zero(((dn_pred1 - g_com(1, 0)) * at_risk_tb) %*% N0 / D0) +
                nan_to_zero(apply(g_com(1, 1) - g_com(1, 0), 1, sum))
            ) %>% as.numeric %>% cumsum
          ),
          M_model = M_model,
          Y_model = Y_model
        )
        
      )
    } else {
      nde <- 1/n * (
        nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% (m_ratio * N1) / D1) -
          nan_to_zero(((dn - dn_pred0) * at_risk_tb) %*% N0 / D0) +
          nan_to_zero(((dn_pred1 - g_com(1, 0) - dn_pred0 + g_com(0, 0)) * at_risk_tb) %*% N0 / D0) +
          nan_to_zero(apply(g_com(1, 0) - g_com(0, 0), 1, sum))
      ) %>% as.numeric %>% cumsum
      nie <- 1/n * (
        nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% ((1 - m_ratio) * N1) / D1) +
          nan_to_zero(((dn_pred1 - g_com(1, 1)) * at_risk_tb) %*% N1 / D1) -
          nan_to_zero(((dn_pred1 - g_com(1, 0)) * at_risk_tb) %*% N0 / D0) +
          nan_to_zero(apply(g_com(1, 1) - g_com(1, 0), 1, sum))
      ) %>% as.numeric %>% cumsum
      return(
        c(stepfun(baseline$time,
                  c(0, nde))(B_points),
          stepfun(baseline$time,
                  c(0, nie))(B_points))
      )
    }
  }
}

f_a_misspecified <- function(data, exposure, mediator, confounders, normal_mediator = TRUE, g_computation = FALSE, B_points = NULL) {
  setDT(data)
  subdata <- data[data[, .I[1], by = id]$V1, mget(c("id", exposure, mediator, confounders))]
  subdata[, conconfound := scale(get(mediator))[,1]]
  n <- nrow(subdata)
  A <- subdata[, get(exposure)] %>% as.matrix
  M <- subdata[, get(mediator)] %>% as.matrix
  C <- subdata[, mget(confounders)] %>% as.matrix
  Y_model <- formula(paste0("Surv(time = t.start, time2 = t.stop, event = event)~",
                            paste0(c(exposure, mediator, confounders, "cluster(id)"), collapse = "+"))) %>%
    coxph(data = data, control = coxph.control(timefix = FALSE))
  baseline <- basehaz(Y_model, centered = F)
  if(normal_mediator) {
    M_model <- lm(formula(paste0(mediator, "~", paste(c(exposure, confounders), collapse = "+"))), data = subdata)
    if(g_computation) {
      eta_hat <- baseline$hazard * exp(M_model$coefficients["(Intercept)"] * Y_model$coefficients[mediator] + 1/2 * Y_model$coefficients[mediator]^2 * summary(M_model)$sigma^2) *
        mean(exp(C %*% M_model$coefficients[confounders] * Y_model$coefficients[mediator] + C %*% Y_model$coefficients[confounders]))
      if(is.null(B_points)) {
        return(
          list(
            estimate = data.frame(
              time = baseline$time,
              NDE = eta_hat * (exp(Y_model$coefficients[exposure])-1),
              NIE = eta_hat * (exp(Y_model$coefficients[exposure] + Y_model$coefficients[mediator] * M_model$coefficients[exposure]) - exp(Y_model$coefficients[exposure]))
            )
          )
        )
      } else {
        nde <- eta_hat * (exp(Y_model$coefficients[exposure])-1)
        nie <- eta_hat * (exp(Y_model$coefficients[exposure] + Y_model$coefficients[mediator] * M_model$coefficients[exposure]) - exp(Y_model$coefficients[exposure]))
        return(
          c(stepfun(baseline$time,
                    c(0, nde))(B_points),
            stepfun(baseline$time,
                    c(0, nie))(B_points))
        )
      }
      
    }
    g_com <- function(a, b) {
      eta_hat <- (diff(c(0, baseline$hazard)) * exp(M_model$coefficients["(Intercept)"] * Y_model$coefficients[mediator] + 1/2 * Y_model$coefficients[mediator]^2 * summary(M_model)$sigma^2)) %*%
        t(exp(C %*% M_model$coefficients[confounders] * Y_model$coefficients[mediator] + C %*% Y_model$coefficients[confounders]))
      eta_hat * exp(Y_model$coefficients[exposure] * a + Y_model$coefficients[mediator] * M_model$coefficients[exposure] * b)
    }
    A_model <- glm(reformulate(paste(c(confounders, "conconfound"), collapse = "+"), exposure), data = subdata, family = "binomial")
    if(identical(data$t.stop %>% unique %>% length, baseline$time)) stop("time points are not all equal to package(coxph) output")
    at_risk_tb <- data[,(function(x) {tmp <- rep(1L, length(baseline$time)) ; tmp[baseline$time >= x[[1]][x[[2]] == 0]] <- 0L ; return(tmp)})(list(t.stop, event)), by = id] %>%
      .[, V1] %>% matrix(nrow = length(baseline$time))
    dn <- data[,(function(x) {tmp <- rep(0L, length(baseline$time)) ; tmp[baseline$time %in% x[[1]][x[[2]] == 1]] <- 1L ; return(tmp)})(list(t.stop, event)), by = id] %>%
      .[, V1] %>% matrix(nrow = length(baseline$time))
    dn_pred1 <- diff(c(0, baseline$hazard)) %*% exp(t(cbind(1, M, C) %*% Y_model$coefficients))
    dn_pred0 <- diff(c(0, baseline$hazard)) %*% exp(t(cbind(0, M, C) %*% Y_model$coefficients))
    m_ratio <- exp(-1/2/summary(M_model)$sigma^2 * ((M - cbind(1, 0, C) %*% M_model$coefficients)^2 - (M - cbind(1, 1, C) %*% M_model$coefficients)^2))
    ps <- cbind(1, C) %*% head(A_model$coefficients, 3) %>% expit
    N1 <- ((A == 1) / ps)
    N0 <- ((A == 0) / (1 - ps))
    D1 <- at_risk_tb %*% N1 %>% as.numeric / n
    D0 <- at_risk_tb %*% N0 %>% as.numeric / n
    nan_to_zero <- function(x) {x[is.nan(x)] <- 0 ; return(x)}
    if(is.null(B_points)) {
      return(
        list(
          estimate = data.frame(
            time = baseline$time,
            NDE = 1/n * (
              nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% (m_ratio * N1) / D1) -
                nan_to_zero(((dn - dn_pred0) * at_risk_tb) %*% N0 / D0) +
                nan_to_zero(((dn_pred1 - g_com(1, 0) - dn_pred0 + g_com(0, 0)) * at_risk_tb) %*% N0 / D0) +
                nan_to_zero(apply(g_com(1, 0) - g_com(0, 0), 1, sum))
            ) %>% as.numeric %>% cumsum,
            NIE = 1/n * (
              nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% ((1 - m_ratio) * N1) / D1) +
                nan_to_zero(((dn_pred1 - g_com(1, 1)) * at_risk_tb) %*% N1 / D1) -
                nan_to_zero(((dn_pred1 - g_com(1, 0)) * at_risk_tb) %*% N0 / D0) +
                nan_to_zero(apply(g_com(1, 1) - g_com(1, 0), 1, sum))
            ) %>% as.numeric %>% cumsum
          ),
          M_model = M_model,
          Y_model = Y_model
        )
        
      )
    } else {
      nde <- 1/n * (
        nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% (m_ratio * N1) / D1) -
          nan_to_zero(((dn - dn_pred0) * at_risk_tb) %*% N0 / D0) +
          nan_to_zero(((dn_pred1 - g_com(1, 0) - dn_pred0 + g_com(0, 0)) * at_risk_tb) %*% N0 / D0) +
          nan_to_zero(apply(g_com(1, 0) - g_com(0, 0), 1, sum))
      ) %>% as.numeric %>% cumsum
      nie <- 1/n * (
        nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% ((1 - m_ratio) * N1) / D1) +
          nan_to_zero(((dn_pred1 - g_com(1, 1)) * at_risk_tb) %*% N1 / D1) -
          nan_to_zero(((dn_pred1 - g_com(1, 0)) * at_risk_tb) %*% N0 / D0) +
          nan_to_zero(apply(g_com(1, 1) - g_com(1, 0), 1, sum))
      ) %>% as.numeric %>% cumsum
      return(
        c(stepfun(baseline$time,
                  c(0, nde))(B_points),
          stepfun(baseline$time,
                  c(0, nie))(B_points))
      )
    }
  }
}

f_y_misspecified <- function(data, exposure, mediator, confounders, normal_mediator = TRUE, g_computation = FALSE, B_points = NULL) {
  setDT(data)
  subdata <- data[data[, .I[1], by = id]$V1, mget(c("id", exposure, mediator, confounders))]
  # subdata[, conconfound := scale(get(mediator))[,1] + get(exposure) * get(confounders[1]) * get(confounders[2])]
  # data[, conconfound := subdata$conconfound[id]]
  n <- nrow(subdata)
  A <- subdata[, get(exposure)] %>% as.matrix
  M <- subdata[, get(mediator)] %>% as.matrix
  C <- subdata[, mget(confounders)] %>% as.matrix
  Y_model <- formula(paste0("Surv(time = t.start, time2 = t.stop, event = event)~",
                            paste0(c(exposure, mediator, confounders, "cluster(id)"), collapse = "+"))) %>%
    coxph(data = data, control = coxph.control(timefix = FALSE))
  baseline <- basehaz(Y_model, centered = F)
  Y_model$coefficients[mediator] <- Y_model$coefficients[mediator] + rnorm(1, 0.05, 0.01)
  Y_model$coefficients[exposure] <- Y_model$coefficients[exposure] + rnorm(1, -0.05, 0.01)
  if(normal_mediator) {
    M_model <- lm(formula(paste0(mediator, "~", paste(c(exposure, confounders), collapse = "+"))), data = subdata)
    if(g_computation) {
      eta_hat <- baseline$hazard * exp(M_model$coefficients["(Intercept)"] * Y_model$coefficients[mediator] + 1/2 * Y_model$coefficients[mediator]^2 * summary(M_model)$sigma^2) *
        mean(exp(C %*% M_model$coefficients[confounders] * Y_model$coefficients[mediator] + C %*% Y_model$coefficients[confounders]))
      if(is.null(B_points)) {
        return(
          list(
            estimate = data.frame(
              time = baseline$time,
              NDE = eta_hat * (exp(Y_model$coefficients[exposure])-1),
              NIE = eta_hat * (exp(Y_model$coefficients[exposure] + Y_model$coefficients[mediator] * M_model$coefficients[exposure]) - exp(Y_model$coefficients[exposure]))
            )
          )
        )
      } else {
        nde <- eta_hat * (exp(Y_model$coefficients[exposure])-1)
        nie <- eta_hat * (exp(Y_model$coefficients[exposure] + Y_model$coefficients[mediator] * M_model$coefficients[exposure]) - exp(Y_model$coefficients[exposure]))
        return(
          c(stepfun(baseline$time,
                    c(0, nde))(B_points),
            stepfun(baseline$time,
                    c(0, nie))(B_points))
        )
      }
      
    }
    g_com <- function(a, b) {
      eta_hat <- (diff(c(0, baseline$hazard)) * exp(M_model$coefficients["(Intercept)"] * Y_model$coefficients[mediator] + 1/2 * Y_model$coefficients[mediator]^2 * summary(M_model)$sigma^2)) %*%
        t(exp(C %*% M_model$coefficients[confounders] * Y_model$coefficients[mediator] + C %*% Y_model$coefficients[confounders]))
      eta_hat * exp(Y_model$coefficients[exposure] * a + Y_model$coefficients[mediator] * M_model$coefficients[exposure] * b)
    }
    A_model <- glm(formula(paste0(exposure, "~", paste(confounders, collapse = "+"))), data = subdata, family = "binomial")
    if(identical(data$t.stop %>% unique %>% length, baseline$time)) stop("time points are not all equal to package(coxph) output")
    at_risk_tb <- data[,(function(x) {tmp <- rep(1L, length(baseline$time)) ; tmp[baseline$time >= x[[1]][x[[2]] == 0]] <- 0L ; return(tmp)})(list(t.stop, event)), by = id] %>%
      .[, V1] %>% matrix(nrow = length(baseline$time))
    dn <- data[,(function(x) {tmp <- rep(0L, length(baseline$time)) ; tmp[baseline$time %in% x[[1]][x[[2]] == 1]] <- 1L ; return(tmp)})(list(t.stop, event)), by = id] %>%
      .[, V1] %>% matrix(nrow = length(baseline$time))
    dn_pred1 <- diff(c(0, baseline$hazard)) %*% exp(t(cbind(1, M, C) %*% head(Y_model$coefficients, 4)))
    dn_pred0 <- diff(c(0, baseline$hazard)) %*% exp(t(cbind(0, M, C) %*% head(Y_model$coefficients, 4)))
    m_ratio <- exp(-1/2/summary(M_model)$sigma^2 * ((M - cbind(1, 0, C) %*% M_model$coefficients)^2 - (M - cbind(1, 1, C) %*% M_model$coefficients)^2))
    ps <- cbind(1, C) %*% A_model$coefficients %>% expit
    N1 <- ((A == 1) / ps)
    N0 <- ((A == 0) / (1 - ps))
    D1 <- at_risk_tb %*% N1 %>% as.numeric / n
    D0 <- at_risk_tb %*% N0 %>% as.numeric / n
    nan_to_zero <- function(x) {x[is.nan(x)] <- 0 ; return(x)}
    if(is.null(B_points)) {
      return(
        list(
          estimate = data.frame(
            time = baseline$time,
            NDE = 1/n * (
              nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% (m_ratio * N1) / D1) -
                nan_to_zero(((dn - dn_pred0) * at_risk_tb) %*% N0 / D0) +
                nan_to_zero(((dn_pred1 - g_com(1, 0) - dn_pred0 + g_com(0, 0)) * at_risk_tb) %*% N0 / D0) +
                nan_to_zero(apply(g_com(1, 0) - g_com(0, 0), 1, sum))
            ) %>% as.numeric %>% cumsum,
            NIE = 1/n * (
              nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% ((1 - m_ratio) * N1) / D1) +
                nan_to_zero(((dn_pred1 - g_com(1, 1)) * at_risk_tb) %*% N1 / D1) -
                nan_to_zero(((dn_pred1 - g_com(1, 0)) * at_risk_tb) %*% N0 / D0) +
                nan_to_zero(apply(g_com(1, 1) - g_com(1, 0), 1, sum))
            ) %>% as.numeric %>% cumsum
          ),
          M_model = M_model,
          Y_model = Y_model
        )
        
      )
    } else {
      nde <- 1/n * (
        nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% (m_ratio * N1) / D1) -
          nan_to_zero(((dn - dn_pred0) * at_risk_tb) %*% N0 / D0) +
          nan_to_zero(((dn_pred1 - g_com(1, 0) - dn_pred0 + g_com(0, 0)) * at_risk_tb) %*% N0 / D0) +
          nan_to_zero(apply(g_com(1, 0) - g_com(0, 0), 1, sum))
      ) %>% as.numeric %>% cumsum
      nie <- 1/n * (
        nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% ((1 - m_ratio) * N1) / D1) +
          nan_to_zero(((dn_pred1 - g_com(1, 1)) * at_risk_tb) %*% N1 / D1) -
          nan_to_zero(((dn_pred1 - g_com(1, 0)) * at_risk_tb) %*% N0 / D0) +
          nan_to_zero(apply(g_com(1, 1) - g_com(1, 0), 1, sum))
      ) %>% as.numeric %>% cumsum
      return(
        c(stepfun(baseline$time,
                  c(0, nde))(B_points),
          stepfun(baseline$time,
                  c(0, nie))(B_points))
      )
    }
  }
}

