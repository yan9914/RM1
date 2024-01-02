#' Mediation Analysis With Recurrent event outcome
#'
#' Perform causal mediation analysis with recurrent event outcomes, using either binary or continuous mediator. The distribution of the mediator is estimated through logistic or linear regression models. The distribution of the exposure and outcome are estimated by logistic regression and proportional mean model, respectively. This estimation method requires only two of the three models - exposure, mediator, and outcome - to be correctly specified in order to obtain an unbiased estimate.
#'
#' @param data a data.frame, containing the variables in the model. Specifically, column names \code{id}, \code{t.start}, \code{t.stop}, and \code{event} are required.
#' @param exposure a character value specifying the column name of exposure.
#' @param mediator a character value specifying the column name of mediator.
#' @param pretreat_confounder character vector specifying the column names of pre-treatment confounders.
#' @param posttreat_confounder character vector specifying the column names of post-treatment confounders.
#' @param mediator_dist character. If \code{normal}, linear regression is used to predict mediator. If \code{bernoulli}, logistic regression is used to predict mediator.
#' @param B_number integer; the times of bootstrapping.
#' @param B_points a numeric vector; the time points that used to compute pointwise 95 percentage confidence intervals by bootstrapping.
#' @param B_seed integer. \code{set.seed(B_seed)}
#' @param g_computation a optional logical value. The default is \code{FALSE}. If \code{TRUE}, only the model of exposure and mediator will be used, but then there is no triple robustness.
#' @param estimation logical.
#'
#'
#' @export
#' @seealso ...
#'
#' @examples
#' \dontrun{
#' RecMed(data = sim_data1, exposure = "A", mediator = "M", pretreat_confounder = c("C1", "C2"),
#'  B_number = 1000 , B_points = 24 * c(0.2,0.4,0.5,0.6,0.8))
#' }
#'

RecMed <- function(data, exposure, mediator, pretreat_confounder = NULL, posttreat_confounder = NULL, mediator_dist = "normal", B_number, B_points, B_seed = 123, g_computation = FALSE, estimation = TRUE) {
  expit <- function(x) exp(x)/(1+exp(x))
  f <- function(data, exposure, mediator, pretreat_confounder = NULL, posttreat_confounder = NULL, mediator_dist = "normal", g_computation = FALSE, B_points = NULL) {
    setDT(data)
    confounders <- union(pretreat_confounder, posttreat_confounder)
    subdata <- data[data[, .I[1], by = id]$V1, mget(c("id", exposure, mediator, confounders))]
    n <- nrow(subdata)
    A <- subdata[, get(exposure)] %>% as.matrix
    M <- subdata[, get(mediator)] %>% as.matrix
    C <- subdata[, mget(confounders)] %>% as.matrix
    C_pre <- subdata[, mget(pretreat_confounder)] %>% as.matrix
    Y_model <- formula(paste0("Surv(time = t.start, time2 = t.stop, event = event)~",
                              paste0(c(exposure, mediator, confounders, "cluster(id)"), collapse = "+"))) %>%
      coxph(data = data, control = coxph.control(timefix = FALSE))
    Y_coef_A <- Y_model$coefficients[exposure]
    Y_coef_M <- Y_model$coefficients[mediator]
    Y_coef_C <- Y_model$coefficients[confounders]
    baseline <- basehaz(Y_model, centered = F)
    n_time <- length(baseline$time)
    baseline_time <- baseline$time
    baseline_hazard <- baseline$hazard
    if(mediator_dist == "normal") {
      M_model <- lm(formula(paste0(mediator, "~", paste(c(exposure, confounders), collapse = "+"))), data = subdata)
      M_coef_int <- M_model$coefficients["(Intercept)"]
      M_coef_A <- M_model$coefficients[exposure]
      M_coef_C <- M_model$coefficients[confounders]
      M_sd <- summary(M_model)$sigma
      if(g_computation) {
        eta_hat <- baseline_hazard * exp(M_coef_int * Y_coef_M + 1/2 * Y_coef_M^2 * M_sd^2) *
          mean(exp(C %*% M_coef_C * Y_coef_M + C %*% Y_coef_C))
        if(is.null(B_points)) {
          return(
            list(
              estimate = data.frame(
                time = baseline_time,
                NDE = eta_hat * (exp(Y_coef_A)-1),
                NIE = eta_hat * (exp(Y_coef_A + Y_coef_M * M_coef_A) - exp(Y_coef_A))
              )
            )
          )
        } else {
          nde <- eta_hat * (exp(Y_coef_A)-1)
          nie <- eta_hat * (exp(Y_coef_A + Y_coef_M * M_coef_A) - exp(Y_coef_A))
          return(
            c(stepfun(baseline_time,
                      c(0, nde))(B_points),
              stepfun(baseline_time,
                      c(0, nie))(B_points))
          )
        }
      }
      g00 <- (diff.default(c(0, baseline_hazard)) * exp(M_coef_int * Y_coef_M + 1/2 * Y_coef_M^2 * M_sd^2)) %*%
        t(exp(C %*% M_coef_C * Y_coef_M + C %*% Y_coef_C))
      g11 <- g00 * exp(Y_coef_A + Y_coef_M * M_coef_A)
      g10 <- g00 * exp(Y_coef_A)
      A_model <- glm(formula(paste0(exposure, "~", paste(pretreat_confounder, collapse = "+"))), data = subdata, family = "binomial")
      if(identical(data$t.stop %>% unique %>% length, baseline_time)) stop("time points are not all equal to package(coxph) output")
      at_risk_tb <- data[,(function(x) {tmp <- rep(1L, n_time) ; tmp[baseline_time >= x[[1]][x[[2]] == 0]] <- 0L ; return(tmp)})(list(t.stop, event)), by = id] %>%
        .[, V1] %>% matrix(nrow = n_time)
      dn <- data[,(function(x) {tmp <- rep(0L, n_time) ; tmp[baseline_time %in% x[[1]][x[[2]] == 1]] <- 1L ; return(tmp)})(list(t.stop, event)), by = id] %>%
        .[, V1] %>% matrix(nrow = n_time)
      dn_pred1 <- diff.default(c(0, baseline_hazard)) %*% t(exp(cbind(1, M, C) %*% Y_model$coefficients))
      dn_pred0 <- diff.default(c(0, baseline_hazard)) %*% t(exp(cbind(0, M, C) %*% Y_model$coefficients))
      m_ratio <- exp(-1/2/M_sd^2 * ((M - cbind(1, 0, C) %*% M_model$coefficients)^2 - (M - cbind(1, 1, C) %*% M_model$coefficients)^2))
      ps <- cbind(1, C_pre) %*% A_model$coefficients %>% expit
      N1 <- ((A == 1) / ps)
      N0 <- ((A == 0) / (1 - ps))
      D1 <- at_risk_tb %*% N1 %>% as.numeric / n
      D0 <- at_risk_tb %*% N0 %>% as.numeric / n
      nan_to_zero <- function(x) {x[is.nan(x)] <- 0 ; return(x)}
      if(is.null(B_points)) {
        return(
          list(
            estimate = data.frame(
              time = baseline_time,
              NDE = 1/n * (
                nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% (m_ratio * N1) / D1) -
                  nan_to_zero(((dn - dn_pred0) * at_risk_tb) %*% N0 / D0) +
                  nan_to_zero(((dn_pred1 - g10 - dn_pred0 + g00) * at_risk_tb) %*% N0 / D0) +
                  nan_to_zero(apply(g10 - g00, 1, sum))
              ) %>% as.numeric %>% cumsum,
              NIE = 1/n * (
                nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% ((1 - m_ratio) * N1) / D1) +
                  nan_to_zero(((dn_pred1 - g11) * at_risk_tb) %*% N1 / D1) -
                  nan_to_zero(((dn_pred1 - g10) * at_risk_tb) %*% N0 / D0) +
                  nan_to_zero(apply(g11 - g10, 1, sum))
              ) %>% as.numeric %>% cumsum
            ),
            A_model = A_model,
            M_model = M_model,
            Y_model = Y_model
          )
        )
      } else {
        nde <- 1/n * (
          nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% (m_ratio * N1) / D1) -
            nan_to_zero(((dn - dn_pred0) * at_risk_tb) %*% N0 / D0) +
            nan_to_zero(((dn_pred1 - g10 - dn_pred0 + g00) * at_risk_tb) %*% N0 / D0) +
            nan_to_zero(apply(g10 - g00, 1, sum))
        ) %>% as.numeric %>% cumsum
        nie <- 1/n * (
          nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% ((1 - m_ratio) * N1) / D1) +
            nan_to_zero(((dn_pred1 - g11) * at_risk_tb) %*% N1 / D1) -
            nan_to_zero(((dn_pred1 - g10) * at_risk_tb) %*% N0 / D0) +
            nan_to_zero(apply(g11 - g10, 1, sum))
        ) %>% as.numeric %>% cumsum
        return(
          c(stepfun(baseline_time,
                    c(0, nde))(B_points),
            stepfun(baseline_time,
                    c(0, nie))(B_points))
        )
      }
    } else if(mediator_dist == "bernoulli"){
      M_model <- glm(reformulate(paste(c(exposure, confounders), collapse = "+"), mediator), data = subdata, family = "binomial")
      # M_coef_int <- M_model$coefficients["(Intercept)"]
      # M_coef_A <- M_model$coefficients[exposure]
      # M_coef_C <- M_model$coefficients[confounders]
      prob1 <- expit(cbind(1, 1, C) %*% M_model$coefficients)
      prob0 <- expit(cbind(1, 0, C) %*% M_model$coefficients)
      if(g_computation) {
        tmp11 <- exp(Y_coef_A) * mean((1 - prob1 + prob1 * exp(Y_coef_M)) * exp(C %*% Y_coef_C))
        tmp10 <- exp(Y_coef_A) * mean((1 - prob0 + prob0 * exp(Y_coef_M)) * exp(C %*% Y_coef_C))
        tmp00 <- mean((1 - prob0 + prob0 * exp(Y_coef_M)) * exp(C %*% Y_coef_C))
        if(is.null(B_points)) {
          return(
            list(
              estimate = data.frame(
                time = baseline_time,
                NDE = baseline_hazard * (tmp10 - tmp00),
                NIE = baseline_hazard * (tmp11 - tmp10)
              )
            )
          )
        } else {
          nde <- baseline_hazard * (tmp10 - tmp00)
          nie <- baseline_hazard * (tmp11 - tmp10)
          return(
            c(stepfun(baseline_time,
                      c(0, nde))(B_points),
              stepfun(baseline_time,
                      c(0, nie))(B_points))
          )
        }
      }
      g00 <- diff.default(c(0, baseline_hazard)) %*% t((1 - prob0 + prob0 * exp(Y_coef_M) * exp(C %*% Y_coef_C)))
      g11 <- diff.default(c(0, baseline_hazard)) %*% t(exp(Y_coef_A) * (1 - prob1 + prob1 * exp(Y_coef_M)) * exp(C %*% Y_coef_C))
      g10 <- g00 * exp(Y_coef_A)
      g00 %<>% as.matrix()
      g11 %<>% as.matrix()
      g10 %<>% as.matrix()
      A_model <- glm(formula(paste0(exposure, "~", paste(pretreat_confounder, collapse = "+"))), data = subdata, family = "binomial")
      if(identical(data$t.stop %>% unique %>% length, baseline_time)) stop("time points are not all equal to package(coxph) output")
      at_risk_tb <- data[,(function(x) {tmp <- rep(1L, n_time) ; tmp[baseline_time >= x[[1]][x[[2]] == 0]] <- 0L ; return(tmp)})(list(t.stop, event)), by = id] %>%
        .[, V1] %>% matrix(nrow = n_time)
      dn <- data[,(function(x) {tmp <- rep(0L, n_time) ; tmp[baseline_time %in% x[[1]][x[[2]] == 1]] <- 1L ; return(tmp)})(list(t.stop, event)), by = id] %>%
        .[, V1] %>% matrix(nrow = n_time)
      dn_pred1 <- diff.default(c(0, baseline_hazard)) %*% t(exp(cbind(1, M, C) %*% Y_model$coefficients))
      dn_pred0 <- diff.default(c(0, baseline_hazard)) %*% t(exp(cbind(0, M, C) %*% Y_model$coefficients))
      m_ratio <- (M * prob0 + (1 - M) * (1 - prob0)) / (M * prob1 + (1 - M) * (1 - prob1))
      ps <- cbind(1, C_pre) %*% A_model$coefficients %>% expit
      N1 <- ((A == 1) / ps)
      N0 <- ((A == 0) / (1 - ps))
      D1 <- at_risk_tb %*% N1 %>% as.numeric / n
      D0 <- at_risk_tb %*% N0 %>% as.numeric / n
      nan_to_zero <- function(x) {x[is.nan(x)] <- 0 ; return(x)}
      if(is.null(B_points)) {
        return(
          list(
            estimate = data.frame(
              time = baseline_time,
              NDE = 1/n * (
                nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% (m_ratio * N1) / D1) -
                  nan_to_zero(((dn - dn_pred0) * at_risk_tb) %*% N0 / D0) +
                  nan_to_zero(((dn_pred1 - g10 - dn_pred0 + g00) * at_risk_tb) %*% N0 / D0) +
                  nan_to_zero(apply(g10 - g00, 1, sum))
              ) %>% as.numeric %>% cumsum,
              NIE = 1/n * (
                nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% ((1 - m_ratio) * N1) / D1) +
                  nan_to_zero(((dn_pred1 - g11) * at_risk_tb) %*% N1 / D1) -
                  nan_to_zero(((dn_pred1 - g10) * at_risk_tb) %*% N0 / D0) +
                  nan_to_zero(apply(g11 - g10, 1, sum))
              ) %>% as.numeric %>% cumsum
            ),
            A_model = A_model,
            M_model = M_model,
            Y_model = Y_model
          )
        )
      } else {
        nde <- 1/n * (
          nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% (m_ratio * N1) / D1) -
            nan_to_zero(((dn - dn_pred0) * at_risk_tb) %*% N0 / D0) +
            nan_to_zero(((dn_pred1 - g10 - dn_pred0 + g00) * at_risk_tb) %*% N0 / D0) +
            nan_to_zero(apply(g10 - g00, 1, sum))
        ) %>% as.numeric %>% cumsum
        nie <- 1/n * (
          nan_to_zero(((dn - dn_pred1) * at_risk_tb) %*% ((1 - m_ratio) * N1) / D1) +
            nan_to_zero(((dn_pred1 - g11) * at_risk_tb) %*% N1 / D1) -
            nan_to_zero(((dn_pred1 - g10) * at_risk_tb) %*% N0 / D0) +
            nan_to_zero(apply(g11 - g10, 1, sum))
        ) %>% as.numeric %>% cumsum
        return(
          c(stepfun(baseline_time,
                    c(0, nde))(B_points),
            stepfun(baseline_time,
                    c(0, nie))(B_points))
        )
      }
    }
  }
  set.seed(B_seed)
  setDT(data)
  n <- data$id %>% unique %>% length
  mx <- matrix(NA, nrow = 2 * length(B_points), ncol = B_number)
  for(i in 1:B_number) {
    mx[,i] <- f(data = data[data.table(id = unique(data$id)[sample.int(n, n, replace = TRUE)], id2 = 1:n), on = "id"][, id := id2][, id2 := NULL],
                exposure, mediator, pretreat_confounder, posttreat_confounder, mediator_dist, g_computation, B_points = B_points)
  }
  TE_mx <- mx[1:length(B_points),] + mx[(1:length(B_points))+length(B_points),]
  result_lower <- apply(mx, 1, quantile, probs = 0.025)
  result_upper <- apply(mx, 1, quantile, probs = 0.975)
  result_sd <- apply(mx, 1, sd)
  TE_lower1 <- apply(TE_mx, 1, quantile, probs = 0.025)
  TE_upper1 <- apply(TE_mx, 1, quantile, probs = 0.975)
  TE_sd1 <- apply(TE_mx, 1, sd)
  if(estimation) {
    tmp <- f(data, exposure, mediator, pretreat_confounder, posttreat_confounder, mediator_dist, g_computation, B_points = NULL)
    return(
      list(
        estimate = tmp$estimate,
        NDE_lower = result_lower[1:length(B_points)],
        NDE_upper = result_upper[1:length(B_points)],
        NIE_lower = result_lower[(1:length(B_points))+length(B_points)],
        NIE_upper = result_upper[(1:length(B_points))+length(B_points)],
        TE_lower = TE_lower1,
        TE_upper = TE_upper1,
        NDE_sd = result_sd[1:length(B_points)],
        NIE_sd = result_sd[(1:length(B_points))+length(B_points)],
        TE_sd = TE_sd1,
        A_model = tmp$A_model,
        M_model = tmp$M_model,
        Y_model = tmp$Y_model
      )
    )
  } else {
    return(
      list(
        NDE_lower = result_lower[1:length(B_points)],
        NDE_upper = result_upper[1:length(B_points)],
        NIE_lower = result_lower[(1:length(B_points))+length(B_points)],
        NIE_upper = result_upper[(1:length(B_points))+length(B_points)],
        TE_lower = TE_lower1,
        TE_upper = TE_upper1,
        NDE_sd = result_sd[1:length(B_points)],
        NIE_sd = result_sd[(1:length(B_points))+length(B_points)],
        TE_sd = TE_sd1
      )
    )
  }
}
