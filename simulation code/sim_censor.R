source("fun_triply_robust.R")
library(parallel)
library(ggplot2)

true_value <- readRDS("true_value_20221011.rds")
true_value_NDE <- stepfun(true_value$time, c(0, true_value$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8))
true_value_NIE <- stepfun(true_value$time, c(0, true_value$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8))

# PC = 0
simf <- function(i) {
  set.seed(i)
  data <- data_generation(n = 1e3,
                          a_par = c(1, 1, -2),
                          m_par = c(3, -1, -1, 1.5),
                          m_sd = 2,
                          y_par = c(1, 0.1, 0.2, -0.22),
                          cum_baseline = function(x) 0.05 * x,
                          tau = 24,
                          censoring = rep(25, 1e3))
  TR <- f_bootstrap(data, "A", "M", c("C1", "C2"), B_number = 500, B_points = 24 * c(0.2, 0.4, 0.5, 0.6, 0.8), g_computation = F)
  TR_NDE <- stepfun(TR$estimate$time, c(0, TR$estimate$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NDE
  TR_NIE <- stepfun(TR$estimate$time, c(0, TR$estimate$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NIE
  cbind(rep(c("Bias", "SD", "LB", "UB"), 2),
        rep(c("NDE", "NIE"), each = 4),
        rbind(TR_NDE, TR$NDE_sd, TR$NDE_lower, TR$NDE_upper,
              TR_NIE, TR$NIE_sd, TR$NIE_lower, TR$NIE_upper))
}

a <- proc.time()
cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores)
clusterExport(cl, c("true_value_NDE",
                    "true_value_NIE",
                    "simf",
                    "f",
                    "f_bootstrap",
                    "data_generation",
                    "expit",
                    "inv",
                    "simGSC"))
clusterEvalQ(cl, c(library(data.table),
                   library(magrittr),
                   library(survival),
                   library(reReg),
                   library(BB)))
result1 <- parLapply(cl, 1:1000, function(x) simf(x))
stopCluster(cl)
proc.time() - a

# PC = 15%
simf <- function(i) {
  set.seed(i)
  data <- data_generation(n = 1e3,
                          a_par = c(1, 1, -2),
                          m_par = c(3, -1, -1, 1.5),
                          m_sd = 2,
                          y_par = c(1, 0.1, 0.2, -0.22),
                          cum_baseline = function(x) 0.05 * x,
                          tau = 24,
                          censoring = runif(1e3, 0, 80))
  TR <- f_bootstrap(data, "A", "M", c("C1", "C2"), B_number = 500, B_points = 24 * c(0.2, 0.4, 0.5, 0.6, 0.8), g_computation = F)
  TR_NDE <- stepfun(TR$estimate$time, c(0, TR$estimate$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NDE
  TR_NIE <- stepfun(TR$estimate$time, c(0, TR$estimate$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NIE
  cbind(rep(c("Bias", "SD", "LB", "UB"), 2),
        rep(c("NDE", "NIE"), each = 4),
        rbind(TR_NDE, TR$NDE_sd, TR$NDE_lower, TR$NDE_upper,
              TR_NIE, TR$NIE_sd, TR$NIE_lower, TR$NIE_upper))
}

a <- proc.time()
cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores)
clusterExport(cl, c("true_value_NDE",
                    "true_value_NIE",
                    "simf",
                    "f",
                    "f_bootstrap",
                    "data_generation",
                    "expit",
                    "inv",
                    "simGSC"))
clusterEvalQ(cl, c(library(data.table),
                   library(magrittr),
                   library(survival),
                   library(reReg),
                   library(BB)))
result2 <- parLapply(cl, 1:1000, function(x) simf(x))
stopCluster(cl)
proc.time() - a

# PC = 30%
simf <- function(i) {
  set.seed(i)
  data <- data_generation(n = 1e3,
                          a_par = c(1, 1, -2),
                          m_par = c(3, -1, -1, 1.5),
                          m_sd = 2,
                          y_par = c(1, 0.1, 0.2, -0.22),
                          cum_baseline = function(x) 0.05 * x,
                          tau = 24,
                          censoring = runif(1e3, 0, 40))
  TR <- f_bootstrap(data, "A", "M", c("C1", "C2"), B_number = 500, B_points = 24 * c(0.2, 0.4, 0.5, 0.6, 0.8), g_computation = F)
  TR_NDE <- stepfun(TR$estimate$time, c(0, TR$estimate$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NDE
  TR_NIE <- stepfun(TR$estimate$time, c(0, TR$estimate$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NIE
  cbind(rep(c("Bias", "SD", "LB", "UB"), 2),
        rep(c("NDE", "NIE"), each = 4),
        rbind(TR_NDE, TR$NDE_sd, TR$NDE_lower, TR$NDE_upper,
              TR_NIE, TR$NIE_sd, TR$NIE_lower, TR$NIE_upper))
}

a <- proc.time()
cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores)
clusterExport(cl, c("true_value_NDE",
                    "true_value_NIE",
                    "simf",
                    "f",
                    "f_bootstrap",
                    "data_generation",
                    "expit",
                    "inv",
                    "simGSC"))
clusterEvalQ(cl, c(library(data.table),
                   library(magrittr),
                   library(survival),
                   library(reReg),
                   library(BB)))
result3 <- parLapply(cl, 1:1000, function(x) simf(x))
stopCluster(cl)
proc.time() - a

# PC = 50%
simf <- function(i) {
  set.seed(i)
  data <- data_generation(n = 1e3,
                          a_par = c(1, 1, -2),
                          m_par = c(3, -1, -1, 1.5),
                          m_sd = 2,
                          y_par = c(1, 0.1, 0.2, -0.22),
                          cum_baseline = function(x) 0.05 * x,
                          tau = 24,
                          censoring = runif(1e3, 0, 24))
  TR <- f_bootstrap(data, "A", "M", c("C1", "C2"), B_number = 500, B_points = 24 * c(0.2, 0.4, 0.5, 0.6, 0.8), g_computation = F)
  TR_NDE <- stepfun(TR$estimate$time, c(0, TR$estimate$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NDE
  TR_NIE <- stepfun(TR$estimate$time, c(0, TR$estimate$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NIE
  cbind(rep(c("Bias", "SD", "LB", "UB"), 2),
        rep(c("NDE", "NIE"), each = 4),
        rbind(TR_NDE, TR$NDE_sd, TR$NDE_lower, TR$NDE_upper,
              TR_NIE, TR$NIE_sd, TR$NIE_lower, TR$NIE_upper))
}

a <- proc.time()
cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores)
clusterExport(cl, c("true_value_NDE",
                    "true_value_NIE",
                    "simf",
                    "f",
                    "f_bootstrap",
                    "data_generation",
                    "expit",
                    "inv",
                    "simGSC"))
clusterEvalQ(cl, c(library(data.table),
                   library(magrittr),
                   library(survival),
                   library(reReg),
                   library(BB)))
result4 <- parLapply(cl, 1:1000, function(x) simf(x))
stopCluster(cl)
proc.time() - a


result21 <- do.call("rbind", result1)
result22 <- do.call("rbind", result2)
result23 <- do.call("rbind", result3)
result24 <- do.call("rbind", result4)

# saveRDS(result21, "sim_censor_0_20221026.rds")
# saveRDS(result22, "sim_censor_15_20221026.rds")
# saveRDS(result23, "sim_censor_30_20221026.rds")
# saveRDS(result24, "sim_censor_50_20221026.rds")

as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V3))]
as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V4))]
as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V5))]
as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V6))]
as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V7))]

as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V3))]
as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V4))]
as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V5))]
as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V6))]
as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V7))]

as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V3))]
as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V4))]
as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V5))]
as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V6))]
as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V7))]

as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V3))]
as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V4))]
as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V5))]
as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V6))]
as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,mean(as.numeric(V7))]

## 

as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V3))]
as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V4))]
as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V5))]
as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V6))]
as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V7))]

as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V3))]
as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V4))]
as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V5))]
as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V6))]
as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V7))]

as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V3))]
as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V4))]
as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V5))]
as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V6))]
as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V7))]

as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V3))]
as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V4))]
as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V5))]
as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V6))]
as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NDE",] %>% .[,sd(as.numeric(V7))]

## NDE empirical standard error

as.data.table(result21) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V3))]
as.data.table(result21) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V4))]
as.data.table(result21) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V5))]
as.data.table(result21) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V6))]
as.data.table(result21) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V7))]

as.data.table(result22) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V3))]
as.data.table(result22) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V4))]
as.data.table(result22) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V5))]
as.data.table(result22) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V6))]
as.data.table(result22) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V7))]

as.data.table(result23) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V3))]
as.data.table(result23) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V4))]
as.data.table(result23) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V5))]
as.data.table(result23) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V6))]
as.data.table(result23) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V7))]

as.data.table(result24) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V3))]
as.data.table(result24) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V4))]
as.data.table(result24) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V5))]
as.data.table(result24) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V6))]
as.data.table(result24) %>% .[V1 == "SD" & V2 == "NDE",] %>% .[,mean(as.numeric(V7))]

## NDE coverage rate when PC = 0%

as.data.table(result21) %>% 
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V3") %>%
  .[, LB < true_value_NDE[1] & true_value_NDE[1] < UB] %>%
  mean
as.data.table(result21) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V4") %>%
  .[, LB < true_value_NDE[2] & true_value_NDE[2] < UB] %>%
  mean
as.data.table(result21) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V5") %>%
  .[, LB < true_value_NDE[3] & true_value_NDE[3] < UB] %>%
  mean
as.data.table(result21) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V6") %>%
  .[, LB < true_value_NDE[4] & true_value_NDE[4] < UB] %>%
  mean
as.data.table(result21) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V7") %>%
  .[, LB < true_value_NDE[5] & true_value_NDE[5] < UB] %>%
  mean

## NDE coverage rate when PC = 15%

as.data.table(result22) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V3") %>%
  .[, LB < true_value_NDE[1] & true_value_NDE[1] < UB] %>%
  mean
as.data.table(result22) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V4") %>%
  .[, LB < true_value_NDE[2] & true_value_NDE[2] < UB] %>%
  mean
as.data.table(result22) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V5") %>%
  .[, LB < true_value_NDE[3] & true_value_NDE[3] < UB] %>%
  mean
as.data.table(result22) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V6") %>%
  .[, LB < true_value_NDE[4] & true_value_NDE[4] < UB] %>%
  mean
as.data.table(result22) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V7") %>%
  .[, LB < true_value_NDE[5] & true_value_NDE[5] < UB] %>%
  mean

## NDE coverage rate when PC = 30%

as.data.table(result23) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V3") %>%
  .[, LB < true_value_NDE[1] & true_value_NDE[1] < UB] %>%
  mean
as.data.table(result23) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V4") %>%
  .[, LB < true_value_NDE[2] & true_value_NDE[2] < UB] %>%
  mean
as.data.table(result23) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V5") %>%
  .[, LB < true_value_NDE[3] & true_value_NDE[3] < UB] %>%
  mean
as.data.table(result23) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V6") %>%
  .[, LB < true_value_NDE[4] & true_value_NDE[4] < UB] %>%
  mean
as.data.table(result23) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V7") %>%
  .[, LB < true_value_NDE[5] & true_value_NDE[5] < UB] %>%
  mean

## NDE coverage rate when PC = 50%

as.data.table(result24) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V3") %>%
  .[, LB < true_value_NDE[1] & true_value_NDE[1] < UB] %>%
  mean
as.data.table(result24) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V4") %>%
  .[, LB < true_value_NDE[2] & true_value_NDE[2] < UB] %>%
  mean
as.data.table(result24) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V5") %>%
  .[, LB < true_value_NDE[3] & true_value_NDE[3] < UB] %>%
  mean
as.data.table(result24) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V6") %>%
  .[, LB < true_value_NDE[4] & true_value_NDE[4] < UB] %>%
  mean
as.data.table(result24) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NDE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V7") %>% 
  .[, LB < true_value_NDE[5] & true_value_NDE[5] < UB] %>%
  mean

# NIE

as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V3))]
as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V4))]
as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V5))]
as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V6))]
as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V7))]

as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V3))]
as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V4))]
as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V5))]
as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V6))]
as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V7))]

as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V3))]
as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V4))]
as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V5))]
as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V6))]
as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V7))]

as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V3))]
as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V4))]
as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V5))]
as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V6))]
as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,mean(as.numeric(V7))]

##
as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V3))]
as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V4))]
as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V5))]
as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V6))]
as.data.table(result21) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V7))]

as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V3))]
as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V4))]
as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V5))]
as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V6))]
as.data.table(result22) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V7))]

as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V3))]
as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V4))]
as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V5))]
as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V6))]
as.data.table(result23) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V7))]

as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V3))]
as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V4))]
as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V5))]
as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V6))]
as.data.table(result24) %>% .[V1 == "Bias" & V2 == "NIE",] %>% .[,sd(as.numeric(V7))]

##

as.data.table(result21) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V3))]
as.data.table(result21) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V4))]
as.data.table(result21) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V5))]
as.data.table(result21) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V6))]
as.data.table(result21) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V7))]

as.data.table(result22) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V3))]
as.data.table(result22) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V4))]
as.data.table(result22) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V5))]
as.data.table(result22) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V6))]
as.data.table(result22) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V7))]

as.data.table(result23) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V3))]
as.data.table(result23) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V4))]
as.data.table(result23) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V5))]
as.data.table(result23) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V6))]
as.data.table(result23) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V7))]

as.data.table(result24) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V3))]
as.data.table(result24) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V4))]
as.data.table(result24) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V5))]
as.data.table(result24) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V6))]
as.data.table(result24) %>% .[V1 == "SD" & V2 == "NIE",] %>% .[,mean(as.numeric(V7))]

##

as.data.table(result21) %>% 
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V3") %>%
  .[, LB < true_value_NIE[1] & true_value_NIE[1] < UB] %>% 
  mean
as.data.table(result21) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V4") %>%
  .[, LB < true_value_NIE[2] & true_value_NIE[2] < UB] %>%
  mean
as.data.table(result21) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V5") %>%
  .[, LB < true_value_NIE[3] & true_value_NIE[3] < UB] %>%
  mean
as.data.table(result21) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V6") %>%
  .[, LB < true_value_NIE[4] & true_value_NIE[4] < UB] %>%
  mean
as.data.table(result21) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V7") %>%
  .[, LB < true_value_NIE[5] & true_value_NIE[5] < UB] %>%
  mean

#

as.data.table(result22) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V3") %>%
  .[, LB < true_value_NIE[1] & true_value_NIE[1] < UB] %>%
  mean
as.data.table(result22) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V4") %>%
  .[, LB < true_value_NIE[2] & true_value_NIE[2] < UB] %>%
  mean
as.data.table(result22) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V5") %>%
  .[, LB < true_value_NIE[3] & true_value_NIE[3] < UB] %>%
  mean
as.data.table(result22) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V6") %>%
  .[, LB < true_value_NIE[4] & true_value_NIE[4] < UB] %>%
  mean
as.data.table(result22) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V7") %>%
  .[, LB < true_value_NIE[5] & true_value_NIE[5] < UB] %>%
  mean

#

as.data.table(result23) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V3") %>%
  .[, LB < true_value_NIE[1] & true_value_NIE[1] < UB] %>%
  mean
as.data.table(result23) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V4") %>%
  .[, LB < true_value_NIE[2] & true_value_NIE[2] < UB] %>%
  mean
as.data.table(result23) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V5") %>%
  .[, LB < true_value_NIE[3] & true_value_NIE[3] < UB] %>%
  mean
as.data.table(result23) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V6") %>%
  .[, LB < true_value_NIE[4] & true_value_NIE[4] < UB] %>%
  mean
as.data.table(result23) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V7") %>%
  .[, LB < true_value_NIE[5] & true_value_NIE[5] < UB] %>%
  mean
#
as.data.table(result24) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V3") %>%
  .[, LB < true_value_NIE[1] & true_value_NIE[1] < UB] %>%
  mean
as.data.table(result24) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V4") %>%
  .[, LB < true_value_NIE[2] & true_value_NIE[2] < UB] %>%
  mean
as.data.table(result24) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V5") %>%
  .[, LB < true_value_NIE[3] & true_value_NIE[3] < UB] %>%
  mean
as.data.table(result24) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V6") %>%
  .[, LB < true_value_NIE[4] & true_value_NIE[4] < UB] %>%
  mean
as.data.table(result24) %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, index := rep(1:1000, each = 8)] %>%
  .[V2 == "NIE",] %>%
  .[V1 == "LB" | V1 == "UB",] %>%
  dcast(index ~ V1, value.var = "V7") %>%
  .[, LB < true_value_NIE[5] & true_value_NIE[5] < UB] %>%
  mean

