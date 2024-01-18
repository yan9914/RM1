source("fun_triply_robust.R")
library(parallel)
library(ggplot2)

true_value <- readRDS("true_value_20221011.rds")
true_value_0.5 <- readRDS("true_value_20221021_0.5.rds")
true_value_2 <- readRDS("true_value_20221021_2.rds")
true_value_3 <- readRDS("true_value_20221021_3.rds")

true_value_NDE <- stepfun(true_value$time, c(0, true_value$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8))
true_value_NIE <- stepfun(true_value$time, c(0, true_value$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8))

true_value_NDE_0.5 <- stepfun(true_value_0.5$time, c(0, true_value_0.5$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8))
true_value_NIE_0.5 <- stepfun(true_value_0.5$time, c(0, true_value_0.5$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8))

true_value_NDE_2 <- stepfun(true_value_2$time, c(0, true_value_2$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8))
true_value_NIE_2 <- stepfun(true_value_2$time, c(0, true_value_2$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8))

true_value_NDE_3 <- stepfun(true_value_3$time, c(0, true_value_3$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8))
true_value_NIE_3 <- stepfun(true_value_3$time, c(0, true_value_3$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8))

# eta = 0.5

simf <- function(i) {
  set.seed(i)
  data <- data_generation(n = 2e2,
                          a_par = c(1, 1, -2),
                          m_par = c(3, -1, -1, 1.5),
                          m_sd = 2,
                          y_par = c(1, 0.1, 0.2, -0.22),
                          cum_baseline = function(x) 0.5 * 0.05 * x,
                          tau = 24)
  TR <- f(data, "A", "M", c("C1", "C2"), g_computation = F)$estimate
  GC <- f(data, "A", "M", c("C1", "C2"), g_computation = T)$estimate
  TR_NDE <- stepfun(TR$time, c(0, TR$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NDE_0.5
  TR_NIE <- stepfun(TR$time, c(0, TR$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NIE_0.5
  GC_NDE <- stepfun(GC$time, c(0, GC$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NDE_0.5
  GC_NIE <- stepfun(GC$time, c(0, GC$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NIE_0.5
  cbind(c("TR", "TR", "GC", "GC"), c("NDE", "NIE", "NDE", "NIE"), rbind(TR_NDE, TR_NIE, GC_NDE, GC_NIE))
}

a <- proc.time()
cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores)
clusterExport(cl, c("true_value_NDE_0.5",
                    "true_value_NIE_0.5",
                    "simf",
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
result1 <- parLapply(cl, 1:1000, function(x) simf(x))
stopCluster(cl)
proc.time() - a

# eta = 1
simf <- function(i) {
  set.seed(i)
  data <- data_generation(n = 2e2,
                          a_par = c(1, 1, -2),
                          m_par = c(3, -1, -1, 1.5),
                          m_sd = 2,
                          y_par = c(1, 0.1, 0.2, -0.22),
                          cum_baseline = function(x) 0.05 * x,
                          tau = 24)
  TR <- f(data, "A", "M", c("C1", "C2"), g_computation = F)$estimate
  GC <- f(data, "A", "M", c("C1", "C2"), g_computation = T)$estimate
  TR_NDE <- stepfun(TR$time, c(0, TR$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NDE
  TR_NIE <- stepfun(TR$time, c(0, TR$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NIE
  GC_NDE <- stepfun(GC$time, c(0, GC$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NDE
  GC_NIE <- stepfun(GC$time, c(0, GC$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NIE
  cbind(c("TR", "TR", "GC", "GC"), c("NDE", "NIE", "NDE", "NIE"), rbind(TR_NDE, TR_NIE, GC_NDE, GC_NIE))
}

a <- proc.time()
cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores)
clusterExport(cl, c("true_value_NDE",
                    "true_value_NIE",
                    "simf",
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
result2 <- parLapply(cl, 1:1000, function(x) simf(x))
stopCluster(cl)
proc.time() - a

# eta = 2
simf <- function(i) {
  set.seed(i)
  data <- data_generation(n = 2e2,
                          a_par = c(1, 1, -2),
                          m_par = c(3, -1, -1, 1.5),
                          m_sd = 2,
                          y_par = c(1, 0.1, 0.2, -0.22),
                          cum_baseline = function(x) 2 * 0.05 * x,
                          tau = 24)
  TR <- f(data, "A", "M", c("C1", "C2"), g_computation = F)$estimate
  GC <- f(data, "A", "M", c("C1", "C2"), g_computation = T)$estimate
  TR_NDE <- stepfun(TR$time, c(0, TR$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NDE_2
  TR_NIE <- stepfun(TR$time, c(0, TR$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NIE_2
  GC_NDE <- stepfun(GC$time, c(0, GC$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NDE_2
  GC_NIE <- stepfun(GC$time, c(0, GC$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NIE_2
  cbind(c("TR", "TR", "GC", "GC"), c("NDE", "NIE", "NDE", "NIE"), rbind(TR_NDE, TR_NIE, GC_NDE, GC_NIE))
}

a <- proc.time()
cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores)
clusterExport(cl, c("true_value_NDE_2",
                    "true_value_NIE_2",
                    "simf",
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
result3 <- parLapply(cl, 1:1000, function(x) simf(x))
stopCluster(cl)
proc.time() - a

# eta = 3
simf <- function(i) {
  set.seed(i)
  data <- data_generation(n = 2e2,
                          a_par = c(1, 1, -2),
                          m_par = c(3, -1, -1, 1.5),
                          m_sd = 2,
                          y_par = c(1, 0.1, 0.2, -0.22),
                          cum_baseline = function(x) 3 * 0.05 * x,
                          tau = 24)
  TR <- f(data, "A", "M", c("C1", "C2"), g_computation = F)$estimate
  GC <- f(data, "A", "M", c("C1", "C2"), g_computation = T)$estimate
  TR_NDE <- stepfun(TR$time, c(0, TR$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NDE_3
  TR_NIE <- stepfun(TR$time, c(0, TR$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NIE_3
  GC_NDE <- stepfun(GC$time, c(0, GC$NDE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NDE_3
  GC_NIE <- stepfun(GC$time, c(0, GC$NIE))(24 * c(0.2, 0.4, 0.5, 0.6, 0.8)) - true_value_NIE_3
  cbind(c("TR", "TR", "GC", "GC"), c("NDE", "NIE", "NDE", "NIE"), rbind(TR_NDE, TR_NIE, GC_NDE, GC_NIE))
}

a <- proc.time()
cpu.cores <- detectCores()
cl <- makeCluster(cpu.cores)
clusterExport(cl, c("true_value_NDE_3",
                    "true_value_NIE_3",
                    "simf",
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
result4 <- parLapply(cl, 1:1000, function(x) simf(x))
stopCluster(cl)
proc.time() - a

### table
result21 <- do.call("rbind", result1) %>% 
  as.data.table %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, .(mean(V3), mean(V4), mean(V5), mean(V6), mean(V7)), by = .(V1, V2)] %>%
  .[, eta := 0.5]
result22 <- do.call("rbind", result2) %>% 
  as.data.table %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, .(mean(V3), mean(V4), mean(V5), mean(V6), mean(V7)), by = .(V1, V2)] %>%
  .[, eta := 1]
result23 <- do.call("rbind", result3) %>% 
  as.data.table %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, .(mean(V3), mean(V4), mean(V5), mean(V6), mean(V7)), by = .(V1, V2)] %>%
  .[, eta := 2]
result24 <- do.call("rbind", result4) %>% 
  as.data.table %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, .(mean(V3), mean(V4), mean(V5), mean(V6), mean(V7)), by = .(V1, V2)] %>%
  .[, eta := 3]
result5 <- rbind(result21, result22, result23, result24)

result21 <- do.call("rbind", result1) %>% 
  as.data.table %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, .(sd(V3), sd(V4), sd(V5), sd(V6), sd(V7)), by = .(V1, V2)] %>%
  .[, eta := 0.5]
result22 <- do.call("rbind", result2) %>% 
  as.data.table %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, .(sd(V3), sd(V4), sd(V5), sd(V6), sd(V7)), by = .(V1, V2)] %>%
  .[, eta := 1]
result23 <- do.call("rbind", result3) %>% 
  as.data.table %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, .(sd(V3), sd(V4), sd(V5), sd(V6), sd(V7)), by = .(V1, V2)] %>%
  .[, eta := 2]
result24 <- do.call("rbind", result4) %>% 
  as.data.table %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>% 
  .[, .(sd(V3), sd(V4), sd(V5), sd(V6), sd(V7)), by = .(V1, V2)] %>%
  .[, eta := 3]
result5 <- rbind(result5, result21, result22, result23, result24)
result5 <- cbind(data.frame(type = rep(c("Mean", "SD"), each = 16)), result5)

# saveRDS(result5, "sim_freq_n200_20221114.rds")
###

result21 <- do.call("rbind", result1) %>% 
  as.data.table %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>%
  .[, eta := 0.5]
result22 <- do.call("rbind", result2) %>% 
  as.data.table %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>%
  .[, eta := 1]
result23 <- do.call("rbind", result3) %>% 
  as.data.table %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>%
  .[, eta := 2]
result24 <- do.call("rbind", result4) %>% 
  as.data.table %>%
  .[, c(.SD[,1:2],lapply(.SD[,3:7], as.numeric))] %>%
  .[, eta := 3]

result5 <- rbind(result21, result22, result23, result24)

theme_set(theme(panel.background = element_rect(fill = NA, colour = "black"),
                panel.grid = element_blank() ,
                legend.position = c(0.05,0.05)  ,
                legend.justification = c(0,0),
                legend.background = element_rect(fill = NA, colour = "black", size = 0.25)))

result5_TR <- melt(result5[V1 == "TR",], id.vars = c("eta", "V2"), measure.vars = c("V3", "V4", "V5", "V6", "V7"))
result5_GC <- melt(result5[V1 == "GC",], id.vars = c("eta", "V2"), measure.vars = c("V3", "V4", "V5", "V6", "V7"))

ggplot(result5_TR) +
  geom_boxplot(aes(factor(eta), value), outlier.size = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25) +
  labs(x = expression(zeta), y = "Bias") +
  ggh4x::facet_grid2(V2~variable, scales = "free",
                     labeller = labeller(.cols = as_labeller(c(V3 = "4.8 months (20%)",
                                                               V4 = "9.6 months (40%)",
                                                               V5 = "12 months (50%)",
                                                               V6 = "14.4 months (60%)",
                                                               V7 = "19.2 months (80%)"))))
ggplot(result5_GC) +
  geom_boxplot(aes(factor(eta), value), outlier.size = 0.75) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25) +
  labs(x = expression(zeta), y = "Bias") +
  ggh4x::facet_grid2(V2~variable, scales = "free",
                     labeller = labeller(.cols = as_labeller(c(V3 = "4.8 months (20%)",
                                                               V4 = "9.6 months (40%)",
                                                               V5 = "12 months (50%)",
                                                               V6 = "14.4 months (60%)",
                                                               V7 = "19.2 months (80%)"))))
