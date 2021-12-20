######### Tempered Expectation-Maximization algorithm #########
## Authors: Luca Brusa, Francesco Bartolucci, Fulvia Pennoni ##

files.sources = list.files(path = "./functions", pattern = ".R")
sapply(paste0("./functions/", files.sources), source)



#### 1.1. Latent class model - Simulated dataset ####
load("./Dataset/simulated_sample_LC.RData")

# Standard EM algorithm (computational time: around 1 second)
set.seed(24)
std_est <- est_LC(S = sample$S, yv = sample$yv, k = 4)

# T-EM algorithm with monotonic profile (computational time: around 1 second)
set.seed(24)
mon_est <- est_LC(S = sample$S, yv = sample$yv, k = 4,
                  algorithm = 2, profile_pars = list(alpha = 15, beta = 1))

# T-EM algorithm with oscillatinge profile (computational time: around 45 seconds)
set.seed(24)
osc_est <- est_LC(S = sample$S, yv = sample$yv, k = 4,
                  algorithm = 3, profile_pars = list(alpha = 0.3, beta = 100, ro = 5, T0 = 10))

# Results: maximized log-likelihood value with the three algorithms
print(c(std_est$lk, mon_est$lk, osc_est$lk))



#### 1.2. Latent class model - Anxiety and depression data ####
load("./Dataset/anxiety_depression.RData")
out <- MultiLCIRT::aggr_data(as.matrix(hads))

# Standard EM algorithm (computational time: around 1 second)
set.seed(24)
std_est <- est_LC(S = out$data_dis, yv = out$freq, k = 3)

# T-EM algorithm with monotonic profile (computational time: around 1 second)
set.seed(24)
mon_est <- est_LC(S = out$data_dis, yv = out$freq, k = 3,
                  algorithm = 2, profile_pars = list(alpha = 42, beta = 1.5))
# T-EM algorithm with oscillatinge profile (computational time: around 1 second)
set.seed(24)
osc_est <- est_LC(S = out$data_dis, yv = out$freq, k = 3,
                  algorithm = 3, profile_pars = list(alpha = 0.8, beta = 20, ro = 90, T0 = 10))

# Results: maximized log-likelihood value with the three algorithms
print(c(std_est$lk, mon_est$lk, osc_est$lk))



#### 2.1. Hidden Markov model (categorical responses) - Simulated dataset ####
load("./Dataset/simulated_sample_HM_cat.RData")

# Standard EM algorithm (computational time: around 5 seconds)
set.seed(24)
std_est <- est_LM(data = sample$Y, index = c("id", "time"), k = 3, modBasic = 1)

# T-EM algorithm with monotonic profile (computational time: around 15 seconds)
set.seed(24)
mon_est <- est_LM(data = sample$Y, index = c("id", "time"), k = 3, modBasic = 1,
                  algorithm = 2, profile_pars = list(alpha = 2, beta = 2.5))

# Results: maximized log-likelihood value with the two algorithms
print(c(std_est$lk, mon_est$lk))



#### 2.2. Hidden Markov model (categorical responses) - Criminal data ####
load("./Dataset/criminal.RData")

# Standard EM algorithm (computational time: around 10 seconds)
set.seed(24)
std_est <- est_LM(data = criminal, index = c("id", "time"), k = 4, modBasic = 0)

# T-EM algorithm with monotonic profile (computational time: around 100 seconds)
set.seed(24)
mon_est <- est_LM(data = criminal, index = c("id", "time"), k = 4, modBasic = 0,
                  algorithm = 2, profile_pars = list(alpha = 2, beta = 1.5))

# Results: maximized log-likelihood value with the two algorithms
print(c(std_est$lk, mon_est$lk))



#### 3. Hidden Markov model (continuous responses) - Simulated dataset ####
load("./Dataset/simulated_sample_HM_cont.RData")

# Standard EM algorithm (computational time: around 5 seconds)
set.seed(24)
std_est <- est_LM_cont(data = Y, index = c("id", "time"), k = 3, modBasic = 0)

# T-EM algorithm with oscillatinge profile (computational time: around 45 seconds)
set.seed(24)
osc_est <- est_LM_cont(data = Y, index = c("id", "time"), k = 3, modBasic = 0, 
                       algorithm = 3, profile_pars = list(alpha = 0.3, beta = 5, ro = 100, T0 = 2))

# Results: maximized log-likelihood value with the two algorithms
print(c(std_est$lk, osc_est$lk))