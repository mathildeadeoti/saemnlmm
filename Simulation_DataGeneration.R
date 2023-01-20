# This code generates and stores data #

set.seed(2022)
library(nlmm)
# Data generation ####
M <- 250 # data sets replications
ni = 5 # measurement replication
fixef = c(200, 700, 350)
sigma2 = .5 # residual scale

#complexity <- c("RIS", "RI")
spSz <- c(100, 500, 1000)  # sample size

folder <- "./Data"# Data folder direction
# RIS ####
compl <- "RIS"
B_i = matrix(c(1,0,0,0,1,0,0,0,1,0,1,0,0,0,1), nrow = 3)
q = ncol(B_i); q
Sigma = diag(q)
Wposition = matrix(c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,
                     FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE), nrow = 3)

rmydist <- function(n) runif(n, min = 0, max = 2000)
rZ1 <- function(n) rbinom(n = n, size = 1, prob = .2)
rZ2 <- function(n) rbinom(n = n, size = 1, prob = .6)
#rZ1(n = 10)

## Normal (symetric) #####
distr <- "Ind"
lambda = 0
lambdab = rep(0, 5)

for (n in spSz) {
  FkData <- lapply(1:M, function(d){
    nlssnmdata(n = n, n_i = ni, n_x = 1, n_z = 0, n_w = 2,
               nlfun = function(phi, X) phi[1]/(1 + exp((phi[2] - X)/phi[3])),
               Xdist = "mydist", A_i = diag(3), #Zdist = rep("unif", n_z), Zposition = NULL,
               B_i = B_i, Wdist = c("Z1", "Z2"), Wposition = Wposition,
               fixef = fixef, omegabar2 = MOmegabar(omega = sigma2, lambda = lambda),
               deltao = Mdelta(omega = sigma2, lambda = lambda),
               Omegabar_b = MOmegabar(omega = Sigma, lambda = lambdab),
               delta_b = Mdelta(omega = Sigma, lambda = lambdab),
               smixvar = distr, smixvar_b = distr)
  })
  nm <- paste(distr, n, compl, ".RData", sep = "")
  save(FkData, file = paste(folder, nm, sep = "/"))
  #load(paste(folder, nm, sep = "/"))
}

## Skew t ####
distr <- "Gamma"
lambda = 5
lambdab = rep(5, 5)

for (n in spSz) {
  FkData <- lapply(1:M, function(d){
    nlssnmdata(n = n, n_i = ni, n_x = 1, n_z = 0, n_w = 2,
               nlfun = function(phi, X) phi[1]/(1 + exp((phi[2] - X)/phi[3])),
               Xdist = "mydist", A_i = diag(3), #Zdist = rep("unif", n_z), Zposition = NULL,
               B_i = B_i, Wdist = c("Z1", "Z2"), Wposition = Wposition,
               fixef = fixef, omegabar2 = nlmm::MOmegabar(omega = sigma2, lambda = lambda),
               deltao = nlmm::Mdelta(omega = sigma2, lambda = lambda),
               Omegabar_b = nlmm::MOmegabar(omega = Sigma, lambda = lambdab),
               delta_b = nlmm::Mdelta(omega = Sigma, lambda = lambdab),
               smixvar = distr, smixvar_b = distr, nu = 5, nu_b = 5)
  })
  nm <- paste(distr, n, compl, ".RData", sep = "")
  save(FkData, file = paste(folder, nm, sep = "/"))
  #load(paste(folder, nm, sep = "/"))
}

# RI ####
compl <- "RI"
q <- length(fixef); q
B_i = diag(q)
#q = ncol(B_i); q
Sigma = diag(q)

rmydist <- function(n) runif(n, min = 0, max = 2000)
#rZ1 <- function(n) rbinom(n = n, size = 1, prob = .2)
#rZ2 <- functi#on(n) rbinom(n = n, size = 1, prob = .6)
#rZ1(n = 10)

## Normal (symetric) #####
distr <- "Ind"
lambda = 0
lambdab = rep(0, q)

for (n in spSz) {
  FkData <- lapply(1:M, function(d){
    nlssnmdata(n = n, n_i = ni, n_x = 1, n_z = 0, n_w = 0,
               nlfun = function(phi, X) phi[1]/(1 + exp((phi[2] - X)/phi[3])),
               Xdist = "mydist", A_i = diag(3), #Zdist = rep("unif", n_z), Zposition = NULL,
               B_i = B_i, #Wdist = c("Z1", "Z2"), Wposition = Wposition,
               fixef = fixef, omegabar2 = MOmegabar(omega = sigma2, lambda = lambda),
               deltao = Mdelta(omega = sigma2, lambda = lambda),
               Omegabar_b = MOmegabar(omega = Sigma, lambda = lambdab),
               delta_b = Mdelta(omega = Sigma, lambda = lambdab),
               smixvar = distr, smixvar_b = distr)
  })
  nm <- paste(distr, n, compl, ".RData", sep = "")
  save(FkData, file = paste(folder, nm, sep = "/"))
}

## Skew t ####
distr <- "Gamma"
lambda = 5
lambdab = rep(5, q)

for (n in spSz) {
  FkData <- lapply(1:M, function(d){
    nlssnmdata(n = n, n_i = ni, n_x = 1, n_z = 0, n_w = 0,
               nlfun = function(phi, X) phi[1]/(1 + exp((phi[2] - X)/phi[3])),
               Xdist = "mydist", A_i = diag(3), #Zdist = rep("unif", n_z), Zposition = NULL,
               B_i = B_i, #Wdist = c("Z1", "Z2"), Wposition = Wposition,
               fixef = fixef, omegabar2 = MOmegabar(omega = sigma2, lambda = lambda),
               deltao = Mdelta(omega = sigma2, lambda = lambda),
               Omegabar_b = MOmegabar(omega = Sigma, lambda = lambdab),
               delta_b = Mdelta(omega = Sigma, lambda = lambdab),
               smixvar = distr, smixvar_b = distr, nu = 5, nu_b = 5)
  })
  nm <- paste(distr, n, compl, ".RData", sep = "")
  save(FkData, file = paste(folder, nm, sep = "/"))
}
