# This code performs simulations, computes performance parameter and
#store data in files
library(nlmm)

#folder <- "D:/Consult/NonLinear/Data" # Data folder direction
folder <- "D:/LABEF/Confident/Prof. GLELE/For redaction/ASTB_manuscrit_accept?/For simulation/Data"## Result folder direction
#resFolder <- "D:/Consult/NonLinear/Results" ## Result folder direction
resFolder <- "D:/LABEF/Confident/Prof. GLELE/For redaction/ASTB_manuscrit_accept?/For simulation/Results"## Result folder direction

fixef = c(200, 700, 350)- 5
spSz <- c(100, 500, 1000) # n = 100, 500, 1000 # sample size
M = 250

# Skew Normal fits on all ####
fitdistr <- "Ind"
##  RIS data ####
compl <- "RIS"
SMSN <- c("Ind", "Gamma")
for (distr in SMSN) {
  cat("\n", "\nOOOO ~", distr, "data fits", "~ OOOO\n")
  for (n in spSz) {
    nm <- paste(distr, n, compl, ".RData", sep = "")
    load(paste(folder, nm, sep = "/"))
    ### SAEM Fits ####

    cat("\n ________ SAEM fitting on dataset____", nm, "\n")
    fakedata <- FkData[[1]]
    #head(fakedata$data)
    # Generate a framework for non linear mixed model fitting
    Nframe = nlmer.frame (y ~ phi[1]/(1 + exp((phi[2] - X)/phi[3])),
                          X.formula = ~ x1, name.mixef = c("phi1", "phi2", "phi3"),
                          fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                          fixed.matrix = fakedata$A, cluster = ~ cluster,
                          random = list(phi1 = ~ 1, phi2 = ~ w1, phi3 = ~ w2),
                          random.matrix = fakedata$B, data = fakedata$data)

    start = nlmer.start(frame = Nframe, smixvar = fitdistr, smixvar_b = fitdistr,
                        start = list(fixef = fixef))

    Desc <- paste(fitdistr, "fit on", nm)
    res <- vector(mode = "list", length = M)
    attr(res, "desc") <- Desc
    attr(res, "Algorithm") <- "SAEM"
    attr(res, "desc")
    attr(res, "Algorithm")

    for (d in 1:M) {
      cat("data", d, "\n")
      fakedata <- FkData[[d]]
      res[[d]] = ssmnnlmer(y ~ phi[1]/(1 + exp((phi[2] - X)/phi[3])),
                           X.formula = ~ x1, name.mixef = c("phi1", "phi2", "phi3"),
                           fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                           fixed.matrix = fakedata$A, cluster = ~ cluster,
                           random = list(phi1 = ~ 1, phi2 = ~ w1, phi3 = ~ w2),
                           random.matrix = fakedata$B, data = fakedata$data,
                           smixvar = fitdistr, smixvar_b = fitdistr,
                           start = start, control = list(maxwarmit = 300,
                                                         warmtol = 1e-1,
                                                         maxsmoothit = 200,
                                                         smoothtol = 1e-4,
                                                         m_b = 10,
                                                         m_mcem = 150,
                                                         m_saem = 15,
                                                         verbose = TRUE))
    }
    save(res, file = paste(resFolder, paste("SAEM", fitdistr, nm, sep = "_") ,sep = "/"))
    Theta <- makeTheta(res, algo = "SAEM")

    #TrueTheta <- Theta[,1] + .5 # REMOVE AND DEFINE IT ABOVE !!!
   # TrueTheta <-c(200,700,350,0.5735393, 0.5735393, 0.5735393,0,0,0,1,1,1,0.5735393,0.5)
    TrueTheta <-c(200,700,350, #fixef
                  0, 0, 0, 0, 0, # It was symmetric data
                  1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, # Omegabar
                  0, #delta_o
                  .5)

    per1 <- performance(result = res, Theta = Theta, TrueTheta = TrueTheta)
    per2 <- efficiency(result = res, algo = "SAEM")
    simResult <- list(Theta = Theta, per1 = per1, per2 = per2)
    save(simResult, file = paste(resFolder, paste("simRes", "SAEM", fitdistr, nm, sep = "_") ,sep = "/"))

    ## EM Fits ####
    cat("\n ________ EM__", fitdistr, "__fitting on dataset____", nm, "\n")
    Desc <- paste(fitdistr, "fit on", nm)
    EMres <- vector(mode = "list", length = M)
    attr(EMres, "desc") <- Desc
    attr(EMres, "Algorithm") <- "EM"
    attr(EMres, "desc")
    attr(EMres, "Algorithm")

    for (d in 1:M) {
      cat("data", d, "\n")
      fakedata <- FkData[[d]]
      EMres[[d]] = pereEM(maxit = 4000, tol = 1e-4,
                          formula = y ~ phi[1]/(1 + exp((phi[2] - X[1])/phi[3])),
                          X.formula = ~ x1, fixef.name = c("phi1", "phi2", "phi3"),
                          cluster = ~ cluster, family = fitdistr, data = fakedata$data,
                          start = list(alpha = fixef, sigma2 = .5, D = diag(5),
                                       lambda =rep(0, 5), nu = 5), #nu = c(.1, .5),
                          skew = TRUE, pverbose = TRUE, st.err = TRUE, nuEstim = FALSE,
                          errorF = "warning")
    }
    save(EMres, file = paste(resFolder, paste("EM", fitdistr, nm, sep = "_") ,sep = "/"))
    Theta <- makeTheta(EMres, algo = "EM")

    #TrueTheta <- Theta[,1] + .5 # REMOVE AND DEFINE IT ABOVE !!!
    TrueTheta <-c(200,700,350, #fixef
                  0, 0, 0, 0, 0, # It was symmetric data
                  .5,
                  1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1 # Omegabar
                  )


    per1 <- performance(result = EMres, Theta = Theta, TrueTheta = TrueTheta)
    per2 <- efficiency(result = EMres, algo = "EM")
    simResult <- list(Theta = Theta, per1 = per1, per2 = per2)
    save(simResult, file = paste(resFolder, paste("simRes", "EM", fitdistr, nm, sep = "_") ,sep = "/"))

  }
}


## RI data ####
compl <- "RI"
SMSN <- c("Ind", "Gamma")
for (distr in SMSN) {
  cat("\n", "\nOOOO ~", distr, "data fits", "~ OOOO\n")
  for (n in spSz) {
    nm <- paste(distr, n, compl, ".RData", sep = "")
    load(paste(folder, nm, sep = "/"))
    #M <- length(FkData)
    ### SAEM Fits ####

    cat("\n ________ SAEM fitting on dataset____", nm, "\n")
    fakedata <- FkData[[1]]
    #head(fakedata$data)
    # Generate a framework for non linear mixed model fitting
    Nframe = nlmer.frame (y ~ phi[1]/(1 + exp((phi[2] - X)/phi[3])),
                          X.formula = ~ x1, name.mixef = c("phi1", "phi2", "phi3"),
                          fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                          fixed.matrix = fakedata$A, cluster = ~ cluster,
                          random = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                          random.matrix = fakedata$B, data = fakedata$data)

    start = nlmer.start(frame = Nframe, smixvar = fitdistr, smixvar_b = fitdistr,
                        start = list(fixef = fixef))

    Desc <- paste(fitdistr, "fit on", nm)
    res <- vector(mode = "list", length = M)
    attr(res, "desc") <- Desc
    attr(res, "Algorithm") <- "SAEM"
    attr(res, "desc")
    attr(res, "Algorithm")

    for (d in 1:M) {
      cat("data", d, "\n")
      fakedata <- FkData[[d]]
      res[[d]] = ssmnnlmer(y ~ phi[1]/(1 + exp((phi[2] - X)/phi[3])),
                           X.formula = ~ x1, name.mixef = c("phi1", "phi2", "phi3"),
                           fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                           fixed.matrix = fakedata$A, cluster = ~ cluster,
                           random = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                           random.matrix = fakedata$B, data = fakedata$data,
                           smixvar = fitdistr, smixvar_b = fitdistr,
                           start = start, control = list(maxwarmit = 300,
                                                         warmtol = 1e-1,
                                                         maxsmoothit = 200,
                                                         smoothtol = 1e-4,
                                                         m_b = 10,
                                                         m_mcem = 150,
                                                         m_saem = 15,
                                                         verbose = TRUE))
    }
    save(res, file = paste(resFolder, paste("SAEM", fitdistr, nm, sep = "_") ,sep = "/"))
    Theta <- makeTheta(res, algo = "SAEM")

   # TrueTheta <- Theta[,1] + .5 # REMOVE AND DEFINE IT ABOVE !!!
    #TrueTheta <-c(200,700,350,0.5735393, 0.5735393, 0.5735393,0,0,0,1,1,1,0.5735393,0.5)
    TrueTheta <-c(200,700,350, #fixef
                  0, 0, 0, # It was symmetric data
                  1, 0, 0, 1, 0, 1, # Omegabar
                  0, #delta_o
                  .5)


    per1 <- performance(result = res, Theta = Theta, TrueTheta = TrueTheta)
    per2 <- efficiency(result = res, algo = "SAEM")
    simResult <- list(Theta = Theta, per1 = per1, per2 = per2)
    save(simResult, file = paste(resFolder, paste("simRes", "SAEM", fitdistr, nm, sep = "_") ,sep = "/"))

    ## EM Fits ####
    cat("\n ________ EM__", fitdistr, "__fitting on dataset____", nm, "\n")
    Desc <- paste(fitdistr, "fit on", nm)
    EMres <- vector(mode = "list", length = M)
    attr(EMres, "desc") <- Desc
    attr(EMres, "Algorithm") <- "EM"
    attr(EMres, "desc")
    attr(EMres, "Algorithm")

    for (d in 1:M) {
      cat("data", d, "\n")
      fakedata <- FkData[[d]]
      EMres[[d]] = pereEM(maxit = 4000, tol = 1e-4,
                          formula = y ~ phi[1]/(1 + exp((phi[2] - X[1])/phi[3])),
                          X.formula = ~ x1, fixef.name = c("phi1", "phi2", "phi3"),
                          cluster = ~ cluster, family = fitdistr, data = fakedata$data,
                          start = list(alpha = fixef, sigma2 = .5, D = diag(3),
                                       lambda =rep(0, 3), nu = 5), #nu = c(.1, .5),
                          skew = TRUE, pverbose = TRUE, st.err = TRUE, nuEstim = FALSE,
                          errorF = "warning")
    }
    save(EMres, file = paste(resFolder, paste("EM", fitdistr, nm, sep = "_") ,sep = "/"))
    Theta <- makeTheta(EMres, algo = "EM")

    #TrueTheta <- Theta[,1] + .5 # REMOVE AND DEFINE IT ABOVE !!!
    #TrueTheta <-c(200,700,350, 5,5,5, 0,0,0,1,1,1,0.5735393,0.5)
    TrueTheta <-c(200,700,350, #fixef
                  0, 0, 0, # It was symmetric data
                  .5,
                  1, 0, 0, 1, 0, 1 # Omegabar
                  )

    per1 <- performance(result = EMres, Theta = Theta, TrueTheta = TrueTheta)
    per2 <- efficiency(result = EMres, algo = "EM")
    simResult <- list(Theta = Theta, per1 = per1, per2 = per2)
    save(simResult, file = paste(resFolder, paste("simRes", "EM", fitdistr, nm, sep = "_") ,sep = "/"))

  }
}


# Skew t fits on all ####
fitdistr <- "Gamma"
##  RIS data ####
compl <- "RIS"
SMSN <- c("Ind", "Gamma")
for (distr in SMSN) {
  cat("\n", "\nOOOO ~", distr, "data fits", "~ OOOO\n")
  for (n in spSz) {
    nm <- paste(distr, n, compl, ".RData", sep = "")
    load(paste(folder, nm, sep = "/"))
    #M <- length(FkData)
    ### SAEM Fits ####

    cat("\n ________ SAEM fitting on dataset____", nm, "\n")
    fakedata <- FkData[[1]]
    #head(fakedata$data)
    # Generate a framework for non linear mixed model fitting
    Nframe = nlmer.frame (y ~ phi[1]/(1 + exp((phi[2] - X)/phi[3])),
                          X.formula = ~ x1, name.mixef = c("phi1", "phi2", "phi3"),
                          fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                          fixed.matrix = fakedata$A, cluster = ~ cluster,
                          random = list(phi1 = ~ 1, phi2 = ~ w1, phi3 = ~ w2),
                          random.matrix = fakedata$B, data = fakedata$data)

    start = nlmer.start(frame = Nframe, smixvar = fitdistr, smixvar_b = fitdistr,
                        start = list(fixef = fixef))

    Desc <- paste(fitdistr, "fit on", nm)
    res <- vector(mode = "list", length = M)
    attr(res, "desc") <- Desc
    attr(res, "Algorithm") <- "SAEM"
    attr(res, "desc")
    attr(res, "Algorithm")

    for (d in 1:M) {
      cat("data", d, "\n")
      fakedata <- FkData[[d]]
      res[[d]] = ssmnnlmer(y ~ phi[1]/(1 + exp((phi[2] - X)/phi[3])),
                           X.formula = ~ x1, name.mixef = c("phi1", "phi2", "phi3"),
                           fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                           fixed.matrix = fakedata$A, cluster = ~ cluster,
                           random = list(phi1 = ~ 1, phi2 = ~ w1, phi3 = ~ w2),
                           random.matrix = fakedata$B, data = fakedata$data,
                           smixvar = fitdistr, smixvar_b = fitdistr,
                           start = start, control = list(maxwarmit = 300,
                                                         warmtol = 1e-1,
                                                         maxsmoothit = 200,
                                                         smoothtol = 1e-4,
                                                         m_b = 10,
                                                         m_mcem = 150,
                                                         m_saem = 15,
                                                         verbose = TRUE))
    }
    save(res, file = paste(resFolder, paste("SAEM", fitdistr, nm, sep = "_") ,sep = "/"))
    Theta <- makeTheta(res, algo = "SAEM")

    #TrueTheta <- Theta[,1] + .5 # REMOVE AND DEFINE IT ABOVE !!!
    TrueTheta <-c(200,700,350, #fixef
                  0.4454354, 0.4454354, 0.4454354, 0.4454354, 0.4454354, # delta
                  0.8015873, -0.1984127, -0.1984127, -0.1984127, -0.1984127, 0.8015873, -0.1984127, -0.1984127, -0.1984127, 0.8015873, -0.1984127, -0.1984127,  0.8015873, -0.1984127, 0.8015873, # Omegabar
                  0.4902903, #delta_o
                  0.2596154)
    #Vec(MOmegabar(diag(5), lambda = rep(5,5)))
    #Mdelta(.5, lambda = 5)
    #MOmegabar(.5, lambda = 5)


    per1 <- performance(result = res, Theta = Theta, TrueTheta = TrueTheta)
    per2 <- efficiency(result = res, algo = "SAEM")
    simResult <- list(Theta = Theta, per1 = per1, per2 = per2)
    save(simResult, file = paste(resFolder, paste("simRes", "SAEM", fitdistr, nm, sep = "_") ,sep = "/"))

    ## EM Fits ####
    cat("\n ________ EM__", fitdistr, "__fitting on dataset____", nm, "\n")
    Desc <- paste(fitdistr, "fit on", nm)
    EMres <- vector(mode = "list", length = M)
    attr(EMres, "desc") <- Desc
    attr(EMres, "Algorithm") <- "EM"
    attr(EMres, "desc")
    attr(EMres, "Algorithm")

    for (d in 1:M) {
      cat("data", d, "\n")
      fakedata <- FkData[[d]]
      EMres[[d]] = pereEM(maxit = 4000, tol = 1e-4,
                          formula = y ~ phi[1]/(1 + exp((phi[2] - X[1])/phi[3])),
                          X.formula = ~ x1, fixef.name = c("phi1", "phi2", "phi3"),
                          cluster = ~ cluster, family = fitdistr, data = fakedata$data,
                          start = list(alpha = fixef, sigma2 = .5, D = diag(5),
                                       lambda =rep(0, 5), nu = 5), #nu = c(.1, .5),
                          skew = TRUE, pverbose = TRUE, st.err = TRUE, nuEstim = FALSE,
                          errorF = "warning")
    }
    save(EMres, file = paste(resFolder, paste("EM", fitdistr, nm, sep = "_") ,sep = "/"))
    Theta <- makeTheta(EMres, algo = "EM")

    #TrueTheta <- Theta[,1] + .5 # REMOVE AND DEFINE IT ABOVE !!!
    TrueTheta <-c(200,700,350, #fixef
                  0.4454354, 0.4454354, 0.4454354, 0.4454354, 0.4454354, # delta
                  0.2596154,
                  0.8015873, -0.1984127, -0.1984127, -0.1984127, -0.1984127, 0.8015873, -0.1984127, -0.1984127, -0.1984127, 0.8015873, -0.1984127, -0.1984127,  0.8015873, -0.1984127, 0.8015873)



    per1 <- performance(result = EMres, Theta = Theta, TrueTheta = TrueTheta)
    per2 <- efficiency(result = EMres, algo = "EM")
    simResult <- list(Theta = Theta, per1 = per1, per2 = per2)
    save(simResult, file = paste(resFolder, paste("simRes", "EM", fitdistr, nm, sep = "_") ,sep = "/"))

  }
}


## RI data ####
compl <- "RI"
SMSN <- c("Ind", "Gamma")
for (distr in SMSN) {
  cat("\n", "\nOOOO ~", distr, "data fits", "~ OOOO\n")
  for (n in spSz) {
    nm <- paste(distr, n, compl, ".RData", sep = "")
    load(paste(folder, nm, sep = "/"))
    #M <- length(FkData)
    ### SAEM Fits ####

    cat("\n ________ SAEM fitting on dataset____", nm, "\n")
    fakedata <- FkData[[1]]
    #head(fakedata$data)
    # Generate a framework for non linear mixed model fitting
    Nframe = nlmer.frame (y ~ phi[1]/(1 + exp((phi[2] - X)/phi[3])),
                          X.formula = ~ x1, name.mixef = c("phi1", "phi2", "phi3"),
                          fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                          fixed.matrix = fakedata$A, cluster = ~ cluster,
                          random = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                          random.matrix = fakedata$B, data = fakedata$data)

    start = nlmer.start(frame = Nframe, smixvar = fitdistr, smixvar_b = fitdistr,
                        start = list(fixef = fixef))

    Desc <- paste(fitdistr, "fit on", nm)
    res <- vector(mode = "list", length = M)
    attr(res, "desc") <- Desc
    attr(res, "Algorithm") <- "SAEM"
    attr(res, "desc")
    attr(res, "Algorithm")

    for (d in 1:M) {
      cat("data", d, "\n")
      fakedata <- FkData[[d]]
      res[[d]] = ssmnnlmer(y ~ phi[1]/(1 + exp((phi[2] - X)/phi[3])),
                           X.formula = ~ x1, name.mixef = c("phi1", "phi2", "phi3"),
                           fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                           fixed.matrix = fakedata$A, cluster = ~ cluster,
                           random = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                           random.matrix = fakedata$B, data = fakedata$data,
                           smixvar = fitdistr, smixvar_b = fitdistr,
                           start = start, control = list(maxwarmit = 300,
                                                         warmtol = 1e-1,
                                                         maxsmoothit = 200,
                                                         smoothtol = 1e-4,
                                                         m_b = 10,
                                                         m_mcem = 150,
                                                         m_saem = 15,
                                                         verbose = TRUE))
    }
    save(res, file = paste(resFolder, paste("SAEM", fitdistr, nm, sep = "_") ,sep = "/"))
    Theta <- makeTheta(res, algo = "SAEM")

    #TrueTheta <- Theta[,1] + .5 # REMOVE AND DEFINE IT ABOVE !!!
    #TrueTheta <-c(200,700,350,0.5735393, 0.5735393, 0.5735393,0,0,0,1,1,1,0.5735393,0.5)
    TrueTheta <-c(200,700,350, #fixed
                  0.5735393, 0.5735393, 0.5735393, #delta
                  0.6710526, -0.3289474, -0.3289474, 0.6710526, -0.3289474, 0.6710526, #Omegabar
                  0.4902903, 0.2596154)
    #Mdelta(omega = diag(3), lambda = rep(5,3))
    #Vec(MOmegabar(omega = diag(3), lambda = rep(5,3)))
    #Mdelta(omega = .5, lambda = 5)
    #MOmegabar(omega = .5, lambda = 5)

    per1 <- performance(result = res, Theta = Theta, TrueTheta = TrueTheta)
    per2 <- efficiency(result = res, algo = "SAEM")
    simResult <- list(Theta = Theta, per1 = per1, per2 = per2)
    save(simResult, file = paste(resFolder, paste("simRes", "SAEM", fitdistr, nm, sep = "_") ,sep = "/"))

    ## EM Fits ####
    cat("\n ________ EM__", fitdistr, "__fitting on dataset____", nm, "\n")
    Desc <- paste(fitdistr, "fit on", nm)
    EMres <- vector(mode = "list", length = M)
    attr(EMres, "desc") <- Desc
    attr(EMres, "Algorithm") <- "EM"
    attr(EMres, "desc")
    attr(EMres, "Algorithm")

    for (d in 1:M) {
      cat("data", d, "\n")
      fakedata <- FkData[[d]]
      EMres[[d]] = pereEM(maxit = 4000, tol = 1e-4,
                          formula = y ~ phi[1]/(1 + exp((phi[2] - X[1])/phi[3])),
                          X.formula = ~ x1, fixef.name = c("phi1", "phi2", "phi3"),
                          cluster = ~ cluster, family = fitdistr, data = fakedata$data,
                          start = list(alpha = fixef, sigma2 = .5, D = diag(3),
                                       lambda =rep(0, 3), nu = 5), #nu = c(.1, .5),
                          skew = TRUE, pverbose = TRUE, st.err = TRUE, nuEstim = FALSE,
                          errorF = "warning")
    }
    save(EMres, file = paste(resFolder, paste("EM", fitdistr, nm, sep = "_") ,sep = "/"))
    Theta <- makeTheta(EMres, algo = "EM")

    #TrueTheta <- Theta[,1] + .5 # REMOVE AND DEFINE IT ABOVE !!!
    #TrueTheta <-c(200,700,350,0.5735393, 0.5735393, 0.5735393,0,0,0,1,1,1,0.5735393,0.5)
    TrueTheta <-c(200,700,350, #fixed
                  0.5735393, 0.5735393, 0.5735393, #delta
                  0.2596154,
                  0.6710526, -0.3289474, -0.3289474, 0.6710526, -0.3289474, 0.6710526 #Omegabar
    )



    per1 <- performance(result = EMres, Theta = Theta, TrueTheta = TrueTheta)
    per2 <- efficiency(result = EMres, algo = "EM")
    simResult <- list(Theta = Theta, per1 = per1, per2 = per2)
    save(simResult, file = paste(resFolder, paste("simRes", "EM", fitdistr, nm, sep = "_") ,sep = "/"))

  }
}

# Skew Slash fits on all ####
fitdistr <- "Beta"
##  RIS data ####
compl <- "RIS"
SMSN <- c("Ind", "Gamma")
for (distr in SMSN) {
  cat("\n", "\nOOOO ~", distr, "data fits", "~ OOOO\n")
  for (n in spSz) {
    nm <- paste(distr, n, compl, ".RData", sep = "")
    load(paste(folder, nm, sep = "/"))
    #M <- length(FkData)
    ### SAEM Fits ####
    cat("\n ______ SAEM__", fitdistr, "__fitting on dataset___", nm, "\n")

    fakedata <- FkData[[1]]
    #head(fakedata$data)
    # Generate a framework for non linear mixed model fitting
    Nframe = nlmer.frame (y ~ phi[1]/(1 + exp((phi[2] - X)/phi[3])),
                          X.formula = ~ x1, name.mixef = c("phi1", "phi2", "phi3"),
                          fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                          fixed.matrix = fakedata$A, cluster = ~ cluster,
                          random = list(phi1 = ~ 1, phi2 = ~ w1, phi3 = ~ w2),
                          random.matrix = fakedata$B, data = fakedata$data)

    start = nlmer.start(frame = Nframe, smixvar = fitdistr, smixvar_b = fitdistr,
                        start = list(fixef = fixef))

    Desc <- paste(fitdistr, "fit on", nm)
    res <- vector(mode = "list", length = M)
    attr(res, "desc") <- Desc
    attr(res, "Algorithm") <- "SAEM"
    attr(res, "desc")
    attr(res, "Algorithm")

    for (d in 1:M) {
      cat("data", d, "\n")
      fakedata <- FkData[[d]]
      res[[d]] = ssmnnlmer(y ~ phi[1]/(1 + exp((phi[2] - X)/phi[3])),
                           X.formula = ~ x1, name.mixef = c("phi1", "phi2", "phi3"),
                           fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                           fixed.matrix = fakedata$A, cluster = ~ cluster,
                           random = list(phi1 = ~ 1, phi2 = ~ w1, phi3 = ~ w2),
                           random.matrix = fakedata$B, data = fakedata$data,
                           smixvar = fitdistr, smixvar_b = fitdistr,
                           start = start, control = list(maxwarmit = 300,
                                                         warmtol = 1e-1,
                                                         maxsmoothit = 200,
                                                         smoothtol = 1e-4,
                                                         m_b = 10,
                                                         m_mcem = 150,
                                                         m_saem = 15,
                                                         verbose = TRUE))
    }
    save(res, file = paste(resFolder, paste("SAEM", fitdistr, nm, sep = "_") ,sep = "/"))
    Theta <- makeTheta(res, algo = "SAEM")

    #TrueTheta <- Theta[,1] + .5 # REMOVE AND DEFINE IT ABOVE !!!
    TrueTheta <-c(200,700,350, #fixef
                  0.4454354, 0.4454354, 0.4454354, 0.4454354, 0.4454354, # delta
                  0.8015873, -0.1984127, -0.1984127, -0.1984127, -0.1984127, 0.8015873, -0.1984127, -0.1984127, -0.1984127, 0.8015873, -0.1984127, -0.1984127,  0.8015873, -0.1984127, 0.8015873, # Omegabar
                  0.4902903, #delta_o
                  0.2596154)


    per1 <- performance(result = res, Theta = Theta, TrueTheta = TrueTheta)
    per2 <- efficiency(result = res, algo = "SAEM")
    simResult <- list(Theta = Theta, per1 = per1, per2 = per2)
    save(simResult, file = paste(resFolder, paste("simRes", "SAEM", fitdistr, nm, sep = "_") ,sep = "/"))

    ## EM Fits ####
    cat("\n ________ EM__", fitdistr, "__fitting on dataset____", nm, "\n")
    Desc <- paste(fitdistr, "fit on", nm)
    EMres <- vector(mode = "list", length = M)
    attr(EMres, "desc") <- Desc
    attr(EMres, "Algorithm") <- "EM"
    attr(EMres, "desc")
    attr(EMres, "Algorithm")

    for (d in 1:M) {
      cat("data", d, "\n")
      fakedata <- FkData[[d]]
      EMres[[d]] = pereEM(maxit = 4000, tol = 1e-4,
                          formula = y ~ phi[1]/(1 + exp((phi[2] - X[1])/phi[3])),
                          X.formula = ~ x1, fixef.name = c("phi1", "phi2", "phi3"),
                          cluster = ~ cluster, family = fitdistr, data = fakedata$data,
                          start = list(alpha = fixef, sigma2 = .5, D = diag(5),
                                       lambda =rep(0, 5), nu = 5), #nu = c(.1, .5),
                          skew = TRUE, pverbose = TRUE, st.err = TRUE, nuEstim = FALSE,
                          errorF = "warning")
    }
    save(EMres, file = paste(resFolder, paste("EM", fitdistr, nm, sep = "_") ,sep = "/"))
    Theta <- makeTheta(EMres, algo = "EM")

    #TrueTheta <- Theta[,1] + .5 # REMOVE AND DEFINE IT ABOVE !!!
    TrueTheta <-c(200,700,350, #fixef
                  0.4454354, 0.4454354, 0.4454354, 0.4454354, 0.4454354, # delta
                  0.2596154,
                  0.8015873, -0.1984127, -0.1984127, -0.1984127, -0.1984127, 0.8015873, -0.1984127, -0.1984127, -0.1984127, 0.8015873, -0.1984127, -0.1984127,  0.8015873, -0.1984127, 0.8015873)


    per1 <- performance(result = EMres, Theta = Theta, TrueTheta = TrueTheta)
    per2 <- efficiency(result = EMres, algo = "EM")
    simResult <- list(Theta = Theta, per1 = per1, per2 = per2)
    save(simResult, file = paste(resFolder, paste("simRes", "EM", fitdistr, nm, sep = "_") ,sep = "/"))

  }
}


## RI data ####
compl <- "RI"
SMSN <- c("Ind", "Gamma")
for (distr in SMSN) {
  cat("\n", "\nOOOO ~", distr, "data fits", "~ OOOO\n")
  for (n in spSz) {
    nm <- paste(distr, n, compl, ".RData", sep = "")
    load(paste(folder, nm, sep = "/"))
    #M <- length(FkData)
    ### SAEM Fits ####

    cat("\n ______ SAEM__", fitdistr, "__fitting on dataset___", nm, "\n")
    fakedata <- FkData[[1]]
    #head(fakedata$data)
    # Generate a framework for non linear mixed model fitting
    Nframe = nlmer.frame (y ~ phi[1]/(1 + exp((phi[2] - X)/phi[3])),
                          X.formula = ~ x1, name.mixef = c("phi1", "phi2", "phi3"),
                          fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                          fixed.matrix = fakedata$A, cluster = ~ cluster,
                          random = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                          random.matrix = fakedata$B, data = fakedata$data)

    start = nlmer.start(frame = Nframe, smixvar = fitdistr, smixvar_b = fitdistr,
                        start = list(fixef = fixef))

    Desc <- paste(fitdistr, "fit on", nm)
    res <- vector(mode = "list", length = M)
    attr(res, "desc") <- Desc
    attr(res, "Algorithm") <- "SAEM"
    attr(res, "desc")
    attr(res, "Algorithm")

    for (d in 1:M) {
      cat("data", d, "\n")
      fakedata <- FkData[[d]]
      res[[d]] = ssmnnlmer(y ~ phi[1]/(1 + exp((phi[2] - X)/phi[3])),
                           X.formula = ~ x1, name.mixef = c("phi1", "phi2", "phi3"),
                           fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                           fixed.matrix = fakedata$A, cluster = ~ cluster,
                           random = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                           random.matrix = fakedata$B, data = fakedata$data,
                           smixvar = fitdistr, smixvar_b = fitdistr,
                           start = start, control = list(maxwarmit = 300,
                                                         warmtol = 1e-1,
                                                         maxsmoothit = 200,
                                                         smoothtol = 1e-4,
                                                         m_b = 10,
                                                         m_mcem = 150,
                                                         m_saem = 15,
                                                         verbose = TRUE))
    }
    save(res, file = paste(resFolder, paste("SAEM", fitdistr, nm, sep = "_") ,sep = "/"))
    Theta <- makeTheta(res, algo = "SAEM")

    #TrueTheta <- Theta[,1] + .5 # REMOVE AND DEFINE IT ABOVE !!!
    #TrueTheta <-c(200,700,350,0.5735393, 0.5735393, 0.5735393,0,0,0,1,1,1,0.5735393,0.5)
    TrueTheta <-c(200,700,350, #fixed
                  0.5735393, 0.5735393, 0.5735393, #delta
                  0.6710526, -0.3289474, -0.3289474, 0.6710526, -0.3289474, 0.6710526, #Omegabar
                  0.4902903, 0.2596154)


    per1 <- performance(result = res, Theta = Theta, TrueTheta = TrueTheta)
    per2 <- efficiency(result = res, algo = "SAEM")
    simResult <- list(Theta = Theta, per1 = per1, per2 = per2)
    save(simResult, file = paste(resFolder, paste("simRes", "SAEM", fitdistr, nm, sep = "_") ,sep = "/"))

    ## EM Fits ####
    cat("\n ________ EM__", fitdistr, "__fitting on dataset____", nm, "\n")
    Desc <- paste(fitdistr, "fit on", nm)
    EMres <- vector(mode = "list", length = M)
    attr(EMres, "desc") <- Desc
    attr(EMres, "Algorithm") <- "EM"
    attr(EMres, "desc")
    attr(EMres, "Algorithm")

    for (d in 1:M) {
      cat("data", d, "\n")
      fakedata <- FkData[[d]]
      EMres[[d]] = pereEM(maxit = 4000, tol = 1e-4,
                          formula = y ~ phi[1]/(1 + exp((phi[2] - X[1])/phi[3])),
                          X.formula = ~ x1, fixef.name = c("phi1", "phi2", "phi3"),
                          cluster = ~ cluster, family = fitdistr, data = fakedata$data,
                          start = list(alpha = fixef, sigma2 = .5, D = diag(3),
                                       lambda =rep(0, 3), nu = 5), #nu = c(.1, .5),
                          skew = TRUE, pverbose = TRUE, st.err = TRUE, nuEstim = FALSE,
                          errorF = "warning")
    }
    save(EMres, file = paste(resFolder, paste("EM", fitdistr, nm, sep = "_") ,sep = "/"))
    Theta <- makeTheta(EMres, algo = "EM")

    #TrueTheta <- Theta[,1] + .5 # REMOVE AND DEFINE IT ABOVE !!!
    #TrueTheta <-c(200,700,350,0.5735393, 0.5735393, 0.5735393,0,0,0,1,1,1,0.5735393,0.5)
    TrueTheta <-c(200,700,350, #fixed
                  0.5735393, 0.5735393, 0.5735393, #delta
                  0.2596154,
                  0.6710526, -0.3289474, -0.3289474, 0.6710526, -0.3289474, 0.6710526 #Omegabar
    )

    per1 <- performance(result = EMres, Theta = Theta, TrueTheta = TrueTheta)
    per2 <- efficiency(result = EMres, algo = "EM")
    simResult <- list(Theta = Theta, per1 = per1, per2 = per2)
    save(simResult, file = paste(resFolder, paste("simRes", "EM", fitdistr, nm, sep = "_") ,sep = "/"))

  }
}

# Skew Contaminated fits on all ####
fitdistr <- "Bin"
##  RIS data ####
compl <- "RIS"
SMSN <- c("Ind", "Gamma")
for (distr in SMSN) {
  cat("\n", "\nOOOO ~", distr, "data fits", "~ OOOO\n")
  for (n in spSz) {
    nm <- paste(distr, n, compl, ".RData", sep = "")
    load(paste(folder, nm, sep = "/"))
    #M <- length(FkData)
    ### SAEM Fits ####
    cat("\n ______ SAEM__", fitdistr, "__fitting on dataset___", nm, "\n")

    fakedata <- FkData[[1]]
    #head(fakedata$data)
    # Generate a framework for non linear mixed model fitting
    Nframe = nlmer.frame (y ~ phi[1]/(1 + exp((phi[2] - X)/phi[3])),
                          X.formula = ~ x1, name.mixef = c("phi1", "phi2", "phi3"),
                          fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                          fixed.matrix = fakedata$A, cluster = ~ cluster,
                          random = list(phi1 = ~ 1, phi2 = ~ w1, phi3 = ~ w2),
                          random.matrix = fakedata$B, data = fakedata$data)

    start = nlmer.start(frame = Nframe, smixvar = fitdistr, smixvar_b = fitdistr,
                        start = list(fixef = fixef))

    Desc <- paste(fitdistr, "fit on", nm)
    res <- vector(mode = "list", length = M)
    attr(res, "desc") <- Desc
    attr(res, "Algorithm") <- "SAEM"
    attr(res, "desc")
    attr(res, "Algorithm")

    for (d in 1:M) {
      cat("data", d, "\n")
      fakedata <- FkData[[d]]
      res[[d]] = ssmnnlmer(y ~ phi[1]/(1 + exp((phi[2] - X)/phi[3])),
                           X.formula = ~ x1, name.mixef = c("phi1", "phi2", "phi3"),
                           fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                           fixed.matrix = fakedata$A, cluster = ~ cluster,
                           random = list(phi1 = ~ 1, phi2 = ~ w1, phi3 = ~ w2),
                           random.matrix = fakedata$B, data = fakedata$data,
                           smixvar = fitdistr, smixvar_b = fitdistr,
                           start = start, control = list(maxwarmit = 300,
                                                         warmtol = 1e-1,
                                                         maxsmoothit = 200,
                                                         smoothtol = 1e-4,
                                                         m_b = 10,
                                                         m_mcem = 150,
                                                         m_saem = 15,
                                                         verbose = TRUE))
    }
    save(res, file = paste(resFolder, paste("SAEM", fitdistr, nm, sep = "_") ,sep = "/"))
    Theta <- makeTheta(res, algo = "SAEM")

    #TrueTheta <- Theta[,1] + .5 # REMOVE AND DEFINE IT ABOVE !!!
    TrueTheta <-c(200,700,350, #fixef
                  0.4454354, 0.4454354, 0.4454354, 0.4454354, 0.4454354, # delta
                  0.8015873, -0.1984127, -0.1984127, -0.1984127, -0.1984127, 0.8015873, -0.1984127, -0.1984127, -0.1984127, 0.8015873, -0.1984127, -0.1984127,  0.8015873, -0.1984127, 0.8015873, # Omegabar
                  0.4902903, #delta_o
                  0.2596154)


    per1 <- performance(result = res, Theta = Theta, TrueTheta = TrueTheta)
    per2 <- efficiency(result = res, algo = "SAEM")
    simResult <- list(Theta = Theta, per1 = per1, per2 = per2)
    save(simResult, file = paste(resFolder, paste("simRes", "SAEM", fitdistr, nm, sep = "_") ,sep = "/"))

    ## EM Fits ####
    cat("\n ________ EM__", fitdistr, "__fitting on dataset____", nm, "\n")
    Desc <- paste(fitdistr, "fit on", nm)
    EMres <- vector(mode = "list", length = M)
    attr(EMres, "desc") <- Desc
    attr(EMres, "Algorithm") <- "EM"
    attr(EMres, "desc")
    attr(EMres, "Algorithm")

    for (d in 1:M) {
      cat("data", d, "\n")
      fakedata <- FkData[[d]]
      EMres[[d]] = pereEM(maxit = 4000, tol = 1e-4,
                          formula = y ~ phi[1]/(1 + exp((phi[2] - X[1])/phi[3])),
                          X.formula = ~ x1, fixef.name = c("phi1", "phi2", "phi3"),
                          cluster = ~ cluster, family = fitdistr, data = fakedata$data,
                          start = list(alpha = fixef, sigma2 = .5, D = diag(5),
                                       lambda =rep(0, 5), nu = c(0.5, 0.5)),
                          skew = TRUE, pverbose = TRUE, st.err = TRUE, nuEstim = FALSE,
                          errorF = "warning")
    }
    save(EMres, file = paste(resFolder, paste("EM", fitdistr, nm, sep = "_") ,sep = "/"))
    Theta <- makeTheta(EMres, algo = "EM")

    #TrueTheta <- Theta[,1] + .5 # REMOVE AND DEFINE IT ABOVE !!!
    TrueTheta <-c(200,700,350, #fixef
                  0.4454354, 0.4454354, 0.4454354, 0.4454354, 0.4454354, # delta
                  0.2596154,
                  0.8015873, -0.1984127, -0.1984127, -0.1984127, -0.1984127, 0.8015873, -0.1984127, -0.1984127, -0.1984127, 0.8015873, -0.1984127, -0.1984127,  0.8015873, -0.1984127, 0.8015873)


    per1 <- performance(result = EMres, Theta = Theta, TrueTheta = TrueTheta)
    per2 <- efficiency(result = EMres, algo = "EM")
    simResult <- list(Theta = Theta, per1 = per1, per2 = per2)
    save(simResult, file = paste(resFolder, paste("simRes", "EM", fitdistr, nm, sep = "_") ,sep = "/"))

  }
}


## RI data ####
compl <- "RI"
SMSN <- c("Ind", "Gamma")
for (distr in SMSN) {
  cat("\n", "\nOOOO ~", distr, "data fits", "~ OOOO\n")
  for (n in spSz) {
    nm <- paste(distr, n, compl, ".RData", sep = "")
    load(paste(folder, nm, sep = "/"))
    #M <- length(FkData)
    ### SAEM Fits ####

    cat("\n ______ SAEM__", fitdistr, "__fitting on dataset___", nm, "\n")
    fakedata <- FkData[[1]]
    #head(fakedata$data)
    # Generate a framework for non linear mixed model fitting
    Nframe = nlmer.frame (y ~ phi[1]/(1 + exp((phi[2] - X)/phi[3])),
                          X.formula = ~ x1, name.mixef = c("phi1", "phi2", "phi3"),
                          fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                          fixed.matrix = fakedata$A, cluster = ~ cluster,
                          random = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                          random.matrix = fakedata$B, data = fakedata$data)

    start = nlmer.start(frame = Nframe, smixvar = fitdistr, smixvar_b = fitdistr,
                        start = list(fixef = fixef))

    Desc <- paste(fitdistr, "fit on", nm)
    res <- vector(mode = "list", length = M)
    attr(res, "desc") <- Desc
    attr(res, "Algorithm") <- "SAEM"
    attr(res, "desc")
    attr(res, "Algorithm")

    for (d in 1:M) {
      cat("data", d, "\n")
      fakedata <- FkData[[d]]
      res[[d]] = ssmnnlmer(y ~ phi[1]/(1 + exp((phi[2] - X)/phi[3])),
                           X.formula = ~ x1, name.mixef = c("phi1", "phi2", "phi3"),
                           fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                           fixed.matrix = fakedata$A, cluster = ~ cluster,
                           random = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                           random.matrix = fakedata$B, data = fakedata$data,
                           smixvar = fitdistr, smixvar_b = fitdistr,
                           start = start, control = list(maxwarmit = 300,
                                                         warmtol = 1e-1,
                                                         maxsmoothit = 200,
                                                         smoothtol = 1e-4,
                                                         m_b = 10,
                                                         m_mcem = 150,
                                                         m_saem = 15,
                                                         verbose = TRUE))
    }
    save(res, file = paste(resFolder, paste("SAEM", fitdistr, nm, sep = "_") ,sep = "/"))
    Theta <- makeTheta(res, algo = "SAEM")

    #TrueTheta <- Theta[,1] + .5 # REMOVE AND DEFINE IT ABOVE !!!
    #TrueTheta <-c(200,700,350,0.5735393, 0.5735393, 0.5735393,0,0,0,1,1,1,0.5735393,0.5)
    TrueTheta <-c(200,700,350, #fixed
                  0.5735393, 0.5735393, 0.5735393, #delta
                  0.6710526, -0.3289474, -0.3289474, 0.6710526, -0.3289474, 0.6710526, #Omegabar
                  0.4902903, 0.2596154)

    per1 <- performance(result = res, Theta = Theta, TrueTheta = TrueTheta)
    per2 <- efficiency(result = res, algo = "SAEM")
    simResult <- list(Theta = Theta, per1 = per1, per2 = per2)
    save(simResult, file = paste(resFolder, paste("simRes", "SAEM", fitdistr, nm, sep = "_") ,sep = "/"))

    ## EM Fits ####
    cat("\n ________ EM__", fitdistr, "__fitting on dataset____", nm, "\n")
    Desc <- paste(fitdistr, "fit on", nm)
    EMres <- vector(mode = "list", length = M)
    attr(EMres, "desc") <- Desc
    attr(EMres, "Algorithm") <- "EM"
    attr(EMres, "desc")
    attr(EMres, "Algorithm")

    for (d in 1:M) {
      cat("data", d, "\n")
      fakedata <- FkData[[d]]
      EMres[[d]] = pereEM(maxit = 4000, tol = 1e-4,
                          formula = y ~ phi[1]/(1 + exp((phi[2] - X[1])/phi[3])),
                          X.formula = ~ x1, fixef.name = c("phi1", "phi2", "phi3"),
                          cluster = ~ cluster, family = fitdistr, data = fakedata$data,
                          start = list(alpha = fixef, sigma2 = .5, D = diag(3),
                                       lambda =rep(0, 3), nu = c(0.5, 0.5)), #nu = c(.1, .5),
                          skew = TRUE, pverbose = TRUE, st.err = TRUE, nuEstim = FALSE,
                          errorF = "warning")
    }
    save(EMres, file = paste(resFolder, paste("EM", fitdistr, nm, sep = "_") ,sep = "/"))
    Theta <- makeTheta(EMres, algo = "EM")

    #TrueTheta <- Theta[,1] + .5 # REMOVE AND DEFINE IT ABOVE !!!
    #TrueTheta <-c(200,700,350,0.5735393, 0.5735393, 0.5735393,0,0,0,1,1,1,0.5735393,0.5)
    TrueTheta <-c(200,700,350, #fixed
                  0.5735393, 0.5735393, 0.5735393, #delta
                  0.2596154,
                  0.6710526, -0.3289474, -0.3289474, 0.6710526, -0.3289474, 0.6710526 #Omegabar
    )

    per1 <- performance(result = EMres, Theta = Theta, TrueTheta = TrueTheta)
    per2 <- efficiency(result = EMres, algo = "EM")
    simResult <- list(Theta = Theta, per1 = per1, per2 = per2)
    save(simResult, file = paste(resFolder, paste("simRes", "EM", fitdistr, nm, sep = "_") ,sep = "/"))

  }
}

