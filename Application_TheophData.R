
# This is the new version of fnlmer, without SNP codes
library(skewmreg)


# Load the Theophylline data from the package "datasets"
data(Theoph, package = "datasets")
attach(Theoph)
head(Theoph,10)
dim(Theoph)
###############################################################
####################     N-NLMEM      #########################
###############################################################
# Expression of the nonlinear function
nform <- ~ Dose * exp(phi1 + phi2 - phi3) * (exp(-exp(phi1) * Time) - exp(-exp(phi2) * Time))/
           (exp(phi2) - exp(phi1)) # Define the nonlinear function
nfun <- deriv(nform, namevec = c("phi1","phi2","phi3"),
              function.arg = c("Dose","Time","phi1","phi2","phi3"),
              hessian = FALSE)

# To find starting values for regression parameters, first fit a nonlinear
# regression model (only fixed effects) using non linear least squares (nls).
# and Pereira et al. (2018): phi1=-2.52; phi2=.45; phi3=-3.26 as starting value in nls
nlsfit0 = nls (formula = conc ~ nfun (Dose, Time, phi1, phi2, phi3),
               data = Theoph,
               control = nls.control(maxiter = 500),
               start = list(phi1 = -2.52, phi2 = .45, phi3 = -3.26))

sumnlsfit0 = summary(nlsfit0)
sumnlsfit0

# Using 'nlme' to fit the Gaussian model
Timenlme = system.time({
  nlmefit = nlme::nlme (model = conc ~ nfun (Dose, Time, phi1, phi2, phi3),
                        data = Theoph,
                        fixed = list(phi1 ~ 1, phi2 ~ 1, phi3 ~ 1),
                        random = list(phi1 ~ 1, phi2 ~ 1, phi3 ~ 1),
                        groups = ~ Subject, correlation = nlme::corSymm(),
                        control = nlme::nlmeControl(maxIter = 500,
                                                    msMaxIter = 100,
                                                    returnObject = TRUE),
                        method = "REML",
                        start = list(fixed = c(-2.52, .4, -3.25)))
})

nlmefit$numIter

# Elapsed time
Timenlme[3]/60

# Fixed effects (estimates)
nlmefit$coefficients$fixed #
# Residual variance (estimate)
nlmefit$sigma^2

# Covariance matrix of random effects (estimates)
sv = nlme::VarCorr(nlmefit)
sds <- as.numeric(sv[1:3,2])
sdc <- vech2mat(as.numeric(c(sv[2:3,3], sv[3,4])), diag = FALSE)
Sigma0 <-  diag(sds) %*% sdc %*% diag(sds) # lme4
Sigma0

# Table of coefficients
sumnlmefit = summary(nlmefit)
sumnlmefit

# Log-like and AIC
logLik(nlmefit)[1]
AIC(nlmefit)
#save(nlmefit, file = "AppN.RData")

###############################################################
####################     SN-NLMEM      ########################
###############################################################
library(skewmreg)

A = diag(3) # binary matrix indicating the design of fixed effects
B = diag(3) # binary matrix indicating the design of random effects


Nframe = nlmer.frame (conc ~ X[1] * exp(phi[1] + phi[2] - phi[3]) *
                        (exp(-exp(phi[1]) * X[2]) - exp(-exp(phi[2]) * X[2]))/
                        (exp(phi[2]) - exp(phi[1])),
                      X.formula = ~ Dose + Time,
                      name.mixef = c("phi1", "phi2", "phi3"),
                      fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                      fixed.matrix = A, cluster = ~ Subject,
                      random = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                      random.matrix = B, data = Theoph)

# Starting values for SN model
nlmeranef <- matrix(0, 12, 3) # A matrix to receive the predicted random effects from the nlme fit
nlmeranef[as.numeric(rownames(nlmefit$coefficients$random$Subject)),] <-
  nlmefit$coefficients$random$Subject # use estimates of random effects from nlme
snstart = nlmer.start(frame = Nframe, family = "ssmn",
                      smixvar = "SNorm", smixvar_b = "SNorm",
                      start = list(fixef = nlmefit$coefficients$fixed,
                                   omegabar2 = nlmefit$sigma^2,
                                   deltao = 0.1,
                                   Omegabar_b = Sigma0,
                                   delta_b = rep(0.1, 3),
                                   ranef = nlmeranef))
snstart$fixef
snstart$omegabar2

TimeSN = system.time({
  theoSN <- ssmnnlmer(conc ~ X[1] * exp(phi[1]+.1 + phi[2] - phi[3]) *
                        (exp(-exp(phi[1]+.1) * X[2]) - exp(-exp(phi[2]) * X[2]))/
                        (exp(phi[2]) - exp(phi[1]+.1)),
                      X.formula = ~ Dose + Time,
                      name.mixef = c("phi1", "phi2", "phi3"),
                      fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                      fixed.matrix = A, cluster = ~ Subject,
                      random = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                      random.matrix = B, data = Theoph,
                      smixvar = "SNorm", smixvar_b = "SNorm",
                      start = snstart,
                      method = "ssmnnlmer.saem",
                      control = list(maxwarmit = 300,
                                     warmtol = 1e-1,
                                     maxsmoothit = 500,
                                     smoothtol = 1e-4,
                                     m_b = 10,
                                     m_mcem = 150,
                                     m_saem = 15,
                                     verbose = TRUE))
})

theoSN$fit$param$omegabar2
#save(theoSN, file = "theoSN.RData")
theoSN$fit$mdata
# Elapsed time
TimeSN[3]/60

# Fixed effects
snstart$fixef # starting value
theoSN$fit$param$fixef # SN MNM

# Residual variance
Rsn = varcovssmn(delta = theoSN$fit$param$deltao,
           Omegabar = theoSN$fit$param$omegabar2,
           mu = theoSN$fit$param$ctau_e1 * theoSN$fit$param$deltao,
           smixvar = "SNorm")$covariance

# (Co)Variance of random effects
Vaesn = varcovssmn(delta = theoSN$fit$param$delta_b,
           Omegabar = theoSN$fit$param$Omegabar_b,
           mu = -theoSN$fit$param$ctau_e1 * theoSN$fit$param$delta_b,
           smixvar = "SNorm")$covariance


###############################################################
####################     ST-NLMEM      ########################
###############################################################

ststart = nlmer.start(frame = Nframe, family = "ssmn",
                      smixvar = "Gamma", smixvar_b = "Gamma",
                      start = list(fixef = theoSN$fit$param$fixef,
                                   nu = 3.49, nu_b = 3.49, # From Pereira et al. (2018)
                                   omegabar2 = theoSN$fit$param$omegabar2,
                                   deltao = theoSN$fit$param$deltao,
                                   Omegabar_b = theoSN$fit$param$Omegabar_b,
                                   delta_b = theoSN$fit$param$delta_b,
                                   ranef = nlmeranef))

Timenlme = system.time({
  theoST <- ssmnnlmer(conc ~ X[1] * exp(phi[1]+.1 + phi[2] - phi[3]) *
                        (exp(-exp(phi[1]+.1) * X[2]) - exp(-exp(phi[2]) * X[2]))/
                        (exp(phi[2]) - exp(phi[1]+.1)),
                      X.formula = ~ Dose + Time,
                      name.mixef = c("phi1", "phi2", "phi3"),
                      fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                      fixed.matrix = A, cluster = ~ Subject,
                      random = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                      random.matrix = B, data = Theoph,
                      smixvar = "Gamma", smixvar_b = "Gamma",
                      method = "ssmnnlmer.saem",
                      start = ststart,
                      control = list(update.df = TRUE,
                                     maxwarmit = 300,
                                     warmtol = 1e-1,
                                     maxsmoothit = 500,
                                     smoothtol = 1e-4,
                                     m_b = 10,
                                     m_mcem = 150,
                                     m_saem = 15,
                                     verbose = TRUE))
})
#save(theoST, file = "theoST.RData")

# Elapsed time
Timenlme[3]/60

# Fixed effects
theoST$fit$param$fixef # ST MNM by SAEM

# Residual variance
varcovssmn(delta = theoST$fit$param$deltao,
           Omegabar = theoST$fit$param$omegabar2,
           mu = theoST$fit$param$ctau_e1 * theoST$fit$param$deltao,
           smixvar = "Gamma",
           df = theoST$fit$param$nu)$covariance # TRUE value

# (Co)Variance of random effects

varcovssmn(delta = theoST$fit$param$delta_b,
           Omegabar = theoST$fit$param$Omegabar,
           mu = theoST$fit$param$ctau_e1 * theoST$fit$param$delta_b,
           smixvar = "Gamma",
           df = theoST$fit$param$nu_b)$covariance # True values




###############################################################
####################     SCN-NLMEM      #######################
###############################################################
scnstart = nlmer.start(frame = Nframe, family = "ssmn",
                       smixvar = "Bin", smixvar_b = "Bin",
                       start = list(fixef = c(-1.538, 0.487, -2), #theoSN$fit$param$fixef
                                    nu = c(0.2, 0.2), nu_b=c(0.2, 0.2),
                                    omegabar2 = 0.467 , #theoSN$fit$param$omegabar2
                                    deltao = 0.075, #ltheoSN$fit$param$deltao
                                   # Omegabar_b = theoSN$fit$param$Omegabar_b,
                                    delta_b = c(0.082, 0.9, -0.007) ,# theoSN$fit$param$delta_b
                                    ranef = nlmeranef))

Timenlme = system.time({
  theoSCN <- ssmnnlmer(conc ~ X[1] * exp(phi[1]+.1 + phi[2] - phi[3]) *
                         (exp(-exp(phi[1]+.1) * X[2]) - exp(-exp(phi[2]) * X[2]))/
                         (exp(phi[2]) - exp(phi[1]+.1)),
                       X.formula = ~ Dose + Time,
                       name.mixef = c("phi1", "phi2", "phi3"),
                       fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                       fixed.matrix = A, cluster = ~ Subject,
                       random = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                       random.matrix = B, data = Theoph,
                       smixvar = "Bin", smixvar_b = "Bin",
                       start = scnstart,
                       control = list(update.df = TRUE,
                                      maxwarmit = 300,
                                      warmtol = 1e-1,
                                      maxsmoothit = 500,
                                      smoothtol = 1e-4,
                                      m_b = 10,
                                      m_mcem = 150,
                                      m_saem = 15,
                                      verbose = TRUE))
})

Timenlme[3]/60

#save(theoSCN, file = "theoSCN.RData")

# Fixed effects
theoSCN$fit$param$fixef # SCN MNM by SAEM

# Residual variance
varcovssmn(delta = theoSCN$fit$param$deltao,
           Omegabar = theoSCN$fit$param$omegabar2,
           mu = theoSCN$fit$param$ctau_e1 * theoSCN$fit$param$deltao,
           smixvar = "Bin",
           df = theoSCN$fit$param$nu)$covariance # TRUE value

# (Co)Variance of random effects

varcovssmn(delta = theoSCN$fit$param$delta_b,
           Omegabar = theoSCN$fit$param$Omegabar,
           mu = theoSCN$fit$param$ctau_e1 * theoSCN$fit$param$delta_b,
           smixvar = "Bin",
           df = theoSCN$fit$param$nu_b)$covariance # True values



###############################################################
####################     SSL-NLMEM      #######################
###############################################################
sslstart = nlmer.start(frame = Nframe, family = "ssmn",
                      smixvar = "Beta", smixvar_b = "Beta",
                      start = list(fixef = theoSN$fit$param$fixef,
                                   nu = 3.49, nu_b = 3.49, # From Pereira et al. (2018)
                                   omegabar2 = theoSN$fit$param$omegabar2,
                                   deltao = theoSN$fit$param$deltao,
                                   Omegabar_b = theoSN$fit$param$Omegabar_b,
                                   delta_b = theoSN$fit$param$delta_b,
                                   ranef = nlmeranef))

Timenlme = system.time({
  theoSSL <- ssmnnlmer(conc ~ X[1] * exp(phi[1]+.1 + phi[2] - phi[3]) *
                        (exp(-exp(phi[1]+.1) * X[2]) - exp(-exp(phi[2]) * X[2]))/
                        (exp(phi[2]) - exp(phi[1]+.1)),
                      X.formula = ~ Dose + Time,
                      name.mixef = c("phi1", "phi2", "phi3"),
                      fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                      fixed.matrix = A, cluster = ~ Subject,
                      random = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                      random.matrix = B, data = Theoph,
                      smixvar = "Beta", smixvar_b = "Beta",
                      method = "ssmnnlmer.saem",
                      start = ststart,
                      control = list(update.df = TRUE,
                                     maxwarmit = 300,
                                     warmtol = 1e-1,
                                     maxsmoothit = 500,
                                     smoothtol = 1e-4,
                                     m_b = 10,
                                     m_mcem = 150,
                                     m_saem = 15,
                                     verbose = TRUE))
})
#save(theoSSL, file = "theoSSL.RData")

# Elapsed time
Timenlme[3]/60

# Fixed effects
theoSSL$fit$param$fixef # SSL MNM by SAEM

# Residual variance
varcovssmn(delta = theoSSL$fit$param$deltao,
           Omegabar = theoSSL$fit$param$omegabar2,
           mu = theoSSL$fit$param$ctau_e1 * theoSSL$fit$param$deltao,
           smixvar = "Beta",
           df = theoSSL$fit$param$nu)$covariance # TRUE value

# (Co)Variance of random effects

varcovssmn(delta = theoSSL$fit$param$delta_b,
           Omegabar = theoSSL$fit$param$Omegabar,
           mu = theoSSL$fit$param$ctau_e1 * theoSSL$fit$param$delta_b,
           smixvar = "Beta",
           df = theoSSL$fit$param$nu_b)$covariance # True values


###############################################################
####################     ST_ SSL-NLMEM      ########################
###############################################################

stslstart = nlmer.start(frame = Nframe, family = "ssmn",
                      smixvar = "Gamma", smixvar_b = "Beta",
                      start = list(fixef = theoSN$fit$param$fixef,
                                   nu = 3.49, nu_b = 3.49, # From Pereira et al. (2018)
                                   omegabar2 = theoSN$fit$param$omegabar2,
                                   deltao = theoSN$fit$param$deltao,
                                   Omegabar_b = theoSN$fit$param$Omegabar_b,
                                   delta_b = theoSN$fit$param$delta_b,
                                   ranef = nlmeranef))

Timenlme = system.time({
  theoSTSL <- ssmnnlmer(conc ~ X[1] * exp(phi[1]+.1 + phi[2] - phi[3]) *
                        (exp(-exp(phi[1]+.1) * X[2]) - exp(-exp(phi[2]) * X[2]))/
                        (exp(phi[2]) - exp(phi[1]+.1)),
                      X.formula = ~ Dose + Time,
                      name.mixef = c("phi1", "phi2", "phi3"),
                      fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                      fixed.matrix = A, cluster = ~ Subject,
                      random = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                      random.matrix = B, data = Theoph,
                      smixvar = "Gamma", smixvar_b = "Beta",
                      method = "ssmnnlmer.saem",
                      start = stslstart,
                      control = list(update.df = TRUE,
                                     maxwarmit = 300,
                                     warmtol = 1e-1,
                                     maxsmoothit = 500,
                                     smoothtol = 1e-4,
                                     m_b = 10,
                                     m_mcem = 150,
                                     m_saem = 15,
                                     verbose = TRUE))
})
#save(theoSTSL, file = "theoSTSL.RData")

# Elapsed time
Timenlme[3]/60

# Fixed effects
theoSTSL$fit$param$fixef # STSL MNM by SAEM

# Residual variance
varcovssmn(delta = theoSTSL$fit$param$deltao,
           Omegabar = theoSTSL$fit$param$omegabar2,
           mu = theoSTSL$fit$param$ctau_e1 * theoSTSL$fit$param$deltao,
           smixvar = "Gamma",
           df = theoSTSL$fit$param$nu)$covariance # TRUE value

# (Co)Variance of random effects

varcovssmn(delta = theoSTSL$fit$param$delta_b,
           Omegabar = theoSTSL$fit$param$Omegabar,
           mu = theoSTSL$fit$param$ctau_e1 * theoSTSL$fit$param$delta_b,
           smixvar = "Gamma",
           df = theoSTSL$fit$param$nu_b)$covariance # True values

###############################################################
####################     SSL_ST-NLMEM      ########################
###############################################################

ssltstart = nlmer.start(frame = Nframe, family = "ssmn",
                        smixvar = "Beta", smixvar_b = "Gamma",
                        start = list(fixef = theoSN$fit$param$fixef,
                                     nu = 3.49, nu_b = 3.49, # From Pereira et al. (2018)
                                     omegabar2 = theoSN$fit$param$omegabar2,
                                     deltao = theoSN$fit$param$deltao,
                                     Omegabar_b = theoSN$fit$param$Omegabar_b,
                                     delta_b = theoSN$fit$param$delta_b,
                                     ranef = nlmeranef))

Timenlme = system.time({
  theoSSLT <- ssmnnlmer(conc ~ X[1] * exp(phi[1]+.1 + phi[2] - phi[3]) *
                          (exp(-exp(phi[1]+.1) * X[2]) - exp(-exp(phi[2]) * X[2]))/
                          (exp(phi[2]) - exp(phi[1]+.1)),
                        X.formula = ~ Dose + Time,
                        name.mixef = c("phi1", "phi2", "phi3"),
                        fixed = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                        fixed.matrix = A, cluster = ~ Subject,
                        random = list(phi1 = ~ 1, phi2 = ~ 1, phi3 = ~ 1),
                        random.matrix = B, data = Theoph,
                        smixvar = "Beta", smixvar_b = "Gamma",
                        method = "ssmnnlmer.saem",
                        start = stslstart,
                        control = list(update.df = TRUE,
                                       maxwarmit = 300,
                                       warmtol = 1e-1,
                                       maxsmoothit = 500,
                                       smoothtol = 1e-4,
                                       m_b = 10,
                                       m_mcem = 150,
                                       m_saem = 15,
                                       verbose = TRUE))
})
#save(theoSSLT, file = "theoSSLT.RData")

# Elapsed time
Timenlme[3]/60

# Fixed effects
theoSSLT$fit$param$fixef # STSL MNM by SAEM

# Residual variance
varcovssmn(delta = theoSSLT$fit$param$deltao,
           Omegabar = theoSSLT$fit$param$omegabar2,
           mu = theoSSLT$fit$param$ctau_e1 * theoSSLT$fit$param$deltao,
           smixvar = "Beta",
           df = theoSSLT$fit$param$nu)$covariance # TRUE value

# (Co)Variance of random effects

varcovssmn(delta = theoSSLT$fit$param$delta_b,
           Omegabar = theoSSLT$fit$param$Omegabar,
           mu = theoSSLT$fit$param$ctau_e1 * theoSSLT$fit$param$delta_b,
           smixvar = "Beta",
           df = theoSSLT$fit$param$nu_b)$covariance # True values

