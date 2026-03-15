############################################################
# HCM with structural equations + NEGATIVE LOGNORMAL COST
# MR, Job, Carbon ~ Normal (Carbon mean depends on Climate LV)
# Cost ~ negative lognormal
# LVs: Transcen, Conserve, Relation, Anthro, Climate
# Data: CE (4 Oct, updated).csv
############################################################

rm(list = ls())
library(apollo)
apollo_initialise()

############################################################
# 1. Control and data
############################################################

setwd("C:/Users/A02369659/Documents/12-2020/Marine Project/Kreg_hybrid/Fw_ R _ Apollo resources")

apollo_control <- list(
  modelName  = "HCM_logNormalCost_CE",
  modelDescr = "Hybrid choice model with LVs and lognormal cost",
  indivID    = "ID",
  panelData  = TRUE,
  nCores     = 10,
  weights    = "Wt_CE",
  mixing     = TRUE    # random coefficients
)

# Read data
database <- read.csv("CE (4 Oct, updated).csv", header = TRUE)

# Exclusion rule
database <- subset(database, Exclall3 == 0)

# OPTIONAL: drop rows with NA in key variables (uncomment if needed)
# vars_to_keep <- c("ID","Choice","Wt_CE",
#                   "Reserve1","Reserve2",
#                   "Jobs1","Jobs2",
#                   "Carbon1","Carbon2",
#                   "Cost1","Cost2",
#                   "Clim_1","Clim_2","Clim_5","Clim_6",
#                   "Rel_1","Rel_2","Rel_3",
#                   "TR_1","TR_2","TR_3","TR_4",
#                   "Schw_1","Schw_2","Schw_10",
#                   "Schw_4","Schw_5","Schw_9")
# database <- database[complete.cases(database[, vars_to_keep]), ]

# Order by ID and create task index
database <- database[order(database$ID), ]
database$task <- ave(database$ID, database$ID, FUN = seq_along)

cat("Rows in database:", nrow(database), "\n")

############################################################
# 2. Parameters
############################################################

# MR, Job, Carbon ~ Normal:
#   MR_i   = mu_MR   + sigma_MR   * draws_MR
#   Job_i  = mu_Job  + sigma_Job  * draws_Job
#   Carb_i = mu_Carb + gamma_carbon_climate * LV_Climate_i
#                     + sigma_Carb * draws_Carb
#
# Cost ~ NEGATIVE LOGNORMAL:
#   Cost_i = -exp( log_mu_Cost + exp(log_sigma_Cost)*draws_Cost )

apollo_beta <- c(
  ## Choice model (3 alternatives: One, Two, Neither)
  ASCSQ                = 0,      # ASC of SQ alternative (Neither)
  
  mu_MR                = 0,
  sigma_MR             = 0,
  mu_Job               = 0,
  sigma_Job            = 0,
  
  mu_Carb              = 0,      # γ10 in Eq (9)
  gamma_carbon_climate = 0,      # γ11 (effect of Climate LV on Carbon coeff)
  sigma_Carb           = 0,      # residual heterogeneity in Carbon slope
  
  log_mu_Cost          = -0.5,   # mean of underlying normal for cost
  log_sigma_Cost       = -1.0,   # log(sd) of underlying normal for cost
  
  ## Structural LV relations (Between)
  # Climate on Relation & Anthro
  beta_clim_rel        = 0,      # b1
  beta_clim_anth       = 0,      # b2
  
  # Anthro on Transcen & Conserve
  beta_anth_trans      = 0,      # b3
  beta_anth_cons       = 0,      # b4
  
  # Relation on Transcen & Conserve
  beta_rel_trans       = 0,      # b5
  beta_rel_cons        = 0,      # b6
  
  ## Measurement model – intercepts for indicators
  # Climate by Clim_1 Clim_2 Clim_5 Clim_6
  int_Clim_1           = 0,
  int_Clim_2           = 0,
  int_Clim_5           = 0,
  int_Clim_6           = 0,
  
  # Relation by Rel_1 Rel_2 Rel_3
  int_Rel_1            = 0,
  int_Rel_2            = 0,
  int_Rel_3            = 0,
  
  # Anthro by TR_1 TR_2 TR_3 TR_4
  int_TR_1             = 0,
  int_TR_2             = 0,
  int_TR_3             = 0,
  int_TR_4             = 0,
  
  # Transcen by Schw_1 Schw_2 Schw_10
  int_Schw_1           = 0,
  int_Schw_2           = 0,
  int_Schw_10          = 0,
  
  # Conserve by Schw_4 Schw_5 Schw_9
  int_Schw_4           = 0,
  int_Schw_5           = 0,
  int_Schw_9           = 0,
  
  ## Measurement model – loadings (first loading per LV fixed to 1). Coeff.” (factor loadings) in Table A-1
  load_Clim_2          = 1,
  load_Clim_5          = 1,
  load_Clim_6          = 1,
  
  load_Rel_2           = 1,
  load_Rel_3           = 1,
  
  load_TR_2            = 1,
  load_TR_3            = 1,
  load_TR_4            = 1,
  
  load_Schw_2          = 1,
  load_Schw_10         = 1,
  
  load_Schw_5          = 1,
  load_Schw_9          = 1,
  
  ## Measurement error standard deviations (log-scale, one per LV). log of the residual SD of the indicators for that LV. Measurement model of LV. 
  log_sigma_Climate    = 0,
  log_sigma_Relation   = 0,
  log_sigma_Anthro     = 0,
  log_sigma_Transcen   = 0,
  log_sigma_Conserve   = 0
  
)

apollo_fixed <- c("gamma_carbon_climate")

############################################################
# 3. Draws
############################################################

apollo_draws <- list(
  interDrawsType = "sobol",
  interNDraws    = 5000,
  interUnifDraws = c(),
  interNormDraws = c(
    "draws_MR", "draws_Job", "draws_Carb", "draws_Cost",
    "eta_climate", "eta_relation", "eta_anthro", "eta_trans", "eta_conserve"
  ),
  intraDrawsType = "sobol",
  intraNDraws    = 0,
  intraUnifDraws = c(),
  intraNormDraws = c()
)

############################################################
# 4. Random coefficients & latent variables
############################################################

apollo_randCoeff <- function(apollo_beta, apollo_inputs){
  b <- as.list(apollo_beta)
  d <- apollo_inputs$draws
  
  randcoeff <- list()
  
  ## Bottom-level LVs: Transcen and Conserve → can simplify as fixed.
  ##In Mplus, they estimate variances v1–v5 for Climate, Relation, Anthro, Transcen, Conserve. In Apollo code, fix the variance of the error term of the structural equation (the part of each LV that cannot be explained by the other LVs) of each LV disturbance to 1 via d$eta_* ~ N(0,1) → this is the structural error.Here each eta_* draw has variance 1 by construction, so fixing the innovation variances to 1 and letting the measurement side carry the scale.
 # and no scale parameter on those disturbances.That’s fine – it just means the scale of each LV is partly set by these innovations plus the measurement parameters, instead of by free v1–v5.
  ##Mplus explicitly has Conserve with Transcen (c) (a free covariance). In Apollo code, d$eta_trans and d$eta_conserve are independent draws, so corr(Transcen, Conserve) = 0 by construction.
  
  #Self-transcendence and Conservation are exogenous latent factors (ξ’s).They are not predicted by other LVs, but they are still random variables: the model estimates their variances and covariance (v4, v5, and c in the Mplus syntax). Climate, Relation, Anthro are endogenous latent factors (η’s), with structural equations and innovation variances (w1–w3).The phrase “can simplify as fixed” in the comment should really be read as “no structural equation; only measurement model and variance” — not “deterministic, same value for everyone”.
  
 # For exogenous LVs (Self-transcendence, Conservation), there is no structural equation, so they don’t have a “linear combination of other LVs” part – they are just random factors with some variance. So in that sense, for exogenous LVs: The LV equals its disturbance term, because there’s no predictor part.
  
  # This line is specifying the prior / structural distribution of the latent factor→ standard normal draw.
  
  #Self-transcendence and Conservation are exogenous LVs
  LV_Transcen <- d$eta_trans
  LV_Conserve <- d$eta_conserve
  
  ## Relation LV: on Transcen and Conserve
  LV_Relation <-
    b$beta_rel_trans * LV_Transcen +
    b$beta_rel_cons  * LV_Conserve +
    d$eta_relation
  
  ## Anthro LV: on Transcen and Conserve
  LV_Anthro <-
    b$beta_anth_trans * LV_Transcen +
    b$beta_anth_cons  * LV_Conserve +
    d$eta_anthro
  
  ## Climate LV: on Relation and Anthro
  LV_Climate <-
    b$beta_clim_rel  * LV_Relation +
    b$beta_clim_anth * LV_Anthro +
    d$eta_climate
  
  randcoeff$LV_Transcen <- LV_Transcen
  randcoeff$LV_Conserve <- LV_Conserve
  randcoeff$LV_Relation <- LV_Relation
  randcoeff$LV_Anthro   <- LV_Anthro
  randcoeff$LV_Climate  <- LV_Climate
  
  ## Preference parameters
  
  # MR and Job ~ Normal
  randcoeff$MR  <- b$mu_MR  + b$sigma_MR  * d$draws_MR
  randcoeff$Job <- b$mu_Job + b$sigma_Job * d$draws_Job
  
  # Carbon ~ Normal with mean depending on Climate LV (Eq. 9 style)
  randcoeff$Carb <-
    b$mu_Carb +
    b$gamma_carbon_climate * LV_Climate +
    b$sigma_Carb * d$draws_Carb
  
  # Cost ~ NEGATIVE LOGNORMAL
  sd_cost   <- exp(b$log_sigma_Cost)
  ln_beta_c <- b$log_mu_Cost + sd_cost * d$draws_Cost
  randcoeff$Cost <- -exp(ln_beta_c)
  
  return(randcoeff)
}

apollo_inputs <- apollo_validateInputs()

############################################################
# 5. Likelihood function
############################################################

apollo_probabilities <- function(apollo_beta, apollo_inputs,
                                 functionality = "estimate"){
  
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  P <- list()
  
  ##########################################################
  # 5.1 Choice model (MNL component)
  ##########################################################
  
  V <- list()
  
  # alt 1 & 2 have attributes; alt 3 (Neither) has ASC only
  V[["One"]] <-
    MR   * Reserve1 +
    Job  * Jobs1    +
    Carb * Carbon1  +
    Cost * Cost1
  
  V[["Two"]] <-
    MR   * Reserve2 +
    Job  * Jobs2    +
    Carb * Carbon2  +
    Cost * Cost2
  
  V[["Neither"]] <- ASCSQ  # SQ, attributes set to 0
  
  mnl_settings <- list(
    alternatives = c(One = 1, Two = 2, Neither = 3),
    avail        = 1,
    choiceVar    = Choice,
    utilities    = V
  )
  
  P[["choice"]] <- apollo_mnl(mnl_settings, functionality)
  
  ##########################################################
  # 5.2 Measurement model (continuous indicators, Normal)
  ##########################################################
  
  sigma_Climate  <- exp(log_sigma_Climate)
  sigma_Relation <- exp(log_sigma_Relation)
  sigma_Anthro   <- exp(log_sigma_Anthro)
  sigma_Transcen <- exp(log_sigma_Transcen)
  sigma_Conserve <- exp(log_sigma_Conserve)
  
  # Use only one row per individual as “anchor” for indicators
  base_row <- (task == 1)
  
  ## -------- Climate indicators: Clim_1, Clim_2, Clim_5, Clim_6 --------
  
  ND_Clim_1 <- list(
    outcomeNormal = Clim_1,
    xNormal       = LV_Climate,      # first loading = 1
    mu            = int_Clim_1,
    sigma         = sigma_Climate,
    rows          = base_row & !is.na(Clim_1)
  )
  P[["indic_Clim_1"]] <- apollo_normalDensity(ND_Clim_1, functionality)
  
  ND_Clim_2 <- list(
    outcomeNormal = Clim_2,
    xNormal       = load_Clim_2 * LV_Climate,
    mu            = int_Clim_2,
    sigma         = sigma_Climate,
    rows          = base_row & !is.na(Clim_2)
  )
  P[["indic_Clim_2"]] <- apollo_normalDensity(ND_Clim_2, functionality)
  
  ND_Clim_5 <- list(
    outcomeNormal = Clim_5,
    xNormal       = load_Clim_5 * LV_Climate,
    mu            = int_Clim_5,
    sigma         = sigma_Climate,
    rows          = base_row & !is.na(Clim_5)
  )
  P[["indic_Clim_5"]] <- apollo_normalDensity(ND_Clim_5, functionality)
  
  ND_Clim_6 <- list(
    outcomeNormal = Clim_6,
    xNormal       = load_Clim_6 * LV_Climate,
    mu            = int_Clim_6,
    sigma         = sigma_Climate,
    rows          = base_row & !is.na(Clim_6)
  )
  P[["indic_Clim_6"]] <- apollo_normalDensity(ND_Clim_6, functionality)
  
  ## -------- Relation indicators: Rel_1 Rel_2 Rel_3 --------
  
  ND_Rel_1 <- list(
    outcomeNormal = Rel_1,
    xNormal       = LV_Relation,
    mu            = int_Rel_1,
    sigma         = sigma_Relation,
    rows          = base_row & !is.na(Rel_1)
  )
  P[["indic_Rel_1"]] <- apollo_normalDensity(ND_Rel_1, functionality)
  
  ND_Rel_2 <- list(
    outcomeNormal = Rel_2,
    xNormal       = load_Rel_2 * LV_Relation,
    mu            = int_Rel_2,
    sigma         = sigma_Relation,
    rows          = base_row & !is.na(Rel_2)
  )
  P[["indic_Rel_2"]] <- apollo_normalDensity(ND_Rel_2, functionality)
  
  ND_Rel_3 <- list(
    outcomeNormal = Rel_3,
    xNormal       = load_Rel_3 * LV_Relation,
    mu            = int_Rel_3,
    sigma         = sigma_Relation,
    rows          = base_row & !is.na(Rel_3)
  )
  P[["indic_Rel_3"]] <- apollo_normalDensity(ND_Rel_3, functionality)
  
  ## -------- Anthro indicators: TR_1 TR_2 TR_3 TR_4 --------
  
  ND_TR_1 <- list(
    outcomeNormal = TR_1,
    xNormal       = LV_Anthro,
    mu            = int_TR_1,
    sigma         = sigma_Anthro,
    rows          = base_row & !is.na(TR_1)
  )
  P[["indic_TR_1"]] <- apollo_normalDensity(ND_TR_1, functionality)
  
  ND_TR_2 <- list(
    outcomeNormal = TR_2,
    xNormal       = load_TR_2 * LV_Anthro,
    mu            = int_TR_2,
    sigma         = sigma_Anthro,
    rows          = base_row & !is.na(TR_2)
  )
  P[["indic_TR_2"]] <- apollo_normalDensity(ND_TR_2, functionality)
  
  ND_TR_3 <- list(
    outcomeNormal = TR_3,
    xNormal       = load_TR_3 * LV_Anthro,
    mu            = int_TR_3,
    sigma         = sigma_Anthro,
    rows          = base_row & !is.na(TR_3)
  )
  P[["indic_TR_3"]] <- apollo_normalDensity(ND_TR_3, functionality)
  
  ND_TR_4 <- list(
    outcomeNormal = TR_4,
    xNormal       = load_TR_4 * LV_Anthro,
    mu            = int_TR_4,
    sigma         = sigma_Anthro,
    rows          = base_row & !is.na(TR_4)
  )
  P[["indic_TR_4"]] <- apollo_normalDensity(ND_TR_4, functionality)
  
  ## -------- Transcen indicators: Schw_1 Schw_2 Schw_10 --------
  
  ND_Schw_1 <- list(
    outcomeNormal = Schw_1, #observed likert response
    xNormal       = LV_Transcen, # the unobserved factor
    mu            = int_Schw_1,
    sigma         = sigma_Transcen,
    rows          = base_row & !is.na(Schw_1)
  )
  P[["indic_Schw_1"]] <- apollo_normalDensity(ND_Schw_1, functionality)
  
  ND_Schw_2 <- list(
    outcomeNormal = Schw_2,
    xNormal       = load_Schw_2 * LV_Transcen,
    mu            = int_Schw_2,
    sigma         = sigma_Transcen,
    rows          = base_row & !is.na(Schw_2)
  )
  P[["indic_Schw_2"]] <- apollo_normalDensity(ND_Schw_2, functionality)
  
  ND_Schw_10 <- list(
    outcomeNormal = Schw_10,
    xNormal       = load_Schw_10 * LV_Transcen,
    mu            = int_Schw_10,
    sigma         = sigma_Transcen,
    rows          = base_row & !is.na(Schw_10)
  )
  P[["indic_Schw_10"]] <- apollo_normalDensity(ND_Schw_10, functionality)
  
  ## -------- Conserve indicators: Schw_4 Schw_5 Schw_9 --------
  
  ND_Schw_4 <- list(
    outcomeNormal = Schw_4,
    xNormal       = LV_Conserve,
    mu            = int_Schw_4,
    sigma         = sigma_Conserve,
    rows          = base_row & !is.na(Schw_4)
  )
  P[["indic_Schw_4"]] <- apollo_normalDensity(ND_Schw_4, functionality)
  
  ND_Schw_5 <- list(
    outcomeNormal = Schw_5,
    xNormal       = load_Schw_5 * LV_Conserve,
    mu            = int_Schw_5,
    sigma         = sigma_Conserve,
    rows          = base_row & !is.na(Schw_5)
  )
  P[["indic_Schw_5"]] <- apollo_normalDensity(ND_Schw_5, functionality)
  
  ND_Schw_9 <- list(
    outcomeNormal = Schw_9,
    xNormal       = load_Schw_9 * LV_Conserve,
    mu            = int_Schw_9,
    sigma         = sigma_Conserve,
    rows          = base_row & !is.na(Schw_9)
  )
  P[["indic_Schw_9"]] <- apollo_normalDensity(ND_Schw_9, functionality)
  
  ##########################################################
  # 5.3 Combine model components + panel, draws, weights
  ##########################################################
  
  P <- apollo_combineModels(P, apollo_inputs, functionality)
  P <- apollo_panelProd(P, apollo_inputs, functionality)
  P <- apollo_weighting(P, apollo_inputs, functionality)
  P <- apollo_avgInterDraws(P, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  
  return(P)
}

############################################################
# 6. Estimation
############################################################

apollo_inputs <- apollo_validateInputs()

estimate_settings <- list(
  estimationRoutine = "BFGS",
  maxIterations     = 400,
  writeIter         = TRUE
)

model <- apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs,
  estimate_settings = estimate_settings
)


apollo_modelOutput(model)
apollo_saveOutput(model)



### Cost conversion
############################################################
# Moments of NEGATIVE LOGNORMAL cost coefficient
#   from Apollo HCM model object `model`
############################################################

library(mvtnorm)

est <- model$estimate
V   <- model$varcov

# Parameters on log scale
par_names <- c("log_mu_Cost", "log_sigma_Cost")
theta_hat <- est[par_names]
V_theta   <- V[par_names, par_names]

mu_c <- theta_hat["log_mu_Cost"]
ls_c <- theta_hat["log_sigma_Cost"]

## ---------- 1) Point estimates (analytic) -----------------

# SD of underlying normal
sd_under <- exp(ls_c)
sd2      <- sd_under^2

# MEDIAN of β_cost (negative lognormal)
beta_cost_med  <- -exp(mu_c)

# MEAN of β_cost (negative lognormal)
beta_cost_mean <- -exp(mu_c + 0.5*sd2)

# VARIANCE and SD of β_cost (same for negative and positive)
var_pos      <- (exp(sd2) - 1) * exp(2*mu_c + sd2)
beta_cost_sd <- sqrt(var_pos)

#beta_cost_med → report median!
beta_cost_mean
beta_cost_sd

## ---------- 2) Standard errors via simulation --------------

set.seed(123)
nsim <- 10000

theta_sim <- rmvnorm(nsim, mean = theta_hat, sigma = V_theta)
mu_sim    <- theta_sim[,1]
ls_sim    <- theta_sim[,2]

sd_under_sim <- exp(ls_sim)
sd2_sim      <- sd_under_sim^2

# Median, mean and SD in each simulation
beta_cost_med_sim  <- -exp(mu_sim)
beta_cost_mean_sim <- -exp(mu_sim + 0.5*sd2_sim)
var_pos_sim        <- (exp(sd2_sim) - 1) * exp(2*mu_sim + sd2_sim)
beta_cost_sd_sim   <- sqrt(var_pos_sim)

# Point estimates (simulation-based, should match analytic closely)
beta_cost_med_hat  <- mean(beta_cost_med_sim)
beta_cost_mean_hat <- mean(beta_cost_mean_sim)
beta_cost_sd_hat   <- mean(beta_cost_sd_sim)

# Standard errors (what report in the table)
se_beta_cost_med   <- sd(beta_cost_med_sim)
se_beta_cost_mean  <- sd(beta_cost_mean_sim)
se_beta_cost_sd    <- sd(beta_cost_sd_sim)

# Nice summary table for paper
cost_moments <- cbind(
  estimate = c(beta_cost_med_hat,  beta_cost_mean_hat,  beta_cost_sd_hat),
  se       = c(se_beta_cost_med,   se_beta_cost_mean,   se_beta_cost_sd)
)
rownames(cost_moments) <- c("median_beta_cost", "mean_beta_cost", "sd_beta_cost")

round(cost_moments, 5)


