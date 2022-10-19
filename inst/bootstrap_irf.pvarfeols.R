library(panelvar)

# Example:
ex1_pvarfeols <-
  pvarfeols(dependent_vars = c("log_sales", "log_price"),
            lags = 1,
            exog_vars = c("cpi"),
            transformation = "demean",
            data = Cigar,
            panel_identifier= c("state", "year"))

summary(ex1_pvarfeols)

ex1_pvarfeols_again <-
  pvarfeols(dependent_vars = c("log_sales", "log_price"),
            lags = 1,
            exog_vars = c("cpi"),
            transformation = "demean",
            data = ex1_pvarfeols$Set_Vars_with_NAs,
            panel_identifier= c("category", "period"))

summary(ex1_pvarfeols_again)


ex1_pvarfeols_boot <- bootstrap_irf(model = ex1_pvarfeols,
                                          typeof_irf = c("GIRF"),
                                          n.ahead = 8,
                                          nof_Nstar_draws = 100,
                                          confidence.band = 0.95,
                                    mc.cores = 100)

plot(girf(ex1_pvarfeols, n.ahead = 8, 10), ex1_pvarfeols_boot)




ex1_dahlberg_data_feols <-
  pvarfeols(dependent_vars = c("expenditures", "revenues", "grants"),
            lags = 2,
            transformation = "demean",
            data = Dahlberg,
            panel_identifier= c("id", "year"))

summary(ex1_dahlberg_data_feols)

ex1_dahlberg_data_feols_bs <-  bootstrap_irf(ex1_dahlberg_data_feols, typeof_irf = c("GIRF"),
                                             n.ahead = 8,
                                             nof_Nstar_draws = 100,
                                             confidence.band = 0.99)

plot(girf(ex1_dahlberg_data_feols, n.ahead = 8, 8), ex1_dahlberg_data_feols_bs)


ex1_pvargmm <- pvargmm(dependent_vars = c("log_sales", "log_price"),
                          lags = 1,
                          #predet_vars = c("log_ndi"),
                          exog_vars = c("cpi"),
                          transformation = "fod",
                          data = Cigar,
                          panel_identifier= c("state", "year"),
                          steps = c("twostep"),
                          system_instruments = TRUE,
                          max_instr_dependent_vars = 10,
                          max_instr_predet_vars = 10,
                          min_instr_dependent_vars = 2L,
                          min_instr_predet_vars = 1L,
                          collapse = TRUE
)

summary(ex1_pvargmm)

ex1_pvargmm_again <- pvargmm(dependent_vars = c("log_sales", "log_price"),
                       lags = 1,
                       #predet_vars = c("log_ndi"),
                       exog_vars = c("cpi"),
                       transformation = "fod",
                       data = ex1_pvargmm$Set_Vars_with_NAs,
                       panel_identifier= c("category", "period"),
                       steps = c("twostep"),
                       system_instruments = TRUE,
                       max_instr_dependent_vars = 10,
                       max_instr_predet_vars = 10,
                       min_instr_dependent_vars = 2L,
                       min_instr_predet_vars = 1L,
                       collapse = TRUE
)

summary(ex1_pvargmm_again)

ex1_pvargmm_boot <- bootstrap_irf(model = ex1_pvargmm,
                                    typeof_irf = c("GIRF"),
                                    n.ahead = 8,
                                    nof_Nstar_draws = 500,
                                    confidence.band = 0.95,
                                  mc.cores = 200)

plot(girf(ex1_pvargmm, n.ahead = 8, 10), ex1_pvargmm_boot)


bootstrap_irf.pvarfeols <- function(model,
                                    typeof_irf = c("OIRF", "GIRF"),
                                    n.ahead,
                                    nof_Nstar_draws,
                                    confidence.band = 0.95
){
  
  #example:
  # model <- ex1_feols
  # data  <- Cigar
  # typeof_irf = c("OIRF")
  # nof_Nstar_draws <- 1
  # n.ahead <- 8
  # confidence.band <- 0.95
  
  
  # Start First Block: Set up parameters of the model
  
  # We use blockwise sampling over individuals. In our GMM frameowork, sampling over time periods makes no sense
  # as it would destroy the GMM instruments.
  # Input arguments
  # 1. Number of N* draws: nof_Nstar_draws
  # 2. Size of each Nstar_draws:
  # 3. Each Nstar_draws: are a random resample with replacement of blocks
  
  # Fix the input argument as a lot of models will be estimated.
  
  Data_initial <- model$Set_Vars
  
  # optionales Argument
  # ???
  size_Nstar_draws <- length(unique(model$Set_Vars$category))
  
  # Nstar_draws:
  Nstar_draws <- list()
  
  # Set up an uniform distribution and sample individuals between 1 and N
  for (i0 in 1:nof_Nstar_draws){
    
    possible.draws <- length(unique(model$Set_Vars[,c("category")]))
    
    x <- unique(model$Set_Vars[,c("category")])
    
    prob.x <- rep(1 / possible.draws, possible.draws)
    Nstar_draws[[i0]] <- sample(x, size = size_Nstar_draws, replace = TRUE, prob = prob.x)
    
  }
  
  # Set up the model reestimation:
  pb <-
    progress::progress_bar$new(format = "Bootstrapping: [:bar] Iteration :current of :total" ,
                               clear = TRUE,
                               total = nof_Nstar_draws)
  pvar_irf            <- list()
  for (i0 in 1:nof_Nstar_draws){
    
    pb$tick()
    # Set up the sampled data for each nof_Nstar_draws estimations
    #data_resampled <- Data_initial[Data_initial[, model$panel_identifier[1]] ==  Nstar_draws[[i0]][1],]
    data_resampled <- Data_initial[Data_initial[, "category"] ==  Nstar_draws[[i0]][1],]
    
    # Overwrite data_resampled-id, requires unique category
    data_resampled[,"category"] <- 1
    
    
    # combine resampled data and rename category
    for (i1 in 2:size_Nstar_draws){
      
      zwischen <- Data_initial[Data_initial[, "category"] ==  Nstar_draws[[i0]][i1],]
      # Overwrite data_resampled-id, requires unique category
      zwischen[,"category"] <- i1
      
      data_resampled <- rbind(data_resampled, zwischen)
      
    }
    
    # Estimate the pvar model with the sampled data and only save the errors.
    # do.call(rbind,model$residuals)
    fct_agr_list <- list(
      dependent_vars = model$dependent_vars,
      lags = model$lags,
      exog_vars = model$exog_vars,
      transformation = c("demean"),
      data = data_resampled,
      panel_identifier = c("category", "period"))
    
    if (is.null(model$exog_vars) == TRUE){fct_agr_list$exog_vars = model$exog_vars}
    
    pvar_zwischen <- suppressWarnings(do.call(pvarfeols, fct_agr_list))
    
    # Give track how the parameters of the estimated model change.
    coef(pvar_zwischen)
    # Save the residuals:
    #Residuals_resampled[[i0]] <- do.call(rbind,pvar_zwischen$residuals)
    
    if (typeof_irf == c("OIRF")){
      pvar_irf[[i0]] <- oirf(model = pvar_zwischen,
                             n.ahead = n.ahead)
    }
    
    if (typeof_irf == c("GIRF")){
      pvar_irf[[i0]] <- girf(model = pvar_zwischen,
                             n.ahead = n.ahead,
                             ma_approx_steps = n.ahead)
    }
    
  }
  
  # Old part from vars package:
  # definite a mistake here:
  #lower <- ci / 2
  #upper <- 1 - ci / 2
  
  
  two_sided_bound <- 1-confidence.band
  
  lower <- two_sided_bound / 2
  upper <- 1 - lower
  
  
  mat.l <- matrix(NA, nrow = n.ahead, ncol = ncol(pvar_irf[[1]][[1]]))
  mat.u <- matrix(NA, nrow = n.ahead, ncol = ncol(pvar_irf[[1]][[1]]))
  
  Lower <- list()
  Upper <- list()
  
  # Impuls: Number of endogenous variables:
  idx1 <- length(pvar_irf[[1]])
  
  # Response: For each Shock all endogenous Variables of the system are shocked.
  idx2 <- ncol(pvar_irf[[1]][[1]])
  
  idx3 <- n.ahead
  
  temp <- rep(NA, length(pvar_irf))
  
  for(j in 1 : idx1){
    for(m in 1 : idx2){
      for(l in 1 : idx3){
        for(i in 1 : nof_Nstar_draws){
          if(idx2 > 1){
            temp[i] <- pvar_irf[[i]][[j]][l, m]
          } else {
            temp[i] <- matrix(pvar_irf[[i]][[j]])[l, m]
          }
        }
        mat.l[l, m] <- quantile(temp, lower, na.rm = TRUE)
        mat.u[l, m] <- quantile(temp, upper, na.rm = TRUE)
      }
    }
    colnames(mat.l) <- model$dependent_vars
    colnames(mat.u) <- model$dependent_vars
    Lower[[j]] <- mat.l
    Upper[[j]] <- mat.u
  }
  
  names(Lower) <- model$dependent_vars
  names(Upper) <- model$dependent_vars
  return(list(Lower = Lower, Upper = Upper, Nstar_draws = Nstar_draws, CI = confidence.band))
}
