# AR(1) Represenatation of a PVAR AR model
# @param Phi Phi
# @param lags lags

pvar1_phi_representation <- function(Phi, lags){

  # ensure that A_endo is a matrix (important for single equation with multiple lags case)
  A_endo <- matrix(Phi[,c(1:( nrow(Phi) * lags ))], ncol = nrow(Phi) * lags)
  Dim_I_k <- nrow(A_endo)

  nof_equations <- dim(Phi)[1]
  I_k <- diag(Dim_I_k)

  if (lags > 1){

    I_0 <- matrix(0, nrow = Dim_I_k, ncol = Dim_I_k)

    J <- cbind(diag(Dim_I_k), matrix(0, nrow = Dim_I_k, ncol = Dim_I_k*(lags-1)))

    lower_part <- cbind(diag(Dim_I_k*(lags-1)), matrix(0, nrow = Dim_I_k*(lags-1), ncol = Dim_I_k))
  }
  if (lags == 1){
    lower_part <- NULL

    J <- diag(Dim_I_k)

  }

  A_full <- rbind(A_endo, lower_part)

  return(A_full)

}

# MA(ma_approx_steps) Approximation of a PVAR AR model
# @param Phi Phi
# @param ma_approx_steps ma_approx_steps
# @param lags lags

ma_phi_representation <- function(Phi, ma_approx_steps, lags){

  # Coefficients of endogenous variables:
  A_endo <- Phi[,c(1:( nrow(Phi) * lags ))]

  Dim_I_k <- nrow(A_endo)

  I_k <- diag(Dim_I_k)

  if (lags > 1){
    I_0 <- matrix(0, nrow = Dim_I_k, ncol = Dim_I_k)

    J <- cbind(diag(Dim_I_k), matrix(0, nrow = Dim_I_k, ncol = Dim_I_k*(lags-1)))

    lower_part <- cbind(diag(Dim_I_k*(lags-1)), matrix(0, nrow = Dim_I_k*(lags-1), ncol = Dim_I_k))
  }
  if (lags == 1){
    lower_part <- NULL
    J <- diag(Dim_I_k)

  }

  A_full <- rbind(A_endo, lower_part)

  MA_Phi <- list()

  MA_Phi[[1]] <- I_k


  # Generate MA_Phi(1) to MA_Phi(ma_approx_steps)
  for (i in 2:(ma_approx_steps)){

    MA_Phi[[i]] <- J %*% matrixcalc::matrix.power(as.matrix(A_full),(i-1)) %*% t(J)


  }

  return(MA_Phi)

}

#' Empirical estimation of PVAR Impulse Response Confidence Bands
#' 
#' @param model A PVAR model
#' @param typeof_irf \code{"OIRF"} or \code{GIRF}
#' @param n.ahead n ahead steps
#' @param nof_Nstar_draws Number of draws
#' @param confidence.band Confidence band
#' @export
#' @description 
#' Uses blockwise sampling of individuals (bootstrapping)
#' @examples 
#' \dontrun{
#' data("ex1_dahlbergdata")
#' bootstrap_irf(ex1_dahlberg_data, 
#'               typeof_irf = "OIRF", 
#'               n.ahead = 12, 
#'               nof_Nstar_draws = 2, 
#'               confidence.band = 0.95)
#' }

bootstrap_irf <- function(model, typeof_irf, n.ahead, nof_Nstar_draws, confidence.band) UseMethod("bootstrap_irf")

#' @rdname bootstrap_irf
#' @export
bootstrap_irf.pvargmm <- function(model,
                                  typeof_irf = c("OIRF", "GIRF"),
                                    n.ahead,
                                    nof_Nstar_draws,
                                    confidence.band = 0.95
                                    ){

  #example:
  #model <- test_Q_for_endo_ex1
  #model <- test_Q_for_endo_ex2B
  #steps = c("onestep")
  #data  <- EmplUK
  #nof_Nstar_draws <- 20
  #n.ahead <- 12
  #convidence.band <- 0.95


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
  size_Nstar_draws <- length(model$instruments)

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
        predet_vars = model$predet_vars,
        exog_vars = model$exog_vars,
        transformation = model$transformation,
        data = data_resampled,
        panel_identifier = c("category", "period"),
        steps = model$steps,
        system_instruments = model$system_instruments,
        max_instr_dependent_vars = model$max_instr_dependent_vars,
        max_instr_predet_vars = model$max_instr_predet_vars,
        min_instr_dependent_vars = model$min_instr_dependent_vars,
        min_instr_predet_vars = model$min_instr_predet_vars,
        collapse = model$collapse,
        system_constant = model$system_constant,
        tol = model$tol,
        progressbar = FALSE
    )

    if (is.null(model$predet_vars) == TRUE){fct_agr_list$predet_vars = model$predet_vars}
    if (is.null(model$exog_vars) == TRUE){fct_agr_list$exog_vars = model$exog_vars}
    if (is.null(model$tol) == TRUE){fct_agr_list$tol = model$tol}


    pvar_zwischen <- suppressWarnings(do.call(pvargmm, fct_agr_list))

    # Give track how the parameters of the estimated model change.
    pvar_zwischen$second_step
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


#' Andrews Lu MMSC Criteria based on Hansen-J-Statistic
#'
#' @param model A PVAR model
#' @param HQ_criterion Hannan Quinn criterion
#' @export
#' @description 
#' ...
#' @return
#' BIC, AIC and HQIC
#' @references 
#' 
#' Andrews, D., Lu, B. (2001) Consistent Model and Momement Selection Procedures for GMM Estimation with Application to Dynamic Panel Data Models, \emph{Journal of Econometrics}, \bold{101}(1), 123--164, \doi{10.1016/S0304-4076(00)00077-4}
#' @examples 
#' data("ex3_abdata")
#' Andrews_Lu_MMSC(ex3_abdata)
#' 
#' 
#' 
Andrews_Lu_MMSC <- function(model, HQ_criterion = 2.1) UseMethod("Andrews_Lu_MMSC")

#' @rdname Andrews_Lu_MMSC
#' @export 
Andrews_Lu_MMSC.pvargmm <- function(model, HQ_criterion = 2.1){

  if (HQ_criterion < 2)
    stop("Hannan Quinn criterion is not defined for less than 2!")

  # Andrews-Lu-2001 notation:
  # c is the number of moment conditions.


  # b is the number of parameters:
  if (model$steps == c("onestep")){
    b <- ncol(model$first_step)
    c <- nrow(model$weighting_matrix_first_step)
  }
  if (model$steps == c("twostep")){ b <- ncol(model$second_step)
  c <- nrow(model$weighting_matrix_second_step)}
  if (model$steps == c("mstep")){ b <- ncol(model$m_step)
  c <- nrow(model$weighting_matrix_m_step)}

  # n used number of observations. N*T but only non NA rows.
  n <- nrow(stats::na.omit(model$Set_Vars))

  # Q for the Hannan Quinn criterion: Q is set as in the Andrews Lu Paper.
  Q <- HQ_criterion

  #J_test_stat <- hansen_j_test(model = model, steps = steps)
  J_test_stat <- hansen_j_test(model = model)
  J_test_stat <- J_test_stat$statistic

  # Minimize these criteria:
  # The lower the Hansen STAT, the lower are the criteria.
  # The higher the number of parmeters the higher the criteria (i.e. number of parameters serve as penalty)
  MMSC_BIC <- J_test_stat - (c-b)*log(n)

  MMSC_AIC <- J_test_stat - (c-b)*2

  MMSC_HQIC <- J_test_stat - Q*(c-b)*log(log(n))

  MMSC <- list(as.numeric(MMSC_BIC), as.numeric(MMSC_AIC), as.numeric(MMSC_HQIC))

  names(MMSC) <- c("MMSC_BIC", "MMSC_AIC", "MMSC_HQIC")

  return(MMSC)

}


# Forcast Error Variance Decomposition for GIRF:
# under construction, not tested yet
#
# @param model A PVAR model
# @param steps steps
# @param n.ahead n ahead
# @param lags lags
# @param normalize normalize
# @export
#
# fevd_generalised <- function(model, steps, n.ahead=10, lags, normalize=TRUE){
# 
#   if (steps == c("onestep")){
# 
#     Phi <- model$first_step
# 
#   }
# 
#   if (steps == c("twostep")){
# 
#     Phi<-model$second_step
# 
#   }
# 
#   MA_Phi <- ma_phi_representation(Phi = Phi,
#                                   ma_approx_steps = n.ahead,
#                                   lags = lags)
# 
#   Sigma_Hat1 <- vcov(model)
# 
#   sigmas <- sqrt(diag(Sigma_Hat1))
# 
#   # Re-structure the code at this position:
# 
#   girf_output_time <- list()
# 
#   for (i0 in 1:n.ahead){
# 
#     girf_output_time[[i0]] <-  t( MA_Phi[[i0]] %*% Sigma_Hat1)  /  sigmas
# 
#   }
#   # as for GIRF:
# 
#   # Adjust from here:
#   # d <- array(0, dim(A)[c(2,3)])
# 
#   MSE_diag <- matrix(NA, nrow = nrow( girf_output_time[[1]] ), ncol = ncol( girf_output_time[[1]]) )
# 
#   # Calcuate the MSE as in L?tkepohl page 64. first equation.
#   # Fill the matrix for each forecasting period.
#   for (i0 in 1:dim(MSE_diag)[2]){
# 
#     MSE_diag[,i0] <- diag( girf_output_time[[i0]] %*% Sigma_Hat1 %*% t(girf_output_time[[i0]]) )
# 
#   }
# 
# 
#   # Calculate the numerator of eq. 2.3.37 in L?tkephol book. page 64.
#   numerator_fevd <- Reduce('+', mapply(function(i) (girf_output_time[[i0]])^2,
#                                        1:length(girf_output_time), SIMPLIFY = FALSE ) )
# 
#   denominator_fevd <- c(apply(MSE_diag, 1, sum))
# 
# 
# 
#   fevd <- numerator_fevd / denominator_fevd
# 
#   if (normalize) {
# 
#     return(fevd/apply(fevd, 1, sum))
# 
#   }
# 
#   else {
# 
#     return(fevd)
# 
#   }
# 
# }

#' Forcast Error Variance Decomposition for PVAR
#'
#' @param model A PVAR model
#' @param n.ahead Number of steps
#' @description 
#' Computes the forecast error variance decomposition of a PVAR(p) model.
#' @return
#' A list with forecast error variances as matrices for each variable.
#' @details
#' The estimation is based on orthogonalised impulse response functions.
#' 
#' @note 
#' A \code{plot} method will be provided in future versions.
#' @references
#' Pfaff, B. (2008) VAR, SVAR and SVEC Models: Implementation Within R Package vars, \emph{Journal of Statistical Software} \bold{27}(4) \url{https://www.jstatsoft.org/v27/i04/}
#' @seealso 
#' \code{\link{pvargmm}} for model estimaion
#' 
#' \code{\link{oirf}} for orthogonal impulse response function
#' 
#' @examples 
#' data("ex1_dahlberg_data")
#' fevd_orthogonal(ex1_dahlberg_data, n.ahead = 8)
#' 
#' @export

fevd_orthogonal <- function(model, n.ahead=10) UseMethod("fevd_orthogonal")

#' @rdname fevd_orthogonal
#' @export
fevd_orthogonal.pvargmm <- function(model, n.ahead=10){

  if (model$steps == c("onestep")){
    Phi <- model$first_step
  }

  if (model$steps == c("twostep")){
    Phi<-model$second_step
  }

  if (model$steps == c("mstep")){
    Phi<-model$m_step
  }
  
  # P <- t(chol(covm))
  # MA Phi representation:
  MA_Phi <- ma_phi_representation(Phi = Phi,
                                  ma_approx_steps = n.ahead,
                                  lags = model$lags)

  # Variance Co-Variance of the model:
  Sigma_Hat1 <- vcov(model)
  
  
  # Sigma.yh results:
  
  
  Sigma.yh <- list()
  Sigma.yh[[1]] <- Sigma_Hat1
 
  if (n.ahead > 1) {
    for (i in 2:n.ahead) {
      temp <- matrix(0, nrow = length(model$dependent_vars), 
                     ncol = length(model$dependent_vars))
      for (j in 2:i) {
        temp <- temp + MA_Phi[[j]] %*% Sigma_Hat1 %*% t(MA_Phi[[j]])
      }
      Sigma.yh[[i]] <- temp + Sigma.yh[[1]]
    }
  }
  
  
  # Set up Psi == MA_Phi_P
  
  P <- chol(Sigma_Hat1)
  MA_Phi_P <- list()

  # Calculate 2.3.27 in Luetkepohl p. 58.
  MA_Phi_P[[1]] <- t(P)
  
  if (n.ahead > 1){
    for (i0 in 2:n.ahead){

      MA_Phi_P[[i0]] <- MA_Phi[[i0]] %*% t(P)

    }
  }
  
  
  # Set up Psi == MA_Phi_P
  
  
  mse <- matrix(NA, nrow = n.ahead, ncol = length(model$dependent_vars))
  Omega <- array(0, dim = c(n.ahead, length(model$dependent_vars), 
                            length(model$dependent_vars)))
  
  for(i in 1 : n.ahead){
    
    mse[i, ] <- diag(Sigma.yh[[i]])
    
    temp <- matrix(0, length(model$dependent_vars), length(model$dependent_vars))
    
    for(l in 1 : length(model$dependent_vars)){
      
      for(m in 1 : length(model$dependent_vars)){
        
        for(j in 1 : i){
          
          temp[l, m] <- temp[l, m] + (MA_Phi_P[[j]][l,m])^2
        
        }
      }
    }
    
    temp <- temp / mse[i, ]
    
    for(j in 1 : length(model$dependent_vars)){
      Omega[i, , j] <- temp[j, ]
    }
    
  }
  
  # Prepare results:
  result <- list()
  for(i in 1 : length(model$dependent_vars)){
    result[[i]] <- matrix(Omega[, , i], 
                          nrow = n.ahead, 
                          ncol = length(model$dependent_vars))
    
    colnames(result[[i]]) <- model$dependent_vars
  }
  names(result) <- model$dependent_vars
  
  return(result)
  
}

#' @rdname fevd_orthogonal
#' @export
fevd_orthogonal.pvarfeols <- function(model, n.ahead=10){
  
  Phi = coef(model)
  
  # P <- t(chol(covm))
  # MA Phi representation:
  MA_Phi <- ma_phi_representation(Phi = Phi,
                                  ma_approx_steps = n.ahead,
                                  lags = model$lags)
  
  # Variance Co-Variance of the model:
  Sigma_Hat1 <- vcov(model)
  
  
  # Sigma.yh results:
  
  
  Sigma.yh <- list()
  Sigma.yh[[1]] <- Sigma_Hat1
  
  if (n.ahead > 1) {
    for (i in 2:n.ahead) {
      temp <- matrix(0, nrow = length(model$dependent_vars), 
                     ncol = length(model$dependent_vars))
      for (j in 2:i) {
        temp <- temp + MA_Phi[[j]] %*% Sigma_Hat1 %*% t(MA_Phi[[j]])
      }
      Sigma.yh[[i]] <- temp + Sigma.yh[[1]]
    }
  }
  
  
  # Set up Psi == MA_Phi_P
  
  P <- chol(Sigma_Hat1)
  MA_Phi_P <- list()
  
  # Calculate 2.3.27 in Luetkepohl p. 58.
  MA_Phi_P[[1]] <- t(P)
  
  if (n.ahead > 1){
    for (i0 in 2:n.ahead){
      
      MA_Phi_P[[i0]] <- MA_Phi[[i0]] %*% t(P)
      
    }
  }
  
  
  # Set up Psi == MA_Phi_P
  
  
  mse <- matrix(NA, nrow = n.ahead, ncol = length(model$dependent_vars))
  Omega <- array(0, dim = c(n.ahead, length(model$dependent_vars), 
                            length(model$dependent_vars)))
  
  for(i in 1 : n.ahead){
    
    mse[i, ] <- diag(Sigma.yh[[i]])
    
    temp <- matrix(0, length(model$dependent_vars), length(model$dependent_vars))
    
    for(l in 1 : length(model$dependent_vars)){
      
      for(m in 1 : length(model$dependent_vars)){
        
        for(j in 1 : i){
          
          temp[l, m] <- temp[l, m] + (MA_Phi_P[[j]][l,m])^2
          
        }
      }
    }
    
    temp <- temp / mse[i, ]
    
    for(j in 1 : length(model$dependent_vars)){
      Omega[i, , j] <- temp[j, ]
    }
    
  }
  
  # Prepare results:
  result <- list()
  for(i in 1 : length(model$dependent_vars)){
    result[[i]] <- matrix(Omega[, , i], 
                          nrow = n.ahead, 
                          ncol = length(model$dependent_vars))
    
    colnames(result[[i]]) <- model$dependent_vars
  }
  names(result) <- model$dependent_vars
  
  return(result)
  
}

# Demeaning
# @param x vector with data
# @examples 
# demean(c(1:5))

demean <- function(x) {
  x - mean(x, na.rm = TRUE)
}

# Panel demeaning
# @param x data
# @param transformation type of transformation
panel_demean <- function(x, transformation) {
  
  res <- apply(x[, -(1:2)], 2, demean)
  colnames(res) <- paste("demeaned", colnames(res), sep = "_")
  
  rownames(res) <- NULL
  res
}

