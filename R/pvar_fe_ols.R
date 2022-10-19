#' Fixed Effects Estimator for PVAR Model
#'
#' @param dependent_vars Dependent variables
#' @param lags Number of lags of dependent variables
#' @param exog_vars Exogenous variables
#' @param transformation Demeaning \code{"demean"}
#' @param data Data set
#' @param panel_identifier Vector of panel identifiers
#'
#' @export
#' @description
#' This function estimates a stationary PVAR with fixed effects.
#' @examples
#' data(Cigar)
#' ex1_feols <-
#' pvarfeols(dependent_vars = c("log_sales", "log_price"),
#'           lags = 1,
#'           exog_vars = c("cpi"),
#'           transformation = "demean",
#'           data = Cigar,
#'           panel_identifier= c("state", "year"))
#'           
#' summary(ex1_feols)         

pvarfeols <- function(dependent_vars,
                           lags,
                           exog_vars,
                           transformation = c("demean"),
                           data,
                           panel_identifier = c(1, 2)){

  # Set up Set_Vars ---------------------------------------------------------

  # Drop factor levels from data
  data <- droplevels(data)

  # Select only required variables from data
  required_vars <- c(dependent_vars)
  if(!missing(exog_vars)) required_vars <- c(required_vars, exog_vars)

  if (is.numeric(panel_identifier)) {
    Set_Vars <-
      data[, c(colnames(data)[panel_identifier], required_vars)]
  } else {
    Set_Vars <-
      data[, c(panel_identifier, required_vars)]
  }

  nof_dependent_vars <- length(dependent_vars)
  if(missing(exog_vars)) {
    nof_exog_vars <- 0
  } else {
    nof_exog_vars <- length(exog_vars)
  }

  # sort categories and periods
  categories <- sort(unique(Set_Vars[, 1]))
  periods <- sort(unique(Set_Vars[, 2]))

  nof_categories <- length(categories)
  nof_periods <- length(periods)


  # keep names of panel identifiers
  name_category <- names(Set_Vars)[1]
  name_period <- names(Set_Vars)[2]
  names(Set_Vars)[1:2] <- c("category", "period")

  Set_Vars <- Set_Vars[order(Set_Vars$category, Set_Vars$period,decreasing = FALSE),]
  Set_Vars$category <- factor(Set_Vars$category)
  Set_Vars$period <- factor(Set_Vars$period)

  # Descriptive statistics --------------------------------------------------

  nof_observations <- dim(data)[1]
  data_panel_identifier <- data[panel_identifier]

  obs_per_group_avg <- mean(table(data_panel_identifier[,1]), na.rm = TRUE)
  obs_per_group_min <- min(table(data_panel_identifier[,1]), na.rm = TRUE)
  obs_per_group_max <- max(table(data_panel_identifier[,1]), na.rm = TRUE)

  nof_groups <- length(unique(data_panel_identifier[,1]))

  # First demean and then lag the data --------------------------------------

  # Add lags of dependent_vars
  Set_Vars <- cbind(Set_Vars,
                    do.call(
                      rbind,
                      mapply(
                        function(i)
                          panel_lag(Set_Vars[Set_Vars$category == categories[i], c("period", dependent_vars)], lags),
                        1:length(categories),
                        SIMPLIFY = FALSE
                      )
                    ))

  # Add forwardlags of exogenous vars to do the trick suggested by Hahn Kuersteiner: page 1640
  # Eq (2) by adding an additional equation for the exogenous variable, but restriciting its
  # coefficients.
  if (!missing(exog_vars)) {
    Set_Vars <- cbind(Set_Vars,
                      do.call(
                        rbind,
                        mapply(
                          function(i)
                            panel_forward_lag(Set_Vars[Set_Vars$category == categories[i],
                                                       c("period", exog_vars)], lags),
                          1:length(categories),
                          SIMPLIFY = FALSE
                        )
                      ))
  }


  # Add demeaned variables, but also from the lagged dependent variables above
  Set_Vars <- cbind(Set_Vars,
                    do.call(
                      rbind,
                      mapply(
                        function(i)
                          panel_demean(Set_Vars[Set_Vars$category == categories[i],], demean),
                        1:length(categories),
                        SIMPLIFY = FALSE
                      )
                    ))


  # Handle exogenous variables (if necessary) -------------------------------

  # If no exogenous variables then proceed in a standard way.
  if (missing(exog_vars)) {
    dependent.vars <- paste0("demeaned_", dependent_vars)

    # Add all chosen lags to the lagged variables.
    lagged.vars <- c()
    for (l1 in 1:lags){
      lagged.vars <- c(lagged.vars, paste0("demeaned_", "lag", l1,"_", dependent_vars))
    }

  }

  if (!missing(exog_vars)) {
    dependent.vars <- c(paste0("demeaned_", dependent_vars))
    
    lagged.vars <- c()
    for (l1 in 1:lags){
      lagged.vars <- c(lagged.vars, paste0("demeaned_", "lag", l1,"_", dependent_vars))
    }
    
    lagged.vars <- c(lagged.vars,
                     paste0("demeaned_", exog_vars))

  }
  
  Set_Vars_with_NAs <- Set_Vars
  Set_Vars <- na.exclude(Set_Vars)

  #n Anzahl der Individuen (Banken)
  n<-length(unique(Set_Vars$category))

  Bankliste<-unique(Set_Vars$category)

  Matrix.Greek.T<- as.matrix(Set_Vars[,lagged.vars])

  ############################################################################

  # Build OLS on demeaned data:
  M1 <- t(as.matrix(Matrix.Greek.T)) %*% as.matrix(Set_Vars[,lagged.vars])
  M2 <- t(as.matrix(Set_Vars[,lagged.vars]))%*% as.matrix(Set_Vars[,lagged.vars])
  M3 <- t(as.matrix(Set_Vars[, lagged.vars])) %*% as.matrix(Set_Vars[,dependent.vars])
  M4 <- t(as.matrix(Matrix.Greek.T))%*% as.matrix(Set_Vars[,dependent.vars])
  M5 <- t(as.matrix(Matrix.Greek.T))%*% as.matrix(Set_Vars[,lagged.vars])

  # Calculate OLS ------------------------------------------------------------
  sum1.theta.hat <-M2
  sum2.theta.hat<-M3

  # here the equation are in the columns:
  theta.hat.t <-  solve(sum1.theta.hat)%*%as.matrix(sum2.theta.hat)


  # Standard Errors for OLS/restricted OLS estimation ------------------------

  ols.res <- as.matrix(Set_Vars[, dependent.vars]) - as.matrix(Set_Vars[, lagged.vars]) %*% theta.hat.t
  ols.res.squ <- ( t(ols.res) %*% ols.res ) / ( dim(ols.res)[1] -  nrow(theta.hat.t))
  ols_variance <- ols.res.squ %x% solve(sum1.theta.hat)
  ols_variance_sde <-  sqrt(matrix(diag(ols_variance), nrow=length(lagged.vars),byrow=F))
  z_values_ols <- theta.hat.t/ols_variance_sde
  p_values_ols <- 2*pnorm(-abs(z_values_ols))


  

  # Input arguments ---------------------------------------------------------

  inputargs <- list(dependent_vars = dependent_vars,
                    lags = lags,
                    transformation = transformation,
                    Set_Vars = Set_Vars,
                    Set_Vars_with_NAs = Set_Vars_with_NAs,
                    panel_identifier = panel_identifier,
                    nof_observations = nof_observations,
                    obs_per_group_avg = obs_per_group_avg,
                    obs_per_group_min = obs_per_group_min,
                    obs_per_group_max = obs_per_group_max,
                    nof_groups = nof_groups
  )
  

  # Results only OLS------------------------------------------------------------------
  OLS.results <- list(coef =t(theta.hat.t), se = t(ols_variance_sde), pvalues = t(p_values_ols))
  results <- append(list(OLS = OLS.results, residuals = ols.res), c(inputargs))
  class(results) <- "pvarfeols"
  results
  
}

