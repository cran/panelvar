#' Hank Kuehrsteiner Estimator for PVAR Model
#'
#' @param dependent_vars Dependent variables
#' @param lags Number of lags of dependent variables
#' @param exog_vars Exogenous variables
#' @param transformation Demeaining \code{"demean"}
#' @param data Data set
#' @param panel_identifier Vector of panel identifiers
#'
#' @export
#' @description
#' This function estimates a stationary PVAR with fixed effects.
#' @examples
#' data(Dahlberg)
#' ex1_hk <-
#' pvarhk(dependent_vars = c("expenditures", "revenues", "grants"),
#'         lags = 1,
#'         transformation = "demean",
#'         data = Dahlberg,
#'         panel_identifier= c("id", "year"))

pvarhk <- function(dependent_vars,
                           lags = 1,
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
    # Generate Matrix R in order to estimate the restricted OLS-estimator
    l <- length(dependent.vars) - 1
    k <- lags

    R <- matrix(0, nrow = l * k, ncol = (l + k) ^ 2)

    R[1:(l * k), ((l + k) * l + 1):(l ^ 2 + 2 * l * k)] <-
      diag(1, nrow = l * k)
  }

  # If exogenous variables then set up the helper additional equation with parameter restrictions.
  if (!missing(exog_vars)) {
    dependent.vars <- c(paste0("demeaned_", dependent_vars),
                        paste0("demeaned_forwardlag1_", exog_vars))
    lagged.vars <- c(paste0("demeaned_", "lag1_", dependent_vars),
                     paste0("demeaned_", exog_vars))

    # Generate Matrix R in order to estimate the restricted OLS-estimator
    l <- length(dependent.vars) - length(exog_vars)
    k <- nof_exog_vars

    # R is the restriction matrix:
    # ncol of R is the number of parameters in the model
    # nrow is the number of restrictions. e.g. 3 endo + 1 exo. ==>
    # 16 parameters in the extended model, 3 equations for the endos,
    # 1 equation for the exogenous variables.
    if (k == 1){
      R <- matrix(0, nrow = l * k, ncol = (l + k) ^ 2)
      R[1:(l * k), ((l + k) * l + 1):(l ^ 2 + 2 * l * k)] <-
        diag(1, nrow = l * k)
    }

    # Only restrict the lagged endo parameters in the exo-equations to zero.
    if (k > 1){

      R <- matrix(0, nrow = l * k, ncol = (l + k) ^ 2)

      for (i0 in 1:nof_exog_vars){
        R[(i0*l-l+1):(l * i0), ((l + nof_exog_vars) * l + (1+(i0-1)*(l+nof_exog_vars) ) ):
            (((l + nof_exog_vars) * l + (1+(i0-1)*(l+nof_exog_vars) ) )+l-1)] <-
          diag(1, nrow = l)
      }


    }

  }


  Set_Vars <- na.exclude(Set_Vars)


  #n Anzahl der Individuen (Banken)
  n<-length(unique(Set_Vars$category))

  Bankliste<-unique(Set_Vars$category)

  Matrix.Greek.T<- as.matrix(Set_Vars[,lagged.vars])

  List.Te <- matrix(0, nrow = length(Bankliste), ncol = 1)

  # OIDA! ####################################################################
  #i l?uft ?ber die ILZ-Nummern
  # sth for unbalanced panels...to calculate an average T.
  for (i1 in 1:length(Bankliste)){

    # Schreibe Vektor welcher f?r jede ILZ-Nummer die Zeitpunkte, an denen Daten erfasst wurden, ausgibt
    Number.of.periods <- Set_Vars$Date_index[Set_Vars$category==Bankliste[i1]]

    Te <- sum(Set_Vars$category==Bankliste[i1])
    Matrix.Greek.T[Set_Vars$category == Bankliste[i1],]<- 1/Te *
                                  Matrix.Greek.T[Set_Vars$category == Bankliste[i1],]

    List.Te[i1,1] <- Te

  }


  ############################################################################

  # Build OLS Sum1, page 1641 Hahn-Kuersteiner 2002 --------------------------
  M1 <- t(as.matrix(Matrix.Greek.T)) %*% as.matrix(Set_Vars[,lagged.vars])
  M2 <- t(as.matrix(Set_Vars[,lagged.vars]))%*% as.matrix(Set_Vars[,lagged.vars])
  M3 <- t(as.matrix(Set_Vars[, lagged.vars])) %*% as.matrix(Set_Vars[,dependent.vars])
  M4 <- t(as.matrix(Matrix.Greek.T))%*% as.matrix(Set_Vars[,dependent.vars])
  M5 <- t(as.matrix(Matrix.Greek.T))%*% as.matrix(Set_Vars[,lagged.vars])

  Greek.T <- 1/n*M1

  # Calculate OLS ------------------------------------------------------------
  sum1.theta.hat <-M2
  sum2.theta.hat<-M3

  # here the equation are in the columns:
  theta.hat.t <-  solve(sum1.theta.hat)%*%as.matrix(sum2.theta.hat)

  # Calculate restricted OLS (in case of exogenous variables) ----------------
  if (!missing(exog_vars)) {
    X <- diag(1, l + k) %x% as.matrix(Set_Vars[, lagged.vars])
    XX_inv <- diag(1, l + nof_exog_vars) %x% solve(M2)

    vec.theta.hat.R <-
      matrix(theta.hat.t, ncol = 1) + XX_inv %*% t(R) %*% solve(R %*% XX_inv %*% t(R)) %*% (-R %*% matrix(theta.hat.t, ncol = 1))

    theta.hat.R <- t(matrix(vec.theta.hat.R, ncol = l + k))


    colnames(theta.hat.R) <- colnames(t(theta.hat.t))
    rownames(theta.hat.R) <- rownames(t(theta.hat.t))
  }

  if (missing(exog_vars)){
    theta.hat.R <- t(theta.hat.t)
  }

  # Standard Errors for OLS/restricted OLS estimation ------------------------

  ols.res <- as.matrix(Set_Vars[, dependent.vars]) - as.matrix(Set_Vars[, lagged.vars]) %*% t(theta.hat.R)
  ols.res.squ <- ( t(ols.res) %*% ols.res ) / ( dim(ols.res)[1] -  nrow(theta.hat.t))
  ols_variance <- ols.res.squ %x% solve(sum1.theta.hat)
  ols_variance_sde <-  sqrt(matrix(diag(ols_variance), nrow=length(lagged.vars),byrow=F))
  z_values_ols <- t(theta.hat.R)/ols_variance_sde
  p_values_ols <- 2*pnorm(-abs(z_values_ols))

  # Calculate Theta-Double-Hat -----------------------------------------------

  sum2.vec.theta.double.hat <- 1/n * M4
  omega.hat <- Greek.T - theta.hat.R%*%Greek.T%*%t(theta.hat.R)

  # Coefficients of an equation are in the rows of the matrix.
  theta.double.hat <- t( solve(Greek.T) %*% (sum2.vec.theta.double.hat +
                                               1/mean(List.Te[,1])*solve(diag(l+k)-theta.hat.R)%*%omega.hat) )

  # Standard Error per Equation-Block for Hahn Kuehrsteiner ------------------

  hk.res <- as.matrix(Set_Vars[, dependent.vars]) - as.matrix(Set_Vars[, lagged.vars]) %*% t(theta.double.hat)
  hk.res.squ <- ( t(hk.res) %*% hk.res ) / ( dim(hk.res)[1] -  nrow(theta.double.hat))
  hk_variance <- hk.res.squ %x% solve(sum1.theta.hat)
  hk_variance_sde <-  sqrt(matrix(diag(hk_variance), nrow=length(lagged.vars),byrow=F))
  z_values_hk <- theta.double.hat/hk_variance_sde
  p_values_hk <- 2*pnorm(-abs(z_values_hk))

  # Results ------------------------------------------------------------------
  OLS.results <- list(coef = t(theta.hat.R), se = ols_variance_sde, pvalues = p_values_ols)
  HK.results  <- list(coef = theta.double.hat, se = hk_variance_sde, pvalues = p_values_hk)
  results <- list(OLS = OLS.results, HK = HK.results)
  results
}



