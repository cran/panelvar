#' GMM Estimation of Panel VAR Models
#'
#' @param dependent_vars Dependent variables
#' @param lags Number of lags of dependent variables
#' @param predet_vars Predetermined variables
#' @param exog_vars Exogenous variables
#' @param transformation First-difference \code{"fd"} or forward orthogonal deviations \code{"fod"}
#' @param data Data set
#' @param panel_identifier Vector of panel identifiers
#' @param steps \code{"onestep"}, \code{"twostep"} or \code{"mstep"}   estimation
#' @param system_instruments System GMM estimator
#' @param system_constant Constant only available with the System GMM estimator in each equation
#' @param pca_instruments Apply PCA to instruments matrix
#' @param pca_eigenvalue Cut-off eigenvalue for PCA analysis
#' @param max_instr_dependent_vars Maximum number of instruments for dependent variables
#' @param max_instr_predet_vars Maximum number of instruments for predetermined variables
#' @param min_instr_dependent_vars Minimum number of instruments for dependent variables
#' @param min_instr_predet_vars Minimum number of instruments for predetermined variables
#' @param collapse Use collapse option
#' @param tol relative tolerance to detect zero singular values in \code{"ginv"}
#' @param progressbar show progress bar
#' @import Matrix
#' @import ggplot2
#' @import texreg
#' @import knitr
#' @importFrom MASS ginv
#' @importFrom stats coef quantile vcov pchisq pnorm na.exclude prcomp
#' @importFrom progress progress_bar
#' @export
#' @description 
#' Estimates a panel vector autoregressive (PVAR) model with fixed effects.
#' @return A \code{pvargmm} object containing the estimation results.
#' @details
#' The first vector autoregressive panel model (PVAR) was introduced by Holtz-Eakin et al. (1988). Binder et al. (2005) extend their equation-by-equation estimator for a PVAR model with only endogenous variables that are lagged by one period. We further improve this model in Sigmund and Ferstl (2021) to allow for \eqn{p} lags of \eqn{m} endogenous variables, \eqn{k} predetermined variables and \eqn{n} strictly exogenous variables.
#' 
#' Therefore, we consider the following stationary PVAR with fixed effects.
#' 
#' 	\ifelse{html}{\out{<b>y</b><sub><i>i,t</i></sub> = <b>&mu;</b><sub><i>i</i></sub> + &sum;<sub><i>l=1</i></sub><sup>p</sup><b>A</b><sub><i>l</i></sub><b>y</b><sub><i>i,t-l</i></sub> + <b>B</b><b>x</b><sub><i>i,t</i></sub> + <b>C</b><b>s</b><sub><i>i,t</i></sub> + <b>&epsilon;</b><sub><i>i,t</i></sub><p></p>}}{\deqn{\mathbf{y}_{i,t} = \mathbf{\mu}_i + \sum_{l = 1}^{p}\mathbf{A}_l\mathbf{y}_{i,t-l} + \mathbf{B}\mathbf{x}_{i, t} + \mathbf{C}\mathbf{s}_{i,t} + \mathbf{\epsilon}_{i,t}}}
#'
#'
#' \ifelse{html}{\out{Let <b>y</b><sub><i>i,t</i></sub> &isin; &real;<sup><i>m</i></sup> be an <i>m&times;1</i> vector of <b>endogenous variables</b> for the <i>i</i>th cross-sectional unit at time <i>t</i>. Let <b>y</b><sub><i>i,t-l</i></sub> &isin; &real;<sup><i>m</i></sup> be an <i>m&times;1</i> vector of <b>lagged endogenous variables</b>. Let <b>x</b><sub><i>i,t</i></sub> &isin; &real;<sup><i>k</i></sup> be an <i>k&times;1</i> vector of <b>predetermined variables</b> that are potentially correlated with past errors. Let <b>s</b><sub><i>i,t</i></sub> &isin; &real;<sup><i>n</i></sup> be an <i>n&times;1</i> vector of <b>strictly exogenous variables</b> that neither depend on <b>&epsilon;</b><sub><i>i,t</i></sub> nor on <b>&epsilon;</b><sub><i>i,t-s</i></sub> for <i>s = 1,&hellip;,T</i>. The idiosyncratic error vector <b>&epsilon;</b><sub><i>i,t</i></sub> &isin; &real;<sup><i>m</i></sup> is assumed to be well-behaved and independent from both the regressors <b>x</b><sub><i>i,t</i></sub>  and <b>s</b><sub><i>i,t</i></sub> and the individual error component <b>&mu;</b><sub><i>i</i></sub>. Stationarity requires that all unit roots of the PVAR model fall inside the unit circle, which therefore places some constraints on the <b>fixed effect</b> <b>&mu;</b><sub><i>i</i></sub>. The cross section <i>i</i> and the time section <i>t</i>  are defined as follows: <i>i = 1,&hellip;,N</i> and <i>t = 1,&hellip;T</i>. In this specification we assume parameter homogeneity for <b>A</b><sub><i>l</i></sub> <i>(m&times;m)</i>, <b>B</b> <i>(m&times;k)</i> and <b>C</b> <i>(m&times;n)</i> for all <i>i</i>.}}{\eqn{\mathbf{I}_m} denotes an \eqn{m\times m} identity matrix. Let \eqn{\mathbf{y}_{i,t} \in \R^m} be an \eqn{m\times 1} vector of endogenous variables for the \eqn{i}th cross-sectional unit at time \eqn{t}. Let \eqn{\mathbf{y}_{i,t-l} \in \R^m} be an \eqn{m\times 1} vector of lagged endogenous variables. Let \eqn{\mathbf{x}_{i,t} \in \R^k} be an \eqn{k \times 1} vector of predetermined variables that are potentially correlated with past errors. Let \eqn{\mathbf{s}_{i,t} \in \R^n} be an \eqn{n \times 1} vector of strictly exogenous variables that neither depend on \eqn{\epsilon_t} nor on \eqn{\epsilon_{t-s}} for \eqn{s = 1,\dots,T}. The idiosyncratic error vector \eqn{\mathbf{\epsilon}_{i,t} \in \R^m} is assumed to be well-behaved and independent from both the regressors \eqn{\mathbf{x}_{i,t}} and \eqn{\mathbf{s}_{i,t}} and the individual error component \eqn{\mathbf{\mu}_i}. Stationarity requires that all unit roots of the PVAR model fall inside the unit circle, which therefore places some constraints on the fixed effect \eqn{\mathbf{\mu}_i}. The cross section \eqn{i} and the time section \eqn{t} are defined as follows: \eqn{i = 1,2,...,N} and \eqn{t = 1,2,...,T}. In this specification we assume parameter homogeneity for \eqn{\mathbf{A}_l (m\times m)}, \eqn{\mathbf{B} (m \times k)} and \eqn{\mathbf{C} (m\times n)} for all \eqn{i}.}
#' 
#' A PVAR model is hence a combination of a single equation dynamic panel model (DPM) and a vector autoregressive model (VAR).
#' 
#' 
#' First difference and system GMM estimators for single equation dynamic panel data models have been implemented in the STATA package \code{xtabond2} by Roodman (2009) and some of the features are also available in the R package \pkg{plm}.
#' 
#' For more technical details on the estimation, please refer to our paper Sigmund and Ferstl (2021).
#' 
#' There we define the first difference moment conditions (see Holtz-Eakin et al., 1988; Arellano and Bond, 1991), formalize the ideas to reduce the number of moment conditions by linear transformations of the instrument matrix and define the one- and two-step GMM estimator. Furthermore, we setup the system moment conditions as defined in Blundell and Bond (1998) and present the extended GMM estimator. In addition to the GMM-estimators we contribute to the literature by providing specification tests (Hansen overidentification test, lag selection criterion and stability test of the PVAR polynomial) and classical structural analysis for PVAR models such as orthogonal and generalized impulse response functions, bootstrapped confidence intervals for impulse response analysis and forecast error variance decompositions. Finally, we implement the first difference and the forward orthogonal transformation to remove the fixed effects.
#' 
#' @references
#' Arellano, M., Bond, S. (1991) Some Tests of Specification for Panel Sata: Monte Carlo Evidence and an Application to Employment Equations \emph{The Review of Economic Studies}, \bold{58}(2), 277--297, \doi{10.2307/2297968}
#'
#' Binder M., Hsiao C., Pesaran M.H. (2005) Estimation and Inference in Short Panel Vector Autoregressions with Unit Roots and Cointegration \emph{Econometric Theory}, \bold{21}(4), 795--837, \doi{10.1017/S0266466605050413}
#'
#' Blundell R., Bond S. (1998). Initial Conditions and Moment Restrictions in Dynamic Panel Data Models \emph{Journal of Econometrics}, \bold{87}(1), 115--143, \doi{10.1016/S0304-4076(98)00009-8}
#'
#' Holtz-Eakin D., Newey W., Rosen H.S. (1988) Estimating Vector Autoregressions with Panel Data, \emph{Econometrica}, \bold{56}(6), 1371--1395, \doi{10.2307/1913103}
#'
#' Roodman, D. (2009) How to Do xtabond2: An Introduction to Difference and System GMM in Stata \emph{The Stata Journal}, \bold{9}(1), 86--136, \url{https://www.stata-journal.com/article.html?article=st0159}
#' 
#' Sigmund, M., Ferstl, R. (2021) Panel Vector Autoregression in R with the Package panelvar \emph{The Quarterly Review of Economics and Finance} \doi{10.1016/j.qref.2019.01.001}
#' 
#' @seealso
#' \code{\link{stability}} for stability tests
#' 
#' \code{\link{oirf}} and \code{\link{girf}} for orthogonal and generalized impulse response functions (including bootstrapped confidence intervals)
#' 
#' \code{\link{coef.pvargmm}}, \code{\link{se}}, \code{\link{pvalue}}, \code{\link{fixedeffects}} for extrator functions for the most important results
#' 
#' \code{\link{fevd_orthogonal}} for forecast error variance decomposition
#' 
#' @examples
#' \dontrun{
#' library(panelvar)
#' data(abdata)
#' ex3_abdata <-pvargmm(
#'  dependent_vars = c("emp"),
#'  lags = 4,
#'  predet_vars = c("wage"),
#'  exog_vars = c("cap"),
#'  transformation = "fd",
#'  data = abdata,
#'  panel_identifier = c("id", "year"),
#'  steps = c("twostep"),
#'  system_instruments = TRUE,
#'  max_instr_dependent_vars = 99,
#'  max_instr_predet_vars = 99,
#'  min_instr_dependent_vars = 2L,
#'  min_instr_predet_vars = 1L,
#'  collapse = FALSE
#' )
#' }
#' data("ex3_abdata")
#' summary(ex3_abdata)
#' 
#' data("Dahlberg")
#' \dontrun{
#' ex1_dahlberg_data <- pvargmm(dependent_vars = c("expenditures", "revenues", "grants"),
#'                              lags = 1,
#'                              transformation = "fod",
#'                              data = Dahlberg,
#'                              panel_identifier=c("id", "year"),
#'                              steps = c("twostep"),
#'                              system_instruments = FALSE,
#'                              max_instr_dependent_vars = 99,
#'                              max_instr_predet_vars = 99,
#'                              min_instr_dependent_vars = 2L,
#'                              min_instr_predet_vars = 1L,
#'                              collapse = FALSE
#' )
#' }
#' data("ex1_dahlberg_data")
#' summary(ex1_dahlberg_data)

pvargmm <-
  function(dependent_vars,
           lags,
           predet_vars,
           exog_vars,
           transformation = "fd",
           data,
           panel_identifier = c(1, 2),
           steps,
           system_instruments = FALSE,
           system_constant = TRUE,
           pca_instruments = FALSE,
           pca_eigenvalue = 1,
           max_instr_dependent_vars,
           max_instr_predet_vars,
           min_instr_dependent_vars = 2L,
           min_instr_predet_vars = 1L,
           collapse = FALSE,
           tol = 1e-09,
           progressbar = TRUE) {


    # Drop factor levels from data --------------------------------------------
    data <- droplevels(data)

    # Select only required variables from data --------------------------------
    required_vars <- c(dependent_vars)
    if(!missing(predet_vars)) required_vars <- c(required_vars, predet_vars)
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
    if(missing(predet_vars)) {
      nof_predet_vars <- 0
    } else {
      nof_predet_vars <- length(predet_vars)
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

    # Use maximum possible instruments ----------------------------------------
    if (missing(max_instr_dependent_vars))
      max_instr_dependent_vars <- nof_periods
    if (missing(max_instr_predet_vars))
      max_instr_predet_vars <- nof_periods

    # Balance the panel filling missing values with NA ------------------------
    Set_Vars <-
      merge(
        expand.grid(period = periods, category = categories)[, c(2, 1)],
        Set_Vars,
        by = c("category", "period"),
        all.x = TRUE
      )

    Set_Vars <- Set_Vars[order(Set_Vars$category, Set_Vars$period,decreasing = FALSE),]
    Set_Vars$category <- factor(Set_Vars$category)
    Set_Vars$period <- factor(Set_Vars$period)
    Set_Vars_with_NAs <- Set_Vars
    
    # Add lags of dependent variables -----------------------------------------
    # - xtabond sets the lagged values to NA if the non-lagged values is NA
    #   in the last period of a specific variable.
    # - nof_dependent_vars >= 1

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

    # FD- or FOD transformation -----------------------------------------------

    if (system_instruments & transformation == "fod") {
      fd_Set_Vars <- do.call(
        rbind,
        mapply(
          function(i)
            panel_fd_fod(Set_Vars[Set_Vars$category == categories[i],],"fd"),
          1:length(categories),
          SIMPLIFY = FALSE
        )
      )
    }

    Set_Vars <- cbind(Set_Vars,
                      do.call(
                        rbind,
                        mapply(
                          function(i)
                            panel_fd_fod(Set_Vars[Set_Vars$category == categories[i],],transformation),
                          1:length(categories),
                          SIMPLIFY = FALSE
                        )
                      ))
    
    # Build stats here: 
    # Descriptive statistics --------------------------------------------------
    nof_observations <- dim(na.exclude(Set_Vars))[1]
    data_panel_identifier <- data[panel_identifier]
    
    obs_per_group_avg <- mean(table(factor(na.exclude(Set_Vars)$category)), na.rm = TRUE)
    obs_per_group_min <- min(table(factor(na.exclude(Set_Vars)$category)), na.rm = TRUE)
    obs_per_group_max <- max(table(factor(na.exclude(Set_Vars)$category)), na.rm = TRUE)
    
    nof_groups <- length(unique(na.exclude(Set_Vars)$category))
    
    
    # -------------------------------------------------------------------------

    # only for System fod instruments
    if (system_instruments & transformation == "fod") {
      Set_Vars <- cbind(Set_Vars, fd_Set_Vars)
      fd_dependent_vars <- paste("fd",c(dependent_vars), sep="_")
      if (missing(predet_vars)) {
        fd_predet_vars <- NULL
      } else {
        fd_predet_vars <- paste("fd",c(predet_vars), sep="_")
      }
    }

    all_variable_names <- names(Set_Vars)[-c(1, 2)] # NOTE: this variable is redundant
    transformed_dependent_vars <- paste(transformation,c(dependent_vars), sep="_")

    # NOTE: not pretty
    transformed_lagged_dependent_vars <- c()
    for (i0 in 1:lags){
      transformed_lagged_dependent_vars <- c(transformed_lagged_dependent_vars,
                                             paste(transformation, "_lag", i0, "_", c(dependent_vars), sep=""))
    }

    if (missing(predet_vars)) {
      transformed_predet_vars <- NULL
    } else {
      transformed_predet_vars <- paste(transformation,c(predet_vars), sep="_")
    }

    if (missing(exog_vars)) {
      transformed_exog_vars <- NULL
    } else {
      transformed_exog_vars <- paste(transformation,c(exog_vars), sep="_")
    }

    names_lags <- paste("lag", 1:lags, sep="")

    names_lagged_vars <- NULL
    for(i in 1:lags) names_lagged_vars <- c( names_lagged_vars, paste(names_lags[i], dependent_vars, sep="_") )

    transformed_lagged_vars <- paste(transformation,names_lagged_vars, sep="_")

    if (system_instruments & transformation == "fd") {
      fd_dependent_vars <- transformed_dependent_vars
      if (missing(predet_vars)) {
        fd_predet_vars <- NULL
      } else {
        fd_predet_vars <- transformed_predet_vars
      }
    }

    
    
    # Onestep estimation (initialization) ------------------------------------

    # initalize Q_i matrix (see Binder eq 6.2)
    # explanation in paper required, more details on dimension

    # Set the dimension of Q unrestricted.
    # Start with the rows: M = M_lagged_endo + M_predetermined + M_exogenous
    M_lagged_endo <- rows_Q_lagged_endogenous_vars(nof_periods,
                                                   nof_dependent_vars,
                                                   max_instr_dependent_vars = nof_periods,
                                                   lags,
                                                   2)

    M_predetermined <- rows_Q_predetermined_vars(nof_periods,
                                                 nof_predet_vars,
                                                 max_instr_predet_vars = nof_periods,
                                                 lags,
                                                 1)

    M_pred <- M_lagged_endo * (!missing(predet_vars))

    M_exogenous <- nof_exog_vars * (!missing(exog_vars))

    M <- M_lagged_endo + M_predetermined + M_exogenous

    #prepare matrix Lambda: Middle Matrix in  D=...
    delta_W<- replicate(nof_categories, Matrix(0,nrow = nof_periods - lags - 1, ncol = length(transformed_dependent_vars)))
    delta_W_NA <- replicate(nof_categories, Matrix(0,nrow = nof_periods - lags - 1, ncol = length(transformed_dependent_vars)))

    delta_W_minus <- vector("list", nof_categories) # RF: TODO improve memory allocation
    delta_W_minus_NA <- vector("list", nof_categories) # RF: TODO improve memory allocation

    Q <- replicate(nof_categories, Matrix(0, nrow = M, ncol = nof_periods - 1 - lags))

    Q_endo <- replicate(nof_categories, Matrix(0, nrow = M_lagged_endo, ncol = nof_periods - 1 - lags))

    Q_predetermined <- replicate(nof_categories, Matrix(0, nrow = M_predetermined, ncol = nof_periods - 1 - lags))

    Q_exogenous <- replicate(nof_categories, Matrix(0, nrow = M_exogenous, ncol = nof_periods - 1 - lags))

    Q_new <- replicate(nof_categories, Matrix(0, nrow = M_lagged_endo + M_predetermined +
                                                        M_exogenous, ncol = nof_periods - 1 - lags))

    # Set up restriction matrices that are applied afterwards.
    # This is faster than applying all the user specified restrictions on Q in the process of setting up Q.
    Endo_Transition_Lmin_Lmax <-  Endo_Lmin_Lmax_transition_matrix(lags = lags,
                                                                   min_instr_dependent_vars = min_instr_dependent_vars,
                                                                   max_instr_dependent_vars = max_instr_dependent_vars,
                                                                   nofperiods = nof_periods,
                                                                   nof_dependent_vars = nof_dependent_vars,
                                                                   collapse = collapse)


    Predetermined_Transition_Lmin_Lmax <-
      Predetermined_Lmin_Lmax_transition_matrix(
        lags = lags,
        min_instr_predet_vars = min_instr_predet_vars,
        max_instr_predet_vars = max_instr_predet_vars,
        nofperiods = nof_periods,
        nof_predet_vars = nof_predet_vars,
        nof_dependent_vars = nof_dependent_vars, # not necessary
        collapse = collapse)

    position_exogenous <- ncol(Endo_Transition_Lmin_Lmax) + 1
    if(!missing(predet_vars)){position_exogenous <- position_exogenous +ncol(Predetermined_Transition_Lmin_Lmax)}
    
    # Begin onestep computation loop -----------------------------------------
    # Output: Q, delta.W, delta_W_minus

    if (progressbar) {
      pb <-
        progress::progress_bar$new(format = "One-step estimation: [:bar] Iteration :current of :total" ,
                                   clear = TRUE,
                                   total = nof_categories)
    }

    for (i in 1:nof_categories) {
      if (progressbar) pb$tick()
      #print(paste0("First step iteration ",i ))

      # build delta_W and delta_W_minus
      Set_Vars_i <- Set_Vars[Set_Vars$category==categories[i],]
      #Set_Vars_i[is.na(Set_Vars_i)]<-0

      delta_W[[i]] <- Matrix(as.matrix(Set_Vars_i[-(1:(lags+1)), transformed_dependent_vars]))

      # Add parts for system instruments.
      if(system_instruments == TRUE & transformation == "fod"){

        delta_W[[i]] <- rbind(0, delta_W[[i]])

      }

      #delta_W[[i]] <- Matrix(Set_Vars_i[, transformed_dependent_vars])[-(1:(lags+1)),]
      # endogenous and exogenous components of delta_W_minus
      delta_W_minus_dep <- as.matrix(Set_Vars_i[, transformed_lagged_dependent_vars])

      #build lags of delta.W
      delta_W_minus_dep_lag<- as.matrix(delta_W_minus_dep[(lags+2):(nrow(delta_W_minus_dep)),])

      delta_W_minus_dep <- delta_W_minus_dep_lag

      # stick exogenous and endogenous components of delta_W_minus together to one matrix
      if(!missing(predet_vars) | !missing(exog_vars)){
        delta_W_minus_exog <- as.matrix(Set_Vars_i[, c(transformed_predet_vars, transformed_exog_vars)])[-(1:(lags+1)),]
        delta_W_minus[[i]] <- Matrix(cbind(delta_W_minus_dep, delta_W_minus_exog))
      }else{delta_W_minus[[i]] <- Matrix(delta_W_minus_dep)}

      # NOTE: fix nulls
      colnames(delta_W_minus[[i]]) <- c(transformed_lagged_vars, transformed_predet_vars, transformed_exog_vars)

      # Add parts for system instruments.
      if(system_instruments == TRUE &  transformation == "fod"){

        delta_W_minus[[i]] <- rbind(0, delta_W_minus[[i]])

      }

      # build matrix Q'

      q <- NULL
      q_pred <- NULL

      l<-0
      l_pred<-0

      #######################################################################################
      # Start Set up W and W.minus.dep for FD System GMM
      #######################################################################################
      if (system_instruments == TRUE){

        W <- as.matrix(Set_Vars_i[Set_Vars_i[,1]==categories[i], names(Set_Vars_i)%in% dependent_vars])

        W <- as.matrix(W[(lags+1):nrow(W),])


        # Build endogenous and exogenous components of the lower part of W.minus
        predet_and_exog_vars <- NULL
        if(!missing(predet_vars)) predet_and_exog_vars <- c(predet_and_exog_vars, predet_vars)
        if(!missing(exog_vars)) predet_and_exog_vars <- c(predet_and_exog_vars, exog_vars)

        if(!missing(predet_vars) | !missing(exog_vars)) {
          W.minus.exog <-
            as.matrix(Set_Vars_i[Set_Vars_i[,1]==categories[i], names(Set_Vars_i) %in% predet_and_exog_vars])
        }


        #build lags of lower part of W.minus
        W.minus.dep <- as.matrix(Set_Vars_i[Set_Vars_i[,1]==categories[i], names(Set_Vars_i) %in% dependent_vars])
        W.minus.dep.lag <- as.matrix(W.minus.dep[lags:(nrow(W.minus.dep)-1),])

        if(lags>=2){
          for(l in 1:(lags-1)){
            W.minus.dep.lag <- cbind(W.minus.dep.lag, as.matrix(W.minus.dep[(lags-l):(nrow(W.minus.dep)-1-l),]) )
          }
        }

        W.minus.dep <- W.minus.dep.lag

        if(!missing(predet_vars) | !missing(exog_vars)) W.minus.exog <- as.matrix(W.minus.exog[(lags+1):nrow(W.minus.exog),])

        #stick exogenous and endogenous components of the lower part of W.minus together to one matrix
        if(!missing(predet_vars) | !missing(exog_vars)){
          W.minus <- cbind(W.minus.dep, W.minus.exog)
        }else{W.minus <- W.minus.dep}


        # Improve next two lines!!!!
        if (system_constant == TRUE){
          cons <- matrix(1,nrow=nrow(W.minus), ncol=1)

          cons <- rbind(matrix(0, nrow = nrow(delta_W_minus[[i]])+nrow(W.minus)-nrow(cons), ncol=1), cons)
        }

        # Be careful here when you combine the two parts:
        delta_W[[i]] <- rbind(delta_W[[i]], W)
        delta_W_minus[[i]] <- rbind(delta_W_minus[[i]], W.minus)

        if (system_constant == TRUE){
          delta_W_minus[[i]] <- cbind(delta_W_minus[[i]], cons)
        }

      }
      #######################################################################################
      # End Set up W and W.minus.dep for FD System GMM
      #######################################################################################
      q <- NULL
      l<-0
      # new z-loop for endogenous variables:
      for (z in 1:(nof_periods-1-lags)){
        #print(z)
        l<-l+length(q)
        #l_pred <- l_pred+length(q_pred)
        #print(l_pred)

        #add endogenous variables: ADDED by MICHAEL
        if (z==1){
          for(j in 1:lags){
            q<-rbind(t(as.matrix(Set_Vars_i[Set_Vars_i[,2]==periods[j],dependent_vars])),q)
          }
        }else{
          q <- rbind(t(as.matrix(Set_Vars_i[Set_Vars_i[,2]==periods[lags+z-1],dependent_vars])),q)
        }

        #add endogenous variables to the matrix Q'
        Q_endo[[i]][(l+1):(l+length(q)),z]<-q

      }
      Q_endo[[i]][is.na(Q_endo[[i]])]<-0

      Q_endo[[i]] <- t( t(Q_endo[[i]]) %*% Endo_Transition_Lmin_Lmax )

      q_pred <- NULL
      l_pred<-0
      #l <- 0
      #length(q) <- 0
      # new z loop for predetermined variables: Improve
      if (!missing(predet_vars)){
        for (z in 1:(nof_periods-1-lags)) {

          l_pred <- l_pred+length(q_pred)

          if (z>1){
            q_pred <- rbind(t(as.matrix(Set_Vars_i[Set_Vars_i[,2]==
                                                     periods[z+lags], predet_vars])), q_pred)
          }

          if (z==1) {
            for (j in 1:(lags+1)){
              q_pred <- rbind(t(as.matrix(Set_Vars_i[Set_Vars_i[,2]==periods[j], predet_vars])),q_pred)
            }
          }

          Q_predetermined[[i]][(1 + l_pred):(l_pred + length(q_pred)),z] <- q_pred
        }
        Q_predetermined[[i]][is.na(Q_predetermined[[i]])]<-0

        Q_predetermined[[i]] <- t( t(Q_predetermined[[i]]) %*% Predetermined_Transition_Lmin_Lmax )

      }

      # Set rows in Q for exogenous variables.
      if (!missing(exog_vars)){
        line_ex <- t(as.matrix(Set_Vars[Set_Vars[,1]==categories[i],transformed_exog_vars]))
        line_ex <- line_ex[,(lags+2):ncol(line_ex)]
        Q_exogenous[[i]][(1:nof_exog_vars), ] <- line_ex
        Q_exogenous[[i]][is.na(Q_exogenous[[i]])]<-0
      }

      # Put the full Q_new matrix together: Improve
      ########################################################
      Q_new[[i]] <- Q_endo[[i]]

      if (!missing(predet_vars)){
        Q_new[[i]] <- rbind(Q_endo[[i]], Q_predetermined[[i]])
      }

      if (!missing(exog_vars)){
        Q_new[[i]] <- rbind(Q_new[[i]], Q_exogenous[[i]])
      }

      Q[[i]] <- Q_new[[i]]
      #######################################################

      if(system_instruments == TRUE & transformation == "fod"){

        Q[[i]] <- cbind(0, Q[[i]])

      }
      ############################################################################

      #######################################################################################
      # End Standard first difference gmm instruments:
      #######################################################################################

      #######################################################################################
      # Start System GMM instruments
      #######################################################################################
      if (system_instruments == TRUE){


        if(collapse==FALSE){

          P2a<-matrix(0, nrow = nof_dependent_vars*(nof_periods-2), ncol = nof_periods-1)

          for (z in 1:(nof_periods-2)){
            P2a[(1+(z-1)*nof_dependent_vars):(z*nof_dependent_vars),1+z] <-
              t(as.matrix(Set_Vars_i[Set_Vars_i[,1]==categories[i] & Set_Vars_i[,2]==periods[z+1], fd_dependent_vars]))
          }

          if(lags==2) P2a <- P2a[,-1]

          if(lags>2) P2a <- P2a[-(1:(nof_dependent_vars*(lags-2))),-(1:(lags-1))]


          # Create matrix with predetermined variables as instruments
          if(!missing(predet_vars)){
            P2b<-matrix(0,nrow = nof_predet_vars*(nof_periods-1),ncol=nof_periods-1)

            for (z in 1:(nof_periods-1)){
              P2b[(1+(z-1)*nof_predet_vars):(z*nof_predet_vars),z] <-
                t(as.matrix(Set_Vars_i[Set_Vars_i[,1]==categories[i] & Set_Vars_i[,2]==periods[z+1], fd_predet_vars]))
            }
            if(lags>1) P2b <- P2b[-(1:(nof_predet_vars*(lags-1))),-(1:(lags-1))]
          }



        } else{

          P2a<-matrix(0, nrow = nof_dependent_vars, ncol = nof_periods-1)


          P2a[,2:(nof_periods-1)] <-
            t(as.matrix(Set_Vars_i[Set_Vars_i[,1]==categories[i] & Set_Vars_i[,2] %in% periods[2:(nof_periods-1)], fd_dependent_vars]))


          if(lags==2) P2a <- P2a[,-1]


          if(lags>2) P2a <- P2a[,-(1:(lags-1))]

          #Create matrix with predetermined variables as instruments
          if(!missing(predet_vars)){

            P2b<-matrix(0,nrow= nof_predet_vars , ncol=nof_periods-1)

            P2b <- t(as.matrix(Set_Vars_i[Set_Vars_i[,1]==categories[i] & Set_Vars_i[,2] %in% periods[2:nof_periods], fd_predet_vars]))


            if(lags>1) P2b <- P2b[,-(1:(lags-1))]
          }
        }

        #Create matrix with 1 as intstruments - for the constant
        if (system_constant == TRUE){
          P2c <- matrix(1, ncol=nof_periods - lags, nrow=1)
        }

        if(!missing(predet_vars) & system_constant == TRUE){
          P2 <- rbind(P2a, P2b, P2c)
        }

        if(!missing(predet_vars) & system_constant == FALSE){
          P2 <- rbind(P2a, P2b)
        }

        if (missing(predet_vars) & system_constant == TRUE){P2 <- rbind(P2a, P2c)}

        if (missing(predet_vars) & system_constant == FALSE){P2 <- P2a}

        #######################################################################################
        # End System GMM instruments
        #######################################################################################

        # Combine fd_instruments with system_instruments

        P1 <- Q[[i]]
        P<-cbind(P1,matrix(0,ncol=ncol(P2),nrow=nrow(P1)))
        P<-rbind(P, cbind(matrix(0,ncol=ncol(P1),nrow=nrow(P2)),P2))

        # add exogenous variables as instruments for the levels
        if(!missing(exog_vars)){
          P[(nrow(P1)+1-nof_exog_vars):nrow(P1),(ncol(P1)+1):(ncol(P1)+ncol(P2))] <-
            t(as.matrix(Set_Vars_i[Set_Vars_i[,1]==categories[i],exog_vars], ncol=nof_exog_vars))[,-(1:lags)]
        }

        Q[[i]] <- P

      }
      #######################################################################################
      # End Combine System GMM instruments
      #######################################################################################



      # NOTE: not pretty
      # uses zero as instruments in the matrix Q' for non available observations
      zero_col <- as.vector(which(as.logical(rowSums(is.na(cbind(as.matrix(delta_W[[i]]), as.matrix(delta_W_minus[[i]]) ))))!=0))
      if (length(zero_col)!=0){
        Q[[i]][,zero_col] <- 0
      }

      #replaces NAs by 0 in Q, delta_W and delta_W_minus
      # reference in paper
      Q[[i]][is.na(Q[[i]])]<-0

      if (sum(rowSums(is.na(delta_W[[i]])))) {
        delta_W[[i]][as.logical(!!rowSums(is.na(delta_W[[i]]))),] <- NA
      }
      delta_W_NA[[i]] <- delta_W[[i]]
      delta_W[[i]][is.na(delta_W[[i]])] <- 0

      if (sum(rowSums(is.na(delta_W_minus[[i]])))) {
        delta_W_minus[[i]][as.logical(rowSums(is.na(delta_W_minus[[i]]))),] <- NA

      }
      delta_W_minus_NA[[i]] <- delta_W_minus[[i]]
      delta_W_minus[[i]][is.na(delta_W_minus[[i]])] <- 0


    }

    #######################################################################################
    # PCA on Instruments ##################################################################
    if(pca_instruments == TRUE){

      Q <- pca_transformation(corrected_instruments_F = Q,
                              nof_categories = nof_categories,
                              nof_predet_vars = nof_predet_vars,
                              nof_exog_vars = nof_exog_vars,
                              system_instruments = system_instruments,
                              system_constant = system_constant,
                              position_exo = position_exogenous,
                              pca_eigenvalue = pca_eigenvalue)

    }
    #######################################################################################

    ############################################################################
    # End onestep computation loop -------------------------------------------#
    ############################################################################

    ############################################################################
    # Construct weighting matrix V/D
    # (TODO: make this separate function)
    ############################################################################

    if (system_instruments == TRUE){

      if (transformation == "fd"){

        M1 <- diag(-1, ncol(P1)+1)

        for(ind in 1:(ncol(P1))){
          M1[ind, ind+1] <- 1
        }

        M1 <- M1[-nrow(M1),]

        M2 <- diag(1, nrow = ncol(P1)+1, ncol = ncol(P1) + 1)

        V <- rbind(M1, M2)

      }


      if (transformation == "fod"){

        M1 <- matrix(0,nrow = nof_periods, ncol=nof_periods)

        for(d in 1:(nrow(M1))){

          diag(M1)[d] <- sqrt(nof_periods-d)/sqrt(nof_periods+1-d)

          if (1+d <= nrow(M1)){
            M1[d,(1+d):ncol(M1)] <- - 1/sqrt((nof_periods+1-d)*(nof_periods-d))
          }

        }

        M1 <- M1[ (1+lags-1) :(nof_periods-1) , 1:nof_periods]

        M2 <- diag(1, nof_periods)

        M2[,1:(lags)] <- 0

        M2 <- M2[ c((1+lags):(nof_periods)) , ]

        V <- rbind(M1, M2)

       }

    }

    # Weighting for only FD-GMM: Adjusted.
    if (system_instruments == FALSE){
      if (transformation == "fd"){

        M1 <- diag(-1, ncol(Q[[1]])+1)

        for(ind in 1:(ncol(Q[[1]]))){
          M1[ind, ind+1] <- 1
        }

        M1 <- M1[-nrow(M1),]

        V <- M1


      }

      if (transformation == "fod"){

        M1 <- matrix(0,nrow = nof_periods, ncol=nof_periods)

        for(d in 1:(nrow(M1))){

          diag(M1)[d] <- sqrt(nof_periods-d)/sqrt(nof_periods+1-d)

          if (1+d <= nrow(M1)){
            M1[d,(1+d):ncol(M1)] <- - 1/sqrt((nof_periods+1-d)*(nof_periods-d))
          }

        }

        M1 <- M1[ (1+lags-1) :(nof_periods-1) , 1:nof_periods]
        M1 <- M1[1:(nof_periods-1-lags), ]

        V <- M1


      }
    }
    ############################################################################
    if (progressbar) {
      pb <- progress::progress_bar$new(format = " One-step estimation: Matrix inversion", total = 2, clear = TRUE)
      pb$tick()
    }

    S_ZX_a <- mapply(function(i) Matrix(Q[[i]] %*% delta_W_minus[[i]], sparse = TRUE), 1:nof_categories)
    sum_S_ZX_a <-Reduce('+', S_ZX_a)

    sum_S_Zy_a <-Reduce('+', mapply(function(i) Matrix( Q[[i]]%*%delta_W[[i]], sparse = TRUE), 1:nof_categories, SIMPLIFY = FALSE))

    sum_Lambda <- Reduce('+', mapply(function(i) Matrix(   Q[[i]] %*%V %*% t(V) %*% t(Q[[i]]), sparse = TRUE), 1:nof_categories, SIMPLIFY = FALSE))

    #sum_Lambda_vec <- sum_Lambda %x% diag(1,nof_dependent_vars)
    sum_Lambda_vec <- sum_Lambda %x% Diagonal(nof_dependent_vars)

    if (min(abs(eigen(sum_Lambda)$values)) < tol){
      #Lambda_inv <-  MASS::ginv(as.matrix(sum_Lambda))
      Lambda_inv_vec <- MASS::ginv(as.matrix(sum_Lambda_vec))
      warning("The matrix Lambda is singular, therefore the general inverse is used")
    } else{
      #Lambda_inv <-  solve(sum_Lambda)
      Lambda_inv_vec <-  solve(as.matrix(sum_Lambda_vec))
    }


    # Estimate Vec-Form of Phi_first_step
    inv.Phi_firststep <- MASS::ginv( as.matrix( t(sum_S_ZX_a %x% Diagonal(nof_dependent_vars)) %*% Lambda_inv_vec %*% (sum_S_ZX_a %x% Diagonal(nof_dependent_vars)) ))
    
    Phi_first_step_vec <- inv.Phi_firststep %*% t(sum_S_ZX_a %x% Diagonal(nof_dependent_vars)) %*% Lambda_inv_vec %*% matrix(t(sum_S_Zy_a), ncol=1)
    
    

    Phi_first_step <- matrix(Phi_first_step_vec,ncol=nof_dependent_vars*lags+nof_predet_vars+nof_exog_vars+
                               (system_instruments == TRUE)*(system_constant == TRUE)*1)

    colnames(Phi_first_step) <- c(transformed_lagged_vars, transformed_predet_vars, transformed_exog_vars,
                                  if(system_instruments== TRUE &  system_constant == TRUE){"const"} else{NULL})

    rownames(Phi_first_step) <- transformed_dependent_vars


    Delta_E_first_step <- mapply(function(i) delta_W[[i]] - delta_W_minus[[i]] %*% t(Phi_first_step), 1:nof_categories)
    Delta_E_first_step_NA <- mapply(function(i) delta_W_NA[[i]] - delta_W_minus_NA[[i]] %*% t(Phi_first_step), 1:nof_categories)
    
    # vectorize residuals
    e_hat <- mapply(function(i) matrix(t(Delta_E_first_step[[i]]), ncol=1), 1:nof_categories, SIMPLIFY = FALSE)

    #print("Before Z")
    Z <- mapply(function(i) Matrix( Q[[i]] %x% Diagonal(nof_dependent_vars), sparse = TRUE), 1:nof_categories)

    S_EQ_first_step <- mapply(function(i) Matrix(t(Delta_E_first_step[[i]]) %*% t(Q[[i]]), sparse = TRUE), 1:nof_categories, SIMPLIFY = FALSE)

    #print("Start sum_D_e")
    sum_D_e <- Reduce('+', mapply(function(i) {
      Matrix( tcrossprod( Z[[i]] %*% e_hat[[i]] ), sparse = TRUE )},
      1:nof_categories, SIMPLIFY = FALSE) )

    if (progressbar) pb$tick()


    # Twostep estimation -----------------------------------------------------

    if (steps == "twostep" | steps == "mstep") {


      if (min(abs(eigen(sum_D_e)$values)) < tol){

        #try sth different than the general inverse.
        D_e.inv <-  MASS::ginv(as.matrix(sum_D_e))

        #D_e.inv <- corpcor::pseudoinverse(sum_D_e)


        warning("The matrix D_e is singular, therefore the general inverse is used")
      } else{
        D_e.inv <-  solve(as.matrix(sum_D_e))

      }

      # Write the first part of the matrix computation for Phi in a seperate matrix in order
      # to use it again and to save m*(m+k) matrix inversion
      inv.Phi1 <- MASS::ginv( as.matrix( t(sum_S_ZX_a %x% Diagonal(nof_dependent_vars)) %*% D_e.inv %*% (sum_S_ZX_a %x% Diagonal(nof_dependent_vars)) ))

      Phi <- inv.Phi1 %*% t(sum_S_ZX_a %x% Diagonal(nof_dependent_vars)) %*% D_e.inv %*% matrix(t(sum_S_Zy_a), ncol=1)

      # secondstep coefficients

      Phi_second_step <- matrix(Phi,ncol=nof_dependent_vars*lags+nof_predet_vars+nof_exog_vars+
                                  (system_instruments == TRUE)*(system_constant == TRUE)*1)

      colnames(Phi_second_step) <- c(transformed_lagged_vars, transformed_predet_vars, transformed_exog_vars,
                                     if(system_instruments == TRUE & system_constant == TRUE){"const"} else{NULL})
      rownames(Phi_second_step) <- transformed_dependent_vars


      Delta_E_second_step <- mapply(function(i) delta_W[[i]] - delta_W_minus[[i]] %*% t(Phi_second_step), 1:nof_categories)
      Delta_E_second_step_NA <- mapply(function(i) delta_W_NA[[i]] - delta_W_minus_NA[[i]] %*% t(Phi_second_step), 1:nof_categories)


      sum_S_EQ <- Reduce('+', mapply(function(i) Matrix( t(Delta_E_second_step[[i]]) %*% t(Q[[i]]), sparse = TRUE), 1:nof_categories, SIMPLIFY = FALSE))

    }

    # M-Step Estimation: -----------------------------------------------------
    # Use repeat until coefficients converage:
    if (steps == "mstep"){

      # Before iterations. Define initial mstep values:
      Delta_E_mstep <- Delta_E_second_step

      # m--number of iterations.
      m <- 2
      Phi_mstep <- list()
      Phi_mstep[[1]] <- Phi_second_step
      convergence_crit <- list()
      convergence_crit[[1]] <- sqrt(sum((Phi_second_step-Phi_first_step)^2)/sum(Phi_first_step^2))
      repeat{
        e_hat <- mapply(function(i) matrix(t(Delta_E_mstep[[i]]), ncol=1), 1:nof_categories, SIMPLIFY = FALSE)

        S_EQ_mstep <- mapply(function(i) Matrix(t(Delta_E_mstep[[i]]) %*% t(Q[[i]]), sparse = TRUE), 1:nof_categories, SIMPLIFY = FALSE)

        #print("Start sum_D_e")
        sum_D_e <- Reduce('+', mapply(function(i) {
          Matrix( tcrossprod( Z[[i]] %*% e_hat[[i]] ), sparse = TRUE )},
          1:nof_categories, SIMPLIFY = FALSE) )


        if (min(abs(eigen(sum_D_e)$values)) < tol){
          D_e.inv <-  MASS::ginv(as.matrix(sum_D_e))

          warning("The matrix D_e is singular, therefore the general inverse is used")
        } else{
          D_e.inv <-  solve(as.matrix(sum_D_e))

        }

        # Write the first part of the matrix computation for Phi in a seperate matrix in order
        # to use it again and to save m*(m+k) matrix inversion
        inv.Phi1 <- MASS::ginv( as.matrix( t(sum_S_ZX_a %x% Diagonal(nof_dependent_vars)) %*% D_e.inv %*% (sum_S_ZX_a %x% Diagonal(nof_dependent_vars)) ))

        Phi <- inv.Phi1 %*% t(sum_S_ZX_a %x% Diagonal(nof_dependent_vars)) %*% D_e.inv %*% matrix(t(sum_S_Zy_a), ncol=1)

        # m-step: coefficients

        Phi_mstep[[m]] <- matrix(Phi,ncol=nof_dependent_vars*lags+nof_predet_vars+nof_exog_vars+
                                   (system_instruments == TRUE)*(system_constant == TRUE)*1)

        colnames(Phi_mstep[[m]]) <- c(transformed_lagged_vars, transformed_predet_vars, transformed_exog_vars,
                                      if(system_instruments == TRUE & system_constant == TRUE){"const"} else{NULL})
        rownames(Phi_mstep[[m]]) <- transformed_dependent_vars


        Delta_E_mstep <- mapply(function(i) delta_W[[i]] - delta_W_minus[[i]] %*% t(Phi_mstep[[m]]), 1:nof_categories)
        Delta_E_mstep_NA <- mapply(function(i) delta_W_NA[[i]] - delta_W_minus_NA[[i]] %*% t(Phi_mstep[[m]]), 1:nof_categories)


        sum_S_EQ <- Reduce('+', mapply(function(i) Matrix( t(Delta_E_mstep[[i]]) %*% t(Q[[i]]), sparse = TRUE), 1:nof_categories, SIMPLIFY = FALSE))


        convergence_crit[[m]] <- sqrt(sum((Phi_mstep[[m]]-Phi_mstep[[m-1]])^2)/sum(Phi_mstep[[m-1]]^2))

        # 1st Criterion to stop further iterations:
        if(all((convergence_crit[[m]])<0.00001)==T){break}

        # 2nd criterion coz of numerical issues with the inverse of second step residuals.
        if(convergence_crit[[m-1]] < convergence_crit[[m]]){
          # in this case, the m-1 iteration had a better fit.
          m <- m-1
          break
        }
        # if Phi does not converge then go to the next iteration.
        print(convergence_crit[[m]])
        m <- m+1
      }
    }


    # Robust standard errors --------------------------------------------------

    # Onestep

    var_first_step <- solve(t(sum_S_ZX_a %x% Diagonal(nof_dependent_vars)) %*% (Lambda_inv_vec) %*%
                              (sum_S_ZX_a %x% Diagonal(nof_dependent_vars))) %*% t(sum_S_ZX_a %x% Diagonal(nof_dependent_vars)) %*%
      (Lambda_inv_vec) %*% sum_D_e %*% (Lambda_inv_vec) %*%
      (sum_S_ZX_a %x% Diagonal(nof_dependent_vars)) %*% solve(t(sum_S_ZX_a %x% Diagonal(nof_dependent_vars)) %*%
                                                                        (Lambda_inv_vec) %*% (sum_S_ZX_a %x% Diagonal(nof_dependent_vars)))
    se_first_step <- matrix(sqrt(diag(var_first_step)), ncol=nof_dependent_vars*lags+nof_predet_vars+nof_exog_vars+
                              (system_instruments == TRUE)*(system_constant == TRUE)*1)
    z_values_first_step <- as.matrix(Phi_first_step/se_first_step)
    p_values_first_step <- 2*pnorm(-abs(z_values_first_step))

    # Twostep

    # Windmeijer finite-sample correction for the twostep covariance matrix
    if (steps == "twostep" | steps == "mstep") {
      if (progressbar) {
        pb <-
          progress::progress_bar$new(
            format = "Windmeijer - Sigmund robust se: [:bar] Iteration :current of :total" ,
            clear = TRUE, total = nof_categories)
      }
      # compute derivative matrices in order to compute the second step standard error

      # create a llist to write in the derivatives of the weighting matrix in order to perform the Windmeijer Correction for the standard errors
      deriv_W <- replicate(nof_dependent_vars*(nof_dependent_vars*lags+nof_predet_vars+nof_exog_vars+(system_instruments == TRUE)*(system_constant == TRUE)*1),
                           Matrix(0, nrow= nrow(Q[[1]]) * nof_dependent_vars, ncol = nrow(Q[[1]]) * nof_dependent_vars))

      # See Windmeijer_2005.pdf Page 33, between equation 3.2 and 3.3.
      # Derivative of W is explained there.
      # Make it faster:
      for(i in 1:nof_categories) {

        if (progressbar) pb$tick()

        mat_Z <- Matrix(S_ZX_a[[i]] %x% Diagonal(nof_dependent_vars), sparse = TRUE)

        if (steps == "twostep"){mat_Z2 <- Matrix(t(matrix(S_EQ_first_step[[i]], ncol = 1)), sparse = TRUE)}

        # Check if this line is mathematically correct. Nobody has ever applied the Windmeijer correction
        # to the m-step estimator!!
        if (steps == "mstep"){mat_Z2 <- Matrix(t(matrix(S_EQ_mstep[[i]], ncol = 1)), sparse = TRUE)}


        deriv_W <- mapply(function(j) deriv_W[[j]] - Matrix(( as.matrix(mat_Z[,j]) %*% mat_Z2 + t(as.matrix(mat_Z[,j]) %*% mat_Z2) ), sparse = TRUE),
                          1:(nof_dependent_vars*(nof_dependent_vars*lags+nof_predet_vars+nof_exog_vars+
                                                   (system_instruments == TRUE)*(system_constant == TRUE)*1)))

      }


      # first derivative of Windmeijers's Taylor expansion
      D_W <- matrix(nrow=nof_dependent_vars*(nof_dependent_vars*lags+nof_predet_vars+nof_exog_vars+
                                               (system_instruments == TRUE)*(system_constant == TRUE)*1),
                    ncol=nof_dependent_vars*(nof_dependent_vars*lags+nof_predet_vars+nof_exog_vars+
                                               (system_instruments == TRUE)*(system_constant == TRUE)*1))


      for(j in 1:ncol(D_W)){
        D_W[,j] <- as.vector(- inv.Phi1 %*% t(sum_S_ZX_a %x% Diagonal(nof_dependent_vars)) %*% D_e.inv %*% deriv_W[[j]] %*%  D_e.inv %*% matrix(sum_S_EQ, ncol=1))
      }

      if (steps == "twostep"){
        var_second_step <- inv.Phi1 + D_W %*% inv.Phi1 + inv.Phi1 %*% t(D_W) + D_W %*% var_first_step %*% t(D_W)
        se_second_step <- matrix(sqrt(diag(var_second_step)), ncol=nof_dependent_vars*lags+nof_predet_vars+nof_exog_vars+
                                   (system_instruments == TRUE)*(system_constant == TRUE)*1)
        z_values_second_step <- Phi_second_step/se_second_step
        p_values_second_step <- 2*pnorm(-abs(z_values_second_step))
      }

      if (steps == "mstep"){
        var_mstep <- inv.Phi1 + D_W %*% inv.Phi1 + inv.Phi1 %*% t(D_W) + D_W %*% var_first_step %*% t(D_W)
        se_mstep <- matrix(sqrt(diag(var_mstep)), ncol=nof_dependent_vars*lags+nof_predet_vars+nof_exog_vars+
                             (system_instruments == TRUE)*(system_constant == TRUE)*1)
        z_values_mstep <- Phi_mstep[[m]]/se_mstep
        p_values_mstep <- 2*pnorm(-abs(z_values_mstep))
      }


    }

    # Results -----------------------------------------------------------------
    # Add a few results for the Sargan-Hansen test
    # in pgmm -> sargan. In xtabond2 it is called hansen-test.

    inputargs <- list(dependent_vars = dependent_vars,
                      lags = lags,
                      transformation = transformation,
                      steps = steps,
                      system_instruments = system_instruments,
                      system_constant = system_constant,
                      collapse = collapse,
                      Set_Vars = Set_Vars,
                      Set_Vars_with_NAs = Set_Vars_with_NAs,
                      panel_identifier = panel_identifier,
                      max_instr_dependent_vars = max_instr_dependent_vars,
                      max_instr_predet_vars = max_instr_predet_vars,
                      min_instr_dependent_vars = min_instr_dependent_vars,
                      min_instr_predet_vars = min_instr_predet_vars,
                      tol = tol,
                      nof_observations = nof_observations,
                      obs_per_group_avg = obs_per_group_avg,
                      obs_per_group_min = obs_per_group_min,
                      obs_per_group_max = obs_per_group_max,
                      nof_groups = nof_groups
    )

    if(missing(predet_vars) == FALSE){inputargs$predet_vars = predet_vars}
    if(missing(exog_vars) == FALSE){inputargs$exog_vars = exog_vars}

    # Construct Output:
    if (steps == "onestep") {
      result <- list(Phi_first_step, se_first_step, p_values_first_step,
                     Delta_E_first_step_NA, Delta_E_first_step,  S_EQ_first_step, Lambda_inv_vec, Q, Z, delta_W,
                     delta_W_minus, V, position_exogenous)
      names(result) <- c("first_step", "standard_error_first_step", "p_values_first_step",
                         "residuals", "residuals_hansen_onestep", "Z_u_Hansen", "weighting_matrix_first_step", "unique_instruments", "instruments", "delta_W",
                         "delta_W_minus", "V", "position_exogenous")
    }

    if (steps == "twostep") {
      result <- list(Phi_first_step, Phi_second_step, se_first_step, se_second_step, p_values_first_step, p_values_second_step,
                     Delta_E_second_step_NA, Delta_E_second_step,  S_EQ_first_step, Lambda_inv_vec,  D_e.inv, Q, Z, delta_W,
                     delta_W_minus
      )
      names(result) <- c("first_step", "second_step", "standard_error_first_step", "standard_error_second_step", "p_values_first_step", "p_values_second_step",
                         "residuals", "residuals_hansen_onestep", "Z_u_Hansen", "weighting_matrix_first_step",  "weighting_matrix_second_step", "unique_instruments", "instruments", "delta_W",
                         "delta_W_minus"
      )
    }

    if (steps == "mstep"){

      result <- list(Phi_first_step, Phi_mstep[[m]], se_first_step, se_mstep, p_values_first_step, p_values_mstep,
                     Delta_E_mstep_NA, Delta_E_mstep,  S_EQ_mstep, Lambda_inv_vec,  D_e.inv, Q, Z, delta_W, 
                     delta_W_minus
      )
      names(result) <- c("first_step", "m_step", "standard_error_first_step", "standard_error_m_step", "p_values_first_step", "p_values_m_step",
                         "residuals", "residuals_hansen_onestep", "Z_u_Hansen", "weighting_matrix_first_step",  "weighting_matrix_m_step", "unique_instruments", "instruments", "delta_W",
                         "delta_W_minus"
      )


    }

    result <- c(result, inputargs)
    class(result) <- "pvargmm"
    result
  }



