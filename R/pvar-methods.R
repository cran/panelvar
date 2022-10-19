#' S3 Print Method for pvargamm
#' @param x object
#' @param ... further arguments
#' @method print pvargmm
#' @export

print.pvargmm <- function(x, ...) {
  res <- extract(x)
  print(texreg::screenreg(res, custom.model.names = names(res), digits = 4))
}

#' Knit Print Method for pvargmm
#' @param x object
#' @param ... further arguments
#' @method knit_print pvargmm
#' @export

knit_print.pvargmm <- function(x, ...) {
  res <- extract(x)
  knitr::asis_output(texreg::htmlreg(res, custom.model.names = names(res), digits = 4))
}

#' S3 Summary Method for pvargmm
#' @param object object
#' @param ... further arguments
#' @method summary pvargmm
#' @export
summary.pvargmm <- function(object, ...)
{
  x <- object
  sumry <- list()

  sumry$model_name <- paste("Dynamic Panel VAR estimation,",ifelse(x$steps=="twostep", "two-step",ifelse(x$steps=="onestep", "one-step","")), "GMM")
  sumry$transformation <-
    switch(x$transformation,
           fd = "First-differences",
           fod = "Forward orthogonal deviations"
           )


  sumry$groupvariable <- x$panel_identifier[1]
  sumry$timevariable <- x$panel_identifier[2]
  sumry$nof_observations <- x$nof_observations
  sumry$obs_per_group_min <- x$obs_per_group_min
  sumry$obs_per_group_avg <- x$obs_per_group_avg
  sumry$obs_per_group_max <- x$obs_per_group_max
  sumry$nof_groups <- x$nof_groups

  sumry$results <- extract(x)
  
  #sumry$transformation <- x$transformation
  sumry$exog_vars <- x$exog_vars
  
  sumry$min_instr_dependent_vars <- x$min_instr_dependent_vars
  sumry$max_instr_dependent_vars <- x$max_instr_dependent_vars
  
  if(is.null(x$predet_vars) == FALSE){sumry$min_instr_predet_vars = x$min_instr_predet_vars}
  if(is.null(x$predet_vars) == FALSE){sumry$max_instr_predet_vars = x$max_instr_predet_vars}
  if(is.null(x$predet_vars) == FALSE){sumry$predet_vars = x$predet_vars}
  
  sumry$collapse <- x$collapse
  
  sumry$instruments_standard <- ifelse(!is.null(x$exog_vars), paste0(toupper(x$transformation), ".(",paste0(x$exog_vars, collapse = " "),")"),"")
  
  #if(object$steps != "mstep")
  sumry$hansen_j_test <- hansen_j_test(object)

  class(sumry) <- "summary.pvargmm"
  sumry
}

#' S3 Print Method for summary.pvargmm
#' @param x object
#' @param ... further arguments
#' @method print summary.pvargmm
#' @export
print.summary.pvargmm <- function(x, ...) {
  cat("---------------------------------------------------\n")
  cat(x$model_name, "\n")
  cat("---------------------------------------------------\n")
  cat("Transformation:", x$transformation,"\n")
  cat("Group variable:", x$groupvariable,"\n")
  cat("Time variable:",x$timevariable,"\n")
  cat("Number of observations =", x$nof_observations,"\n")
  cat("Number of groups =", x$nof_groups, "\n")
  cat("Obs per group: min =",x$obs_per_group_min,"\n")
  cat("               avg =",x$obs_per_group_avg,"\n")
  cat("               max =",x$obs_per_group_max,"\n")
  cat("Number of instruments =", as.numeric(x$hansen_j_test$nof_instruments), "\n")
  #cat("---------------------------------------------------\n")
  print(texreg::screenreg(x$results, custom.model.names = names(x$results), digits = 4))
  cat("\n")
  cat("---------------------------------------------------\n")
  cat("Instruments for ", if(x$transformation=="fd"){"first differences"},if(x$transformation=="fod"){"orthogonal deviations"}," equation\n Standard\n", sep = "")
  cat(" ",x$instruments_standard)
  cat("\n")
  cat(" GMM-type\n")
  cat("  Dependent vars: L(", x$min_instr_dependent_vars, ", ", min(x$max_instr_dependent_vars, x$obs_per_group_max),")",  sep="")
  cat("\n")
  if(is.null(x$predet_vars) == FALSE){cat("  Predet vars: L(", x$min_instr_predet_vars, ", ",min(x$max_instr_predet_vars, x$obs_per_group_max),")\n", sep="")}
  cat("  Collapse = ", as.character(x$collapse),"\n")
  cat("---------------------------------------------------\n")
  cat("\n")
  cat("Hansen test of overid. restrictions: chi2(",as.numeric(x$hansen_j_test$parameter),") = ",round(as.numeric(x$hansen_j_test$statistic),2)," Prob > chi2 = ", round(as.numeric(x$hansen_j_test$p.value),3),"\n",sep = "")
  cat("(Robust, but weakened by many instruments.)")
  #invisible(x)
}

#' Knit Print summary Method
#' @param x object
#' @param ... further arguments
#' @method knit_print summary.pvargmm
#' @export
# @importFrom knitr knit_print
knit_print.summary.pvargmm <- function(x, ...) {
  res <- paste0(c("<p><b>",x$model_name,"</b></p>",
                  "<p>Transformation: <em>",x$transformation,"</em><br>",
                  "Group variable: <em>",x$groupvariable,"</em><br>",
                  "Time variable: <em>",x$timevariable,"</em><br>",
                  "Number of observations = <em>",x$nof_observations,"</em><br>",
                  "Number of groups = <em>",x$nof_groups,"</em><br>",
                  "Obs per group: min = <em>",x$obs_per_group_min,"</em><br>",
                  "Obs per group: avg = <em>",x$obs_per_group_avg,"</em><br>",
                  "Obs per group: max = <em>",x$obs_per_group_max,"</em><br>",
                  "Number of instruments = <em>",as.numeric(x$hansen_j_test$nof_instruments),"</em><br>",
"</p>",
    texreg::htmlreg(x$results, custom.model.names = names(x$results), digits = 4, caption = "", center = FALSE),
"<p><b>", "Instruments for ", if(x$transformation=="fd"){"first differences"},if(x$transformation=="fod"){"orthogonal deviations"}," equation </b><p>Standard<br><em>",paste0(x$exog_vars, collapse = ", "),"</em></p>",
"<p>GMM-type<br><em>","Dependent vars: L(", x$min_instr_dependent_vars, ",",min(x$max_instr_dependent_vars, x$obs_per_group_max),")", 
if(is.null(x$predet_vars) == FALSE){c("<br>Predet vars: L(", x$min_instr_predet_vars, ", ",min(x$max_instr_predet_vars, x$obs_per_group_max))},")<br>Collapse = ",as.character(x$collapse),"</em></p>",
"<p><b>Hansen test of overid. restrictions:</b> <em>chi2(",as.numeric(x$hansen_j_test$parameter),") = ",round(as.numeric(x$hansen_j_test$statistic),2)," Prob > chi2 = ", round(as.numeric(x$hansen_j_test$p.value),3),"<br>",
"(Robust, but weakened by many instruments.)","</em><br>",
"</p>"
)
    )
  knitr::asis_output(res)
}


#' Extract PVAR(p) Model Coefficients
#' @param object object
#' @param ... further arguments
#' @export
#' @examples
#' data("ex1_dahlberg_data")
#' coef(ex1_dahlberg_data)

coef.pvargmm <- function(object, ...) {
  if(object$steps == "onestep") {object$first_step
  } else if (object$steps == "twostep") {
    object$second_step
  }
  else object$m_step
}

#' Standard Error S3 Method
#' @param object Object
#' @param ... Further arguments
#' @export
#' @examples
#' data("ex1_dahlberg_data")
#' se(ex1_dahlberg_data)

se <- function(object, ...) UseMethod("se")

#' @rdname se
#' @export
se.pvargmm <- function(object, ...) {
  if(object$steps == "onestep") {object$standard_error_first_step
    } else if (object$steps == "twostep") {
      object$standard_error_second_step
    }
  else object$standard_error_m_step
}

#' P-value S3 Method
#' @param object Object
#' @param ... Further arguments
#' @export
#' @examples 
#' data("ex1_dahlberg_data")
#' pvalue(ex1_dahlberg_data)

pvalue <- function(object, ...) UseMethod("pvalue")

#' @rdname pvalue
#' @export
pvalue.pvargmm <- function(object, ...) {
  if(object$steps == "onestep") {object$p_values_first_step
  } else if (object$steps == "twostep") {
    object$p_values_second_step
  }
  else object$p_values_m_step
}

#' Extract Coefficients and GOF Measures from a Statistical Object
#' @param model Model
#' @param ... Further arguments passed to or from other methods
#' @export
#' @examples 
#' data("ex1_dahlberg_data")
#' extract(ex1_dahlberg_data)

extract <- function(model, ...) UseMethod("extract")

#' @rdname extract
#' @export

extract.pvargmm <- function(model, ...) {
  co.m <- coef(model)
  modelnames <- row.names(co.m)
  # remove fod_ and fd_
  modelnames <- substr(modelnames,nchar(model$transformation)+2,max(nchar(modelnames)))
  se.m <- se(model)
  pval.m <- pvalue(model)
  equationList <- list()
  coef.names <- colnames(co.m)
  
  # remove fd_ and fod_, handel possible system constant
  if(model$system_instruments & model$system_constant) {
    coef.names <- c(substr(coef.names[1:(length(coef.names)-1)],nchar(model$transformation)+2,max(nchar(coef.names))),"const")
  } else {
    coef.names <- substr(coef.names,nchar(model$transformation)+2,max(nchar(coef.names)))
  }
  
  for (eq in 1:length(model$dependent_vars)) {
    tr <- texreg::createTexreg(
      coef.names = coef.names, 
      coef = as.vector(co.m[eq,]),
      se = as.vector(se.m[eq,]),
      pvalues = as.vector(pval.m[eq,]),
      model.name = modelnames[eq]
    )
    equationList[[eq]] <- tr
  }
  names(equationList) <- modelnames
  equationList
}
#'
#'
#' Extracting Fixed Effects
#' @param model Model
#' @param Only_Non_NA_rows Filter NA rows
#' @param ... Further arguments passed to or from other methods
#' @export
#' @examples
#' data("ex1_dahlberg_data")
#' fixedeffects(ex1_dahlberg_data)

fixedeffects <- function(model, ...) UseMethod("fixedeffects")

#' @rdname fixedeffects
#' @export


fixedeffects.pvargmm <- function(model, Only_Non_NA_rows = TRUE, ...){

  #model <-  test_Q_for_endo
  lags <- c(1:model$lags)

  # Determine LHS and RHS variables that for which afterwards the mean is taken.
  Y_mean <- c(model$dependent_vars)

  # Improve next line:
  X_mean_1 <- unique(matrix(mapply(function(i) mapply(function(j)
    paste("lag", lags[j], "_", model$dependent_vars, sep = ""),  1:length(lags)), 1:length(model$dependent_vars)), ncol = 1))

  X_mean_2 <- c(if(!is.null(model$predet_vars)){model$predet_vars}, if(!is.null(model$exog_vars)){model$exog_vars})

  X_mean <- c(X_mean_1, if(!is.null(X_mean_2)){X_mean_2}, if(model$system_instruments==TRUE){"const"})

  # Add a column with 1 if system instruments are TRUE.
  if (model$system_instruments==TRUE){
    model$Set_Vars$const <- 1
  }

  Set_Vars_reduced <- subset(model$Set_Vars, select = c("category", "period" , Y_mean, X_mean))

  # Two options:
  # 1. Set each row to NA if there is one NA element
  # 2. Calculate mean based on what is there for each variable.
  # STATA uses Option 1. which is implemented for xtreg then type predict somename, u.
  if (Only_Non_NA_rows == TRUE){

    for (i0 in 1:dim(Set_Vars_reduced)[1]){

      if (is.na(rowSums(Set_Vars_reduced[i0,3:ncol(Set_Vars_reduced)]))){

        Set_Vars_reduced[i0,3:ncol(Set_Vars_reduced)] <- NA

      }

    }

  }

  if (model$steps == "onestep"){model$beta <- model$first_step}
  if (model$steps == "twostep"){model$beta <- model$second_step}


  Unique_categories <- unique(Set_Vars_reduced$category)

  # Prepare Set_Vars from the model:
  fixef_pvargmm <- list() #matrix(NA, nrow = length(model$dependent_vars), ncol = length(unique(model$Set_Vars$category)))

  for (i0 in 1:length(unique(model$Set_Vars$category))){

    # Calculation based on Wooldridge Introduction to Econometrics. Part3 Chapter 14 Page 486 equation 14.6
    fixef_pvargmm[[i0]] <- as.matrix(mapply(function(i) mean(Set_Vars_reduced[Set_Vars_reduced$category ==  Unique_categories[i0],
                                                                           model$dependent_vars[i]], na.rm = TRUE), 1:length(model$dependent_vars)) ) -
      model$beta %*% as.matrix( mapply(function(i) mean(Set_Vars_reduced[Set_Vars_reduced$category ==  Unique_categories[i0], X_mean[i]], na.rm = TRUE), 1:length(X_mean)))

  }

  return(fixef_pvargmm)

}

#' Extracting Level Residuals
#' @param model Model
#' @param ... Further arguments passed to or from other methods
#' @export
#' @examples
#' data("ex1_dahlberg_data")
#' residuals_level(ex1_dahlberg_data)

residuals_level <- function(model, ...) UseMethod("residuals_level")

#' @rdname residuals_level
#' @export

residuals_level.pvargmm <- function(model, ...){
  
  #model <-  test_Q_for_endo
  lags <- c(1:model$lags)
  
  Y_mean <- c(model$dependent_vars)
  
  # Improve next line:
  X_mean_1 <- unique(matrix(mapply(function(i) mapply(function(j)
    paste("lag", lags[j], "_", model$dependent_vars, sep = ""),  1:length(lags)), 1:length(model$dependent_vars)), ncol = 1))
  
  X_mean_2 <- c(if(!is.null(model$predet_vars)){model$predet_vars}, if(!is.null(model$exog_vars)){model$exog_vars})
  
  X_mean <- c(X_mean_1, if(!is.null(X_mean_2)){X_mean_2}, if(model$system_instruments==TRUE){"const"})
  
  # Add a column with 1 if system instruments are TRUE.
  if (model$system_instruments==TRUE){
    model$Set_Vars$const <- 1
  }
  
  Set_Vars_reduced <- subset(model$Set_Vars, select = c("category", "period" , Y_mean, X_mean))
  
  
  if (model$steps == "onestep"){Beta <- model$first_step}
  if (model$steps == "twostep"){Beta <- model$second_step}
  
  # Get fixed effect:
  fixed_eff_list <- fixedeffects(model)
  
  # Calculate vector of residuals for each panel category separately:
  
  list.residuals <- list()
  
  vector.categories <- unique(Set_Vars_reduced$category)
  
  for (i0 in 1:length(vector.categories)){
    
    # We accept NAs in the first p-lags:
    # Get data for i
    Set_Vars_reduced.i <- subset(Set_Vars_reduced, 
                                 Set_Vars_reduced$category == vector.categories[i0])
    
    # Get fixed effects for i
    fixed_eff_list.i <- fixed_eff_list[[i0]]
    
    # Dimension: [nof_dependent] x [T]
    big.Y <- t(Set_Vars_reduced.i[,Y_mean])
    
    # Dimension: [nof_dependent] x [nof RHS variables plus 1 row for the fixed effect]
    Big.Beta <- cbind(fixed_eff_list.i, Beta)
    
    # Dimension: [nof RHS variables plus 1 row for the fixed effect] x [T]
    Big.X <- rbind(1, t(Set_Vars_reduced.i[,X_mean]))
    
    # Dimension: [nof RHS variables] x [T]
    list.residuals[[i0]] <- big.Y - Big.Beta %*% Big.X
    colnames(list.residuals[[i0]]) <- Set_Vars_reduced.i$period
  }
  
  categories <- unique(Set_Vars_reduced$category)
  # Name the list entries with the category identifier:
  names(list.residuals) <- categories
  
  # Make a table with first two columns as category, period similar to Set_Vars
  df_residuals <- mapply(function(i) t(list.residuals[[i]]), 1:length(list.residuals), SIMPLIFY = FALSE)
  
  df_residuals <- do.call("rbind", mapply(function(i) data.frame(category = categories[i], 
                                                                 period = row.names(df_residuals[[i]]), df_residuals[[i]]), 1:length(df_residuals), SIMPLIFY = FALSE))
  row.names(df_residuals) <- NULL
  
  return(df_residuals)  
}  

#' @export
vcov.pvargmm <- function(object, ...){

  pvar_model <- object

  residuals <- pvar_model$residuals[[1]]

  for (i0 in 2:length(pvar_model$residuals) ){

    residuals <- rbind(residuals, pvar_model$residuals[[i0]])

  }

  residuals <- na.exclude( as.matrix(residuals) )

  Sigma_Hat1 <- ( t(residuals) %*% residuals ) / (nrow(residuals) - ncol(pvar_model$first_step))


  return(Sigma_Hat1)
}


#' Sargan-Hansen-J-Test for Overidentification
#' @param model A PVAR model
#' @param ... Further arguments passed to or from other methods
#' @export
#' @examples 
#' data("ex1_dahlberg_data")
#' hansen_j_test(ex1_dahlberg_data)
hansen_j_test <- function(model, ...) UseMethod("hansen_j_test")

#' @rdname hansen_j_test
#'@export
hansen_j_test.pvargmm <- function(model, ...){

  steps <- model$steps
  Z_u <- 0

  # Calculate Sum Z %*% u (Instruments times residuals of the first step, where NA are still set to 0)
  for (i0 in 1:length(model$unique_instruments)){

    Z_u <- Z_u + as.matrix(model$instruments[[i0]]) %*% matrix(t(as.matrix(model$residuals_hansen_onestep[[i0]])), ncol=1)

  }

  # Hansen J test for the first step makes no sense. But still is mathematically possible.
  if(steps == c("onestep")){

    system_over_id <- t(Z_u) %*% ( model$weighting_matrix_first_step ) %*% Z_u
    warning("Hansen J test for the first step makes no sense. But still is mathematically possible")
    # Number of parameters:
    Ktot_fdgmm <- ncol(model$first_step) * nrow(model$first_step)

  }

  # Hansen J test for the second step.
  if (steps == c("twostep")){

    system_over_id <- t(Z_u) %*% (model$weighting_matrix_second_step) %*% Z_u

    # Number of parameters:
    Ktot_fdgmm <- ncol(model$second_step) * nrow(model$second_step)

  }

  # Hansen J test for the m step.
  if (steps == c("mstep")){

    system_over_id <- t(Z_u) %*% (model$weighting_matrix_m_step) %*% Z_u

    # Number of parameters:
    Ktot_fdgmm <- ncol(model$m_step) * nrow(model$m_step)

  }

  # Number of instruments:
  p_fdgmm <- nrow(model$unique_instruments[[1]])

  # alternative:
  p_fdgmm <- nrow(model$instruments[[1]])

  # Number of overidentification restrictions:
  parameter_fdgmm <- (p_fdgmm - Ktot_fdgmm)
  nof_instruments <- p_fdgmm 
  
  pval_fdgmm <-  pchisq(as.numeric(system_over_id), df = parameter_fdgmm, lower.tail = FALSE)

  #browser()
  hansen_j_test <- list(statistic = as.numeric(system_over_id), p.value = pval_fdgmm, parameter = parameter_fdgmm,
                        nof_instruments = nof_instruments,
                        method = "Hansen-J-Test")

  return(hansen_j_test)


}


#' Orthogonal Impulse Response Function
#'
#' @param model A PVAR model
#' @param n.ahead Any stable AR() model has an infinite MA representation. Hence any shock can be simulated infinitely into the future. For each forecast step t you need an addtional MA term.
#' @export
oirf <- function(model, n.ahead) UseMethod("oirf")

#' @export
oirf.pvargmm <- function(model, n.ahead){
# check if model is one-dimensional => error

  if (model$steps == c("onestep")){

    Phi <- model$first_step

  }

  if (model$steps == c("twostep")){

    Phi <- model$second_step

  }

  # P <- t(chol(covm))
  MA_Phi <- ma_phi_representation(Phi = Phi,
                                  ma_approx_steps = n.ahead,
                                  lags = model$lags)

  Sigma_Hat1 <- vcov(model)

  P <- chol(Sigma_Hat1)

  irf_output <- list()
  MA_Phi_P <- list()

  # Calculate 2.3.27 in Luetkepohl p. 58.
  for (i0 in 1:n.ahead){

    MA_Phi_P[[i0]] <- MA_Phi[[i0]] %*% t(P)

  }

  # Redistribute the MA_Phi_P into the correct list form.
  for (i0 in 1:nrow(Phi)){

    zwischen <- matrix(NA, nrow = length(MA_Phi_P), ncol = nrow(Phi))

    for (i1 in 1:length(MA_Phi_P)){

      zwischen[i1,] <- MA_Phi_P[[i1]][,i0]

    }

    irf_output[[i0]] <- zwischen
    colnames(irf_output[[i0]]) <- model$dependent_vars

  }
  names(irf_output) <- model$dependent_vars
  class(irf_output) <- "pvaroirf"
  return(irf_output)

}

#' @export
plot.pvaroirf <- function(x, cibootstrap_results, ...) {
  dependent_vars <- names(x)
  steps <- nrow(x[[1]])

  plotdata <- do.call(rbind,lapply(x, reshape2::melt))
  names(plotdata) <- c("periods", "impulse", "value")
  plotdata$impulse <-  paste(rep(dependent_vars, each = steps*length(dependent_vars)), "on", plotdata$impulse)

  p <-
    ggplot2::ggplot(plotdata, ggplot2::aes_(x = ~periods, y = ~value)) +
    ggplot2::geom_line(colour = "orangered4") + ggplot2::facet_wrap(~impulse) +
    ggplot2::labs(title = "Orthogonalized impulse response function") +
    ggplot2::xlab("steps") + ggplot2::ylab("") +
    ggplot2::scale_x_continuous(breaks = pretty(plotdata$periods)) +
    ggplot2::theme_bw()

  if (!missing(cibootstrap_results)) {
    ci_data_upper <-
      data.frame(do.call(rbind,lapply(cibootstrap_results$Upper, reshape2::melt)), Bound = "Upper")

    names(ci_data_upper)[1:3] <- c("periods", "impulse", "value")
    row.names(ci_data_upper) <- NULL
    ci_data_upper$impulse <- paste(rep(dependent_vars,
                                       each = steps*length(dependent_vars)), "on", ci_data_upper$impulse)

    ci_data_lower <-
      data.frame(do.call(rbind,lapply(cibootstrap_results$Lower, reshape2::melt)), Bound = "Lower")

    names(ci_data_lower)[1:3] <- c("periods", "impulse", "value")
    row.names(ci_data_lower) <- NULL
    ci_data_lower$impulse <- paste(rep(dependent_vars,
                                       each = steps*length(dependent_vars)), "on", ci_data_lower$impulse)

    ci_data <- reshape2::dcast(rbind(ci_data_upper, ci_data_lower), periods + impulse ~ Bound, value.var = "value")

    plotdata <-
      merge(plotdata, ci_data, all.x = TRUE, by = c("periods", "impulse"), sort = FALSE)

    p <-
      ggplot2::ggplot(plotdata, ggplot2::aes_(x = ~periods, y = ~value)) +
      ggplot2::geom_line(colour = "orangered4") +
      ggplot2::geom_ribbon(ggplot2::aes_(x = ~periods,ymin = ~Lower, ymax = ~Upper), alpha = 0.3,fill = "steelblue2") +
      ggplot2::facet_wrap(~impulse) +
      ggplot2::labs(title = "Orthogonalized impulse response function") +
      ggplot2::xlab("steps") + ggplot2::ylab("") +
      ggplot2::scale_x_continuous(breaks = pretty(plotdata$periods)) +
      ggplot2::theme_bw() +
      ggplot2::labs(subtitle = paste0("OIRF and ",cibootstrap_results$CI*100,"%", " confidence bands"))
  }
  p
}


#' Generalized Impulse Response Function
#'
#' @param model A PVAR model
#' @param n.ahead Any stable AR() model has an infinite MA representation. Hence any shock can be simulated infinitely into the future. For each forecast step t you need an addtional MA term.
#' @param ma_approx_steps MA approximation steps
#'
#' @export
#' @examples
#' data("ex1_dahlberg_data")
#' girf(ex1_dahlberg_data, n.ahead = 8, ma_approx_steps= 8)
girf <- function(model, n.ahead, ma_approx_steps) UseMethod("girf")


#' @rdname girf
#' @export

girf.pvargmm <- function(model, n.ahead, ma_approx_steps) {

  if (model$steps == c("onestep")){

    Phi <- model$first_step

  }

  if (model$steps == c("twostep")){

    Phi<-model$second_step

  }

  # Phi is a function for the MA representaton in the package vars.
  MA_Phi <- ma_phi_representation(Phi = Phi,
                                  ma_approx_steps = ma_approx_steps,
                                  lags = model$lags)

  Sigma_Hat1 <- vcov(model)

  sigmas <- sqrt(diag(Sigma_Hat1))

  girf_output_time <- list()

  for (i0 in 1:n.ahead){

    girf_output_time[[i0]] <-  t( MA_Phi[[i0]] %*% Sigma_Hat1)  /  sigmas

  }

  girf_output_shock <- list()

  for (i0 in 1:nrow(Sigma_Hat1)){

    zwischen <- matrix(NA, nrow = n.ahead, ncol = nrow(Sigma_Hat1))
    for (i1 in 1:length(girf_output_time)){

      zwischen[i1,] <- girf_output_time[[i1]][i0,]

    }

    girf_output_shock[[i0]] <- zwischen
    colnames(girf_output_shock[[i0]]) <- model$dependent_vars
  }

names(girf_output_shock) <- model$dependent_vars

  class(girf_output_shock) <- "pvargirf"
  return(girf_output_shock)
}


#' @export
plot.pvargirf <- function(x, cibootstrap_results, ...) {
  dependent_vars <- names(x)
  steps <- nrow(x[[1]])

  plotdata <- do.call(rbind,lapply(x, reshape2::melt))
  names(plotdata) <- c("periods", "impulse", "value")
  plotdata$impulse <-  paste(rep(dependent_vars, each = steps*length(dependent_vars)), "on", plotdata$impulse)

  p <-
    ggplot2::ggplot(plotdata, ggplot2::aes_(x = ~periods, y = ~value)) +
    ggplot2::geom_line(colour = "orangered4") + ggplot2::facet_wrap(~impulse) +
    ggplot2::labs(title = "Generalized impulse response function") +
    ggplot2::xlab("steps") + ggplot2::ylab("") +
    ggplot2::scale_x_continuous(breaks = pretty(plotdata$periods)) +
    ggplot2::theme_bw()

  if (!missing(cibootstrap_results)) {
    ci_data_upper <-
      data.frame(do.call(rbind,lapply(cibootstrap_results$Upper, reshape2::melt)), Bound = "Upper")

    names(ci_data_upper)[1:3] <- c("periods", "impulse", "value")
    row.names(ci_data_upper) <- NULL
    ci_data_upper$impulse <- paste(rep(dependent_vars,
                                       each = steps*length(dependent_vars)), "on", ci_data_upper$impulse)

    ci_data_lower <-
      data.frame(do.call(rbind,lapply(cibootstrap_results$Lower, reshape2::melt)), Bound = "Lower")

    names(ci_data_lower)[1:3] <- c("periods", "impulse", "value")
    row.names(ci_data_lower) <- NULL
    ci_data_lower$impulse <- paste(rep(dependent_vars,
                                       each = steps*length(dependent_vars)), "on", ci_data_lower$impulse)

    ci_data <- reshape2::dcast(rbind(ci_data_upper, ci_data_lower), periods + impulse ~ Bound, value.var = "value")

    plotdata <-
      merge(plotdata, ci_data, all.x = TRUE, by = c("periods", "impulse"), sort = FALSE)

    p <-
      ggplot2::ggplot(plotdata, ggplot2::aes_(x = ~periods, y = ~value)) +
      ggplot2::geom_line(colour = "orangered4") +
      ggplot2::geom_ribbon(ggplot2::aes_(x = ~periods,ymin = ~Lower, ymax = ~Upper), alpha = 0.3,fill = "steelblue2") +
      ggplot2::facet_wrap(~impulse) +
      ggplot2::labs(title = "Generalized impulse response function") +
      ggplot2::xlab("steps") + ggplot2::ylab("") +
      ggplot2::scale_x_continuous(breaks = pretty(plotdata$periods)) +
      ggplot2::theme_bw() +
      ggplot2::labs(subtitle = paste0("GIRF and ",cibootstrap_results$CI*100,"%", " confidence bands"))
  }
  p
}

#' Stability of PVAR(p) model
#' @param model PVAR model
#' @param ... Further arguments
#' @return A \code{pvarstability} object containing eigenvalue stability conditions
#' 
#' @export
#' @examples 
#' data("ex1_dahlberg_data")
#' stability_info <- stability(ex1_dahlberg_data)
#' print(stability_info)
#' plot(stability_info)

stability <- function(model, ...) UseMethod("stability")

#' @rdname stability
#' @export
stability.pvargmm <- function(model, ...) {

  if(model$steps == "twostep") {
    Phi <- model$second_step
  }
  if(model$steps == "onestep") {
    Phi <- model$first_step
  }
  if(model$steps == "mstep") {
    Phi <- model$m_step
  }

  # Get the PVAR(1) represenation
  A_full <- pvar1_phi_representation(Phi = Phi, lags = model$lags)

  # Calculate the eigenvalues of the matrix A_full.
  A_full_eigen <- eigen(A_full, only.values = TRUE)

  res <- data.frame(Eigenvalue = A_full_eigen$values, Modulus = abs(A_full_eigen$values))
  class(res) <- c("pvarstability", "data.frame")
  return(res)
}

#' S3 print method for pvarstability object
#' @param x object
#' @param ... further arguments
#' @method print pvarstability
#' @export
print.pvarstability <- function(x, ...) {
  # TODO: html version for markdown
  cat("Eigenvalue stability condition:\n\n")
  print.data.frame(x)

  if(sum(x$Modulus>1)==0) {
      cat("\nAll the eigenvalues lie inside the unit circle.\n")
      cat("PVAR satisfies stability condition.\n")
  } else {
    cat("Not all the eigenvalues lie inside the unit circle.\n")
  cat("PVAR does not satisfies stability condition.\n")
  }
}

#' S3 plot method for pvarstability object, returns a \code{ggplot} object
#' @param x object
#' @param ... further arguments
#' @export
plot.pvarstability <- function(x, ...) {
  x1 <- seq(-1,1,length=1000)
  y1 <- sqrt(1 - x1^2)
  y2 <- -sqrt(1 - x1^2)

  dfroots <- data.frame(Real=Re(x$Eigenvalue),
                        Imaginary=Im(x$Eigenvalue))

  limits <- max(abs(dfroots))

  if (limits < 1) {
    limits <- c(-1,1)
  } else {
    limits <- c(-ceiling(limits), ceiling(limits))
  }

  ggplot2::ggplot(data.frame(x=x1,y1=y1,y2=y2)) + #Using the results from earlier
    ggplot2::geom_hline(ggplot2::aes(yintercept=0)) +
    ggplot2::geom_vline(ggplot2::aes(xintercept=0)) +
    ggplot2::geom_line(ggplot2::aes(x=x,y=y1)) + # top of circle
    ggplot2::geom_line(ggplot2::aes(x=x,y=y2)) +
    ggplot2::geom_point(ggplot2::aes_(x=~Real, y=~Imaginary),
               dfroots, size = 2) +
    ggplot2::scale_y_continuous("Imaginary", limits=c(-1,1)) +
    ggplot2::scale_x_continuous("Real", limits=c(-1,1)) + ggplot2::ggtitle("Roots of the companion matrix")

  }
