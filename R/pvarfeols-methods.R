#' Extract PVARFEOLS(p) Model Coefficients
#' @param object object
#' @param ... further arguments
#' @export



coef.pvarfeols <- function(object, ...) {
  object$OLS$coef
}

#' @rdname se
#' @export
se.pvarfeols <- function(object, ...) {
  object$OLS$se
}

#' @rdname pvalue
#' @export
pvalue.pvarfeols <- function(object, ...) {
  object$OLS$pvalues
}

#' @rdname extract
#' @export

extract.pvarfeols <- function(model, ...) {
  co.m <- coef(model)
  modelnames <- row.names(co.m)
  se.m <- se(model)
  pval.m <- pvalue(model)
  equationList <- list()
  for (eq in 1:length(model$dependent_vars)) {
    tr <- texreg::createTexreg(
      coef.names = colnames(co.m),
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

#' S3 Print Method for pvarfeols
#' @param x object
#' @param ... further arguments
#' @method print pvarfeols
#' @export

print.pvarfeols <- function(x, ...) {
  res <- extract(x)
  print(texreg::screenreg(res, custom.model.names = names(res), digits = 4))
}

#' Knit Print Method for pvarfeols
#' @param x object
#' @param ... further arguments
#' @export

knit_print.pvarfeols <- function(x, ...) {
  res <- extract(x)
  knitr::asis_output(texreg::htmlreg(res, custom.model.names = names(res), digits = 4))
}

#' S3 Summary Method for pvarfeols
#' @param object object
#' @param ... further arguments
#' @method summary pvarfeols
#' @export
summary.pvarfeols <- function(object, ...)
{
  x <- object
  sumry <- list()
  
  sumry$model_name <- "Fixed Effects OLS Panel VAR estimation"
  sumry$transformation <- x$transformation
  
  
  sumry$groupvariable <- x$panel_identifier[1]
  sumry$timevariable <- x$panel_identifier[2]
  sumry$nof_observations <- x$nof_observations
  sumry$obs_per_group_min <- x$obs_per_group_min
  sumry$obs_per_group_avg <- x$obs_per_group_avg
  sumry$obs_per_group_max <- x$obs_per_group_max
  sumry$nof_groups = x$nof_groups
  
  sumry$results <- extract(x)
  
  class(sumry) <- "summary.pvarfeols"
  sumry
}

#' S3 Print Method for summary.pvarfeols
#' @param x object
#' @param ... further arguments
#' @method print summary.pvarfeols
#' @export
print.summary.pvarfeols <- function(x, ...) {
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
  #cat("---------------------------------------------------\n")
  print(texreg::screenreg(x$results, custom.model.names = names(x$results), digits = 4))
  cat("\n")
}

#' Knit Print summary Method
#' @param x object
#' @param ... further arguments
#' @export
knit_print.summary.pvarfeols <- function(x, ...) {
  res <- paste0(c("<p><b>",x$model_name,"</b></p>",
                  "<p>Transformation: <em>",x$transformation,"</em><br>",
                  "Group variable: <em>",x$groupvariable,"</em><br>",
                  "Time variable: <em>",x$timevariable,"</em><br>",
                  "Number of observations = <em>",x$nof_observations,"</em><br>",
                  "Number of groups = <em>",x$nof_groups,"</em><br>",
                  "Obs per group: min = <em>",x$obs_per_group_min,"</em><br>",
                  "Obs per group: avg = <em>",x$obs_per_group_avg,"</em><br>",
                  "Obs per group: max = <em>",x$obs_per_group_max,"</em><br>",
                  "</p>",
                  texreg::htmlreg(x$results, custom.model.names = names(x$results), digits = 4, caption = "", center = FALSE),
                  "</p>"
  )
  )
  knitr::asis_output(res)
}

#' @export
vcov.pvarfeols <- function(object, ...){
  
  residuals <- na.exclude( as.matrix(object$residuals) )
  
  Sigma_Hat1 <- ( t(residuals) %*% residuals ) / (nrow(residuals) - ncol(object$OLS$coef))

  return(Sigma_Hat1)
}

#' @export
oirf.pvarfeols <- function(model, n.ahead){

  Phi = coef(model)
  
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
  class(irf_output) <- "pvarfeolsoirf"
  return(irf_output)
  
}

#' @export
plot.pvarfeolsoirf <- function(x, cibootstrap_results, ...) {
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


#' @export
girf.pvarfeols <- function(model, n.ahead, ma_approx_steps) {
  
  Phi = coef(model)
  
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
  
  class(girf_output_shock) <- "pvarfeolsgirf"
  return(girf_output_shock)
}


#' @export
plot.pvarfeolsgirf <- function(x, cibootstrap_results, ...) {
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

#' @rdname stability
#' @export
stability.pvarfeols <- function(model, ...) {
  
  Phi <- coef(model)
  
  # Get the PVAR(1) represenation
  A_full <- pvar1_phi_representation(Phi = Phi, lags = model$lags)
  
  # Calculate the eigenvalues of the matrix A_full.
  A_full_eigen <- eigen(A_full, only.values = TRUE)
  
  res <- data.frame(Eigenvalue = A_full_eigen$values, Modulus = abs(A_full_eigen$values))
  class(res) <- c("pvarstability", "data.frame")
  return(res)
}
