#' Extract PVARHK(p) Model Coefficients
#' @param object object
#' @param ... further arguments
#' @export

coef.pvarhk <- function(object, ...) {
  object$HK$coef
}

#' @rdname se
#' @export
se.pvarhk <- function(object, ...) {
  object$HK$se
}

#' @rdname pvalue
#' @export
pvalue.pvarhk <- function(object, ...) {
  object$HK$pvalues
}


#' @rdname extract
#' @export

extract.pvarhk <- function(model, ...) {
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

#' S3 Print Method for pvarhk
#' @param x object
#' @param ... further arguments
#' @method print pvarhk
#' @export

print.pvarhk <- function(x, ...) {
  res <- extract(x)
  print(texreg::screenreg(res, custom.model.names = names(res), digits = 4))
}

#' Knit Print Method for pvarhk
#' @param x object
#' @param ... further arguments
#' @method knit_print pvarhk
#' @export

knit_print.pvarhk <- function(x, ...) {
  res <- extract(x)
  knitr::asis_output(texreg::htmlreg(res, custom.model.names = names(res), digits = 4))
}

#' S3 Summary Method for pvarhk
#' @param object object
#' @param ... further arguments
#' @method summary pvarhk
#' @export
summary.pvarhk <- function(object, ...)
{
  x <- object
  sumry <- list()
  
  sumry$model_name <- "Hahn Kuehrsteiner Panel VAR estimation"
  sumry$transformation <- x$transformation
  
  
  sumry$groupvariable <- x$panel_identifier[1]
  sumry$timevariable <- x$panel_identifier[2]
  sumry$nof_observations <- x$nof_observations
  sumry$obs_per_group_min <- x$obs_per_group_min
  sumry$obs_per_group_avg <- x$obs_per_group_avg
  sumry$obs_per_group_max <- x$obs_per_group_max
  sumry$nof_groups = x$nof_groups
  
  sumry$results <- extract(x)
  
  class(sumry) <- "summary.pvarhk"
  sumry
}

#' S3 Print Method for summary.pvarhk
#' @param x object
#' @param ... further arguments
#' @method print summary.pvarhk
#' @export
print.summary.pvarhk <- function(x, ...) {
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
#' @method knit_print summary.pvarhk
#' @importFrom knitr knit_print
knit_print.summary.pvarhk <- function(x, ...) {
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
