######################################################################################
######################################################################################
# Predict residuals for the level equation:
######################################################################################
######################################################################################

# Should work for pvarfeols and pvargmm:

# Need fct fixedeffects:
ex1_dahlberg_data <- pvargmm(dependent_vars = c("expenditures", "revenues", "grants"),
                             lags = 1,
                             transformation = "fod",
                             data = Dahlberg,
                             panel_identifier=c("id", "year"),
                             steps = c("twostep"),
                             system_instruments = FALSE,
                             max_instr_dependent_vars = 99,
                             max_instr_predet_vars = 99,
                             min_instr_dependent_vars = 2L,
                             min_instr_predet_vars = 1L,
                             collapse = FALSE
)

ex2_dahlberg_data <- pvargmm(dependent_vars = c("expenditures", "revenues", "grants"),
                             lags = 2,
                             transformation = "fod",
                             data = Dahlberg,
                             panel_identifier=c("id", "year"),
                             steps = c("twostep"),
                             system_instruments = FALSE,
                             max_instr_dependent_vars = 99,
                             max_instr_predet_vars = 99,
                             min_instr_dependent_vars = 2L,
                             min_instr_predet_vars = 1L,
                             collapse = FALSE
)

ex3_dahlberg_data <- pvargmm(dependent_vars = c("expenditures", "revenues", "grants"),
                             lags = 2,
                             transformation = "fod",
                             data = Dahlberg,
                             panel_identifier=c("id", "year"),
                             steps = c("twostep"),
                             system_instruments = TRUE,
                             max_instr_dependent_vars = 99,
                             max_instr_predet_vars = 99,
                             min_instr_dependent_vars = 2L,
                             min_instr_predet_vars = 1L,
                             collapse = FALSE
)

ex1_dahlberg_data_residuals <- residuals_level.pvargmm(model = ex1_dahlberg_data)
ex2_dahlberg_data_residuals <- residuals_level.pvargmm(model = ex2_dahlberg_data)
ex3_dahlberg_data_residuals <- residuals_level.pvargmm(model = ex3_dahlberg_data)


# model <- ex3_dahlberg_data


#' 
#' 
#' Calculate fixed effect
#' 
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
