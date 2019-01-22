#################################################################################
# Function to transform Instruments of the endogenous variables in the PVAR code.
#################################################################################


Endo_Lmin_Lmax_transition_matrix<- function(lags, min_instr_dependent_vars, max_instr_dependent_vars, nofperiods, nof_dependent_vars,
                                            collapse){


  max_instr_dependent_vars_original <- max_instr_dependent_vars
  if (min_instr_dependent_vars > 2){
    max_instr_dependent_vars <- max_instr_dependent_vars - min_instr_dependent_vars + 2
  }

  if (min_instr_dependent_vars == 2){
    max_instr_dependent_vars <- max_instr_dependent_vars - min_instr_dependent_vars + 1
  }

  # Build now the identity matrices.
  identity_matrices <- list()

  for (i in 1:(nofperiods-2)){

    identity_matrices[[i]] <- diag(i)

    # L_max cut-offs.
    if (i > (max_instr_dependent_vars)){

      # L_max cuf-off columns.
      identity_matrices[[i]] <-  identity_matrices[[i]][, -((max_instr_dependent_vars+1):ncol(identity_matrices[[i]]))]

      if (is.null(dim(identity_matrices[[i]]))){

        identity_matrices[[i]]  <- matrix(identity_matrices[[i]], nrow = 1, ncol = length(identity_matrices[[i]]))


      }
      # Fix problems with max_instr_dependent_vars_original == min_instr_dependent_vars
      if (max_instr_dependent_vars_original == min_instr_dependent_vars){
        identity_matrices[[i]] <- matrix(identity_matrices[[i]], byrow = FALSE)
      }
      # L_max set rows to 0.
      if(nrow(identity_matrices[[i]]) > 1) {
        identity_matrices[[i]][(max_instr_dependent_vars+1):i,] <- 0
      }

    }

    # L_min operations.
    if (min_instr_dependent_vars > 2){

      # Cut of columns that at the beginning of each matrix. 2 would be the starting point of instruments
      # if there were no restrictions on min_instr_dependent_vars
      if (i == 1){

        identity_matrices[[i]] <- 0

      }

      if (i > 1){
        identity_matrices[[i]] <- as.matrix(identity_matrices[[i]][,-(1:(min_instr_dependent_vars-2))])
      }

    }

  }

  # Adjust the identity matices if there are more lags: Cut out the first few identity matices.
  if (lags > 1){

    identity_matrices <- identity_matrices[-(1:(lags-1))]

  }


  if (collapse == TRUE){
    for (i in 1:length(identity_matrices)){

      if (!is.matrix(identity_matrices[[i]])){identity_matrices[[i]] <- matrix(identity_matrices[[i]])}

      identity_matrices[[i]] <- cbind(identity_matrices[[i]], matrix(0, nrow = nrow(identity_matrices[[i]]),
                                                                     ncol = ncol(identity_matrices[[length(identity_matrices)]]) -
                                                                       ncol(identity_matrices[[i]])  ))

    }

    results <- do.call(rbind,identity_matrices)

  }

  if (collapse == FALSE){

    results <-  Matrix::Matrix(Matrix::bdiag(identity_matrices), sparse = TRUE)

  }


  if (nof_dependent_vars > 1){

    results <- results %x% diag(nof_dependent_vars)

  }

  return(results)
}

