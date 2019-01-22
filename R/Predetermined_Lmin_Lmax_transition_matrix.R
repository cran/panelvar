####################################################################################
# Function to transform Instruments of the predetermined variables in the PVAR code.
####################################################################################



Predetermined_Lmin_Lmax_transition_matrix <- function(lags, min_instr_predet_vars, max_instr_predet_vars, nofperiods, nof_predet_vars, nof_dependent_vars, collapse){


  identity_matrices <- list()

  for (i in 2:(nofperiods-1)){

    identity_matrices[[i]] <- diag(i)

    # L_max cut-offs.
    if (i > (max_instr_predet_vars)){

      # L_max cuf-off columns.
      identity_matrices[[i]] <- identity_matrices[[i]][, -((max_instr_predet_vars+1):ncol(identity_matrices[[i]]))]

      # L_max set rows to 0.
      if (is.null(dim(identity_matrices[[i]]))){

        identity_matrices[[i]]  <- matrix(identity_matrices[[i]], nrow = 1, ncol = length(identity_matrices[[i]]))


      }

      if (min_instr_predet_vars == max_instr_predet_vars){
        identity_matrices[[i]] <- matrix(identity_matrices[[i]], byrow = FALSE)
      }

      if(nrow(identity_matrices[[i]]) > 1) {

        identity_matrices[[i]][(max_instr_predet_vars+1):i,] <- 0

      }

    }


    if (min_instr_predet_vars > 1){

      identity_matrices[[i]] <- as.matrix(identity_matrices[[i]][,-(1:(min_instr_predet_vars-1))])

    }

  }

  # Adjust the identity matices if there are more lags: Cut out the first few identity matrices.
  if (lags > 1){

    identity_matrices <- identity_matrices[-(1:(lags-1))]

  }


  identity_matrices<- identity_matrices[-1]


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
    results <- Matrix::Matrix(Matrix::bdiag(identity_matrices), sparse = TRUE)
  }



  if (nof_predet_vars > 1){

    results <- results %x% diag(nof_predet_vars)

  }

  return(results)


}