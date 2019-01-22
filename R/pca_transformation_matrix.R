# PCA Transformation Matrix
#
#@param corrected_instruments_F 
#@param nof_categories Number of categories
#@param nof_predet_vars Number of predetermined variables
#@param nof_exog_vars Number of exogeneous variables
#@param system_instruments System GMM estimator
#@param system_constant Constant only available with the System GMM estimator in each equation
#@param position_exo TODO
#@param pca_eigenvalue Cut-off eigenvalue for pca analysis
#
#@description PCA transformation of instruments matrix
#@references Roodman

pca_transformation <- function(corrected_instruments_F,
                               nof_categories,
                               nof_predet_vars,
                               nof_exog_vars,
                               system_instruments,
                               system_constant,
                               position_exo,
                               pca_eigenvalue){
  `%ni%`<- Negate(`%in%`)
  original_instruments <- corrected_instruments_F
  # Remove rows the exogenous instruments:
  # Be careful if there are system instruments (the row of 1 in the last.)
  if (nof_exog_vars > 0){

    # Make sure that "FALSE" is the default value:
    if (system_instruments == FALSE){

      corrected_instruments_F <-
        mapply(function(i) head(original_instruments[[i]], -nof_exog_vars), 1:nof_categories, SIMPLIFY = FALSE)

    }

    if (system_instruments == TRUE){

      # remove the row of 1.
      if (system_constant == TRUE){
        corrected_instruments_F <-
          mapply(function(i) head(original_instruments[[i]], -1),1:nof_categories, SIMPLIFY = FALSE)
      }

      # now remove the rows of the FD-Exogenous variables:
      corrected_instruments_F <-
        mapply(function(i) original_instruments[[i]][-(position_exo:(position_exo+nof_exog_vars-1)),],1:nof_categories, SIMPLIFY = FALSE)

    }


  }

  corrected_instruments_F_T <- list()

  for (i0 in 1:length(corrected_instruments_F)){

    corrected_instruments_F_T[[i0]] <- t(corrected_instruments_F[[i0]])

    corrected_instruments_F_T[[i0]] <- corrected_instruments_F_T[[i0]][rowSums(corrected_instruments_F_T[[i0]] != 0) > 0,]

  }

  corrected_instruments_F_T_pca <- corrected_instruments_F_T[[1]]

  if (length(corrected_instruments_F) > 1){
    for (i0 in 2:length(corrected_instruments_F)){

      corrected_instruments_F_T_pca <- rbind(corrected_instruments_F_T_pca, corrected_instruments_F_T[[i0]])

    }
  }
  # save zero cols:
  zero_rows <- unique(which(colSums(corrected_instruments_F_T_pca) == 0))

  t_corrected_instruments_F_pca <- corrected_instruments_F_T_pca[,colSums(corrected_instruments_F_T_pca != 0) > 0]

  # Remove the row with the constant for the pca analysis if sytsem = TRUE:s
  if (system_instruments == TRUE){

    zero_rows <- c(zero_rows, ncol(corrected_instruments_F_T_pca))

    t_corrected_instruments_F_pca <- t_corrected_instruments_F_pca[,-ncol(t_corrected_instruments_F_pca)]

  }

  # Use the new t_corrected_instruments_F_pca as input:
  ir.pca <- prcomp(t_corrected_instruments_F_pca,
                   center = TRUE,
                   scale. = TRUE,
                   tol = 1e-15)

  ir.pca.components <- ir.pca$sdev*ir.pca$sdev
  #ir.pca.components

  nof_pca_rotation_cols <- length(ir.pca.components[ir.pca.components > pca_eigenvalue ])


  # Take care of rows that carry no information in the pca:
  # Add them for the matrix multiplication:

  if (dim(t_corrected_instruments_F_pca)[2] < dim(corrected_instruments_F_T_pca)[2]){

    zero_cols <- dim(corrected_instruments_F_T_pca)[2] - dim(t_corrected_instruments_F_pca)[2]

    added.test.rotation <- matrix(0, nrow = nrow(ir.pca$rotation)+zero_cols, ncol = ncol(ir.pca$rotation))

    # Add one column because the matrices do not fit otherwise.
    not_zero_rows <- c(c(1:dim(corrected_instruments_F_T_pca)[2]) %ni% zero_rows)

    # Add row to added.test.rotation at certain positions to solve dimension problem:
    added.test.rotation[not_zero_rows,] <- ir.pca$rotation
    #added.test.rotation[1:nrow(ir.pca$rotation), 1:ncol(ir.pca$rotation)] <- ir.pca$rotation

  }

  if (dim(t_corrected_instruments_F_pca)[2] == dim(corrected_instruments_F_T_pca)[2]){

    added.test.rotation <- ir.pca$rotation

  }


  z_test <- list()

  for (i0 in 1:length(corrected_instruments_F)){

    z_test[[i0]]   <- t(as.matrix(corrected_instruments_F[[i0]])) %*% added.test.rotation[,1:nof_pca_rotation_cols]

  }


  Q <- mapply(function(i) t(z_test[[i]]), 1:nof_categories, SIMPLIFY = FALSE)

  # Add instrument for exogenous variables if included:
  if (nof_exog_vars > 0){

    instruments_exogenous <- mapply(function(i) original_instruments[[i]][(position_exo:(position_exo+nof_exog_vars-1)),],
                                    1:nof_categories, SIMPLIFY = FALSE)

    # Add these variables in the correct position.
    Q <-  mapply(function(i) rbind(t(z_test[[i]]), instruments_exogenous[[i]]), 1:nof_categories, SIMPLIFY = FALSE)

  }

  # Add the removed colum of the 1s:
  if (system_instruments == TRUE){

    add_system_instruments_constant <- mapply(function(i) tail(original_instruments[[i]], 1),
                                              1:nof_categories, SIMPLIFY = FALSE)

    #Q <-  mapply(function(i) rbind(t(z_test[[i]]), add_system_instruments_constant[[i]]), 1:nof_categories, SIMPLIFY = FALSE)
    Q <-  mapply(function(i) rbind(Q[[i]], add_system_instruments_constant[[i]]), 1:nof_categories, SIMPLIFY = FALSE)

  }

  return(Q)

}