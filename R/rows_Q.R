#################################################################################################
#################################################################################################
# Rows of Q based on lagged endogenous variables
#################################################################################################
#################################################################################################
# Basic idea:
# Arithmetric progression that is constant after
# a certain element divide it into two parts
# 1.Part: Normal arithmetric progression (1) value of starting element, (2) number of elements
# 2.Part: Constant time series: (3) value of starting element, (4) number of elements


rows_Q_lagged_endogenous_vars <- function(nof_periods,
                                          nof_dependent_vars,
                                          max_instr_dependent_vars = nof_periods,
                                          lags,
                                          min_instr_dependent_vars = 2){

  # Set up the four parts of the double arithmetic progression and put them together at the end.
  #(1) ###########################################################################
  if(max_instr_dependent_vars+1 == min_instr_dependent_vars){
    a_1 <- 1
  }

  if (max_instr_dependent_vars+1 != min_instr_dependent_vars){

    if (lags > min_instr_dependent_vars){

      a_1 <- min(min(max_instr_dependent_vars + 2, nof_periods-2)- min_instr_dependent_vars, lags)
    }

    if (lags <= min_instr_dependent_vars){
      a_1 <- lags - min(min_instr_dependent_vars - 2, lags - 1)
    }
  }
  ################################################################################

  # (2) ##########################################################################
  n_1 <- max( min(max_instr_dependent_vars, nof_periods-2) - max(min_instr_dependent_vars - 2, lags - 1) ,1)
  ################################################################################

  # (3) ##########################################################################
  a_2 <- a_1 + (n_1-1)*1
  ################################################################################

  # (4) ##########################################################################
  n_2 <- nof_periods - n_1 -(min(lags, 2) + max(min_instr_dependent_vars-lags, lags-1))
  ################################################################################

  # Combine a_1, n_1, a_2, n_2:
  M_lagged_endo <- nof_dependent_vars*((n_1/2 * (2*a_1 + (n_1-1))) + a_2*n_2)



  return(M_lagged_endo)
}

#################################################################################################
#################################################################################################
# Rows of Q based on lagged endogenous variables
#################################################################################################
#################################################################################################
# Basic idea:
# Arithmetric progression that is constant after
# a certain element divide it into two parts
# 1.Part: Normal arithmetric progression (1) value of starting element, (2) number of elements
# 2.Part: Constant time series: (3) value of starting element, (4) number of elements


rows_Q_predetermined_vars <- function(nof_periods,
                                      nof_predet_vars,
                                      max_instr_predet_vars = nof_periods,
                                      lags,
                                      min_instr_predet_vars = 1){



  if(max_instr_predet_vars+1 == min_instr_predet_vars){
    b_1 <- 1
  }

  if (max_instr_predet_vars+1 != min_instr_predet_vars){

    if (lags > min_instr_predet_vars){
      b_1 <- max(min(max_instr_predet_vars + 2- min_instr_predet_vars,lags+1), lags-min_instr_predet_vars)
    }

    if (lags <= min_instr_predet_vars){
      b_1 <- lags + 1 - max(min_instr_predet_vars - 1, lags - 1)
    }

  }
  ################################################################################

  # (2) ##########################################################################

  o_1_old <- max( min(nof_periods-lags, max_instr_predet_vars) -
                (max(min_instr_predet_vars - 1, lags - 1)), 
              max(min(nof_periods-lags,max_instr_predet_vars)  + 1 - b_1), 1)
  
  o_1 <- max(min(nof_periods-2, max_instr_predet_vars)-
                       (max(min_instr_predet_vars - 1, lags - 1)), 
                     max(min(nof_periods-lags,max_instr_predet_vars) + 1 - b_1),1)

  if (lags > 2){
    o_1 <- nof_periods-lags-1
    
  }


  ################################################################################

  # (3) ##########################################################################
  b_2 <- b_1 + (o_1 - 1)*1
  ################################################################################

  # (4) ##########################################################################
  o_2 <- nof_periods - o_1 -(min(lags, 1) + max(min_instr_predet_vars - lags, lags))
  ################################################################################

  # Combine b_1, o_1, b_2, o_2:
  M_predetermined <- nof_predet_vars*((o_1/2 * (2*b_1 + (o_1-1))) + b_2*o_2)




  return(M_predetermined)

}
