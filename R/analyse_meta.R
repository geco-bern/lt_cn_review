analyse_meta <- function(df, nam_target){
  
  ## df: data frame of the data, containing a column 'var' which holds the names of the target variable
  ## nam_target: a character string specifying the name of the target variable
  
  if (nrow(filter(df, var == !!nam_target)) < 3){
    
    warning(paste("Too little data for meta analysis for", nam_target))
    
    modl <- NA
    df_box <- tibble()
    
  } else {
   
    # adjusted step size (see http://www.metafor-project.org/doku.php/tips:convergence_problems_rma)
    modl <- try(metafor::rma.mv( 
      logr, 
      logr_var,
      method = "REML", 
      random = ~ 1 | exp, 
      slab = exp,
      control = list(stepadj = 0.3), 
      data = df %>% 
        filter(var == !!nam_target)
    ))
    
    if (class(modl)[[1]] == "try-error"){
      
      modl <- NA
      df_box <- tibble()
      
    } else {
      
      # get estimate and confidence intervals of scaled variable 
      modl_scaled <- predict( modl, transf = exp )
      
      df_box <- tibble(
        var = nam_target,
        middle = modl$b[1,1],
        ymin   = modl$ci.lb,
        ymax   = modl$ci.ub,
        middle_scaled = modl_scaled$pred,
        ymin_scaled   = modl_scaled$ci.lb,
        ymax_scaled   = modl_scaled$ci.ub
      )
      
    } 
  }
  
  return(list(df_box = df_box, modl = modl))  
}