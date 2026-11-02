assess_EDC1 = function(pars, fixed_pars = list("sow_day" = 274, "harvest_length" = 342), 
                       parnames = vars_DALEC.A3.C3.H2.M1){
  # pars = pars[,vars_DALEC.A3.C3.H2.M1] 
  # check if we are given a full dataset including fixed variables
  if(ncol(pars) < length(parnames)){
    missing_pars = parnames[!parnames %in% colnames(pars)]
    # find out which parameters are missing
    cat("the following pars are missing", missing_pars, "\n")
    inputs_par_train = list()
    for(i in missing_pars){
      inputs_par_train[[i]] = rep(fixed_pars[[i]], fixed_parstimes = nrow(pars))
      
    }
    fixed_df = do.call("cbind", inputs_par_train)
    pars = cbind(pars, fixed_df)[parnames] # sort parameters in right order
  }

  # Turnover of litter faster than turnover of som
  fail1 = (pars[,"min_som"] > pars[,"min_litter"])

  # decomposition of litter to SOM greater than SOM to air
  fail2 = (pars[,"min_som"] > pars[,"decomp_rate"])
  
  # pre_DR should be greater than post_DR, this is consistent across currently
  # available SPA crop parameter files
  fail3 = (pars[,"DR_post"] > pars[,"DR_pre"])
  
  # for development: Tmin should be < topt and topt should be < tmax
  fail4 = ((pars[,"tmin"] > pars[,"topt"]) | (pars[,"topt"] > pars[,"tmax"])) 
  
  
  # for development: the difference between each Tmin,Topt,Tmax > 1.
  fail5 = ((abs(pars[,"tmin"]-pars[,"topt"]) < 1) | (abs(pars[,"topt"]-pars[,"tmax"]) < 1)| (abs(pars[,"tmin"]-pars[,"tmax"]) < 1)) 
  
  
  # for vernalisation: Tmin < Topt < Tmax
  fail6 = ((pars[,"tmin_v"] > pars[,"topt_v"]) | (pars[,"topt_v"] > pars[,"tmax_v"])) 
  
  
  # for vernalisation: the difference between each Tmin, Topt, Tmax
  fail7 = ( (abs(pars[,"tmin_v"]-pars[,"topt_v"]) < 1) | (abs(pars[,"tmin_v"]-pars[,"tmax_v"]) < 1 ) ) 
  
  
  # development temperature value should be larger corresponding vernalisation
  fail8 = ((pars[,"tmin_v"] > pars[,"tmin"]) | (pars[,"topt_v"] > pars[,"topt"]) | (pars[,"tmax_v"] > pars[,"tmax"]))
  
  
  all_fails = cbind(fail1, fail2,  fail3, fail4,  fail5,  fail6, fail7,  fail8)
  row_has_true = apply(all_fails, 1, any)
  cat("EDC1 reject ", sum(row_has_true), "out of ", length(fail1), " points \n")
  all_trues = !row_has_true
  
  return(pars[all_trues,])
}