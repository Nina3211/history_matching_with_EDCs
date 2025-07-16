generate_input_output_DALEC_Crop = function(input_pars, fixed_pars = list("sow_day" = 274, "harvest_length" = 342), 
                          parname = vars_DALEC.A3.C3.H2.M1, obs_dates, obs_dates_names,
                          year = 2017, is_final_wave = FALSE){
  # check if we are given a full dataset including fixed variables
  if(ncol(input_pars) < length(vars_DALEC.A3.C3.H2.M1)){
    missing_pars = vars_DALEC.A3.C3.H2.M1[!vars_DALEC.A3.C3.H2.M1 %in% colnames(input_pars)]
    # find out which parameters are missing
    cat("the following pars are missing", missing_pars, "\n")
    inputs_par_train = list()
    for(i in missing_pars){
      inputs_par_train[[i]] = rep(fixed_pars[[i]], fixed_parstimes = nrow(input_pars))

    }
    fixed_df = do.call("cbind", inputs_par_train)
    input_pars = cbind(input_pars, fixed_df)[vars_DALEC.A3.C3.H2.M1] # sort parameters in right order
  }
  
  parameters_in = t(input_pars)
  
  ### Code to run DALEC Crop 
  ## Manipulate the parameters if you want 
  n = 1 # one site
  # Fit into format of (npar,ensembmeterle_nos,chain)
  parameters_in = array(as.matrix(parameters_in), dim = c(PROJECT$model$nopars[n], dim(parameters_in)[2], 1))
  
  soil_info = c(drivers$top_sand,drivers$bot_sand,drivers$top_clay,drivers$bot_clay)
  C_cycle = simulate_all(n,new_PROJECT,PROJECT$model$name,drivers$met,parameters_in,
                         drivers$lat,PROJECT$ctessel_pft[n],PROJECT$parameter_type,
                         PROJECT$exepath,soil_info)
  
  # only plotting "new" model runs, not the previously calibrated ones
  par(mfrow=c(1,3))
  plot(apply(C_cycle$gpp_gCm2day,2,median), type="l", lwd=3, ylim=c(0,15), main="GPP")
  plot(apply(C_cycle$lai_m2m2,2,median), type="l", lwd=3, ylim=c(0,10), main="LAI")
  lines(apply(C_cycle$lai_m2m2,2,mean), col="blue", lwd=3)
  lines(apply(C_cycle$lai_m2m2,2,mean) + 3*apply(C_cycle$lai_m2m2,2,sd), col="grey", lwd=3)
  lines(apply(C_cycle$lai_m2m2,2,mean) - 3*apply(C_cycle$lai_m2m2,2,sd), col="grey", lwd=3)
  plot(apply(C_cycle$StorageOrgan_gCm2,2,median), type="l", lwd=3, ylim=c(0,350), main="StorageOrgan")
  lines(apply(C_cycle$StorageOrgan_gCm2,2,mean), col="blue", lwd=3)
  lines(apply(C_cycle$StorageOrgan_gCm2,2,mean) + 3*apply(C_cycle$StorageOrgan_gCm2,2,sd), col="grey", lwd=3)
  lines(apply(C_cycle$StorageOrgan_gCm2,2,mean) -3*apply(C_cycle$StorageOrgan_gCm2,2,sd), col="grey", lwd=3)
  
  # we are interested in the following outputs: LAI, harvest, yield_GPP_ratio, ra_GPP_ratio, EDC_SOM1, EDC_SOM2, DS_max
  # LAI
  LAI_out = as.data.frame(t(as.matrix(C_cycle$foliage_gCm2)))
  print(dim(LAI_out))
  LAI_out = as.data.frame(lapply(LAI_out, function(col) {
    zoo::rollmean(col, k = 2, align = "left", fill = NA)
  }))
  LAI_out = sweep(LAI_out, 2, input_pars$LMA, FUN = "/")
  LAI_out$date = seq(as.Date("2016-01-01"), by = "day", length.out = dim(LAI_out)[1])
  
  if(year == 2017){
    current_LAI = LAI_out %>% dplyr::filter(between(date, as.Date("2016-09-01"), as.Date("2017-09-30")))
  } else{
    current_LAI = LAI_out %>% dplyr::filter(between(date, as.Date("2017-09-01"), as.Date("2018-09-30")))
  }

  # now subset according to the days that we have observations
  LAI_subset = as.data.frame(t(as.matrix(subset(LAI_out, date %in% obs_dates) %>% dplyr::select(-date))))
  colnames(LAI_subset) = obs_dates_names
  
    
  # yield at harvest output
  harvest_day = 250 + 366 # 2016 is leap year
  if(year != 2017){harvest_day = 237+365+365}
  harvest_vector = C_cycle$harvest_gCm2day[, harvest_day]
  print(which(colSums(C_cycle$harvest_gCm2day) != 0))
  print(harvest_vector)
  
  
  outputs_df = cbind(LAI_subset, harvest_vector, "EDC_SOM1" = C_cycle$EDC_SOM1, 
                     "EDC_SOM2" = C_cycle$EDC_SOM2, 
                     "ra_GPP_ratio" = C_cycle$ra_GPP_ratio, 
                     "yield_GPP_ratio" =C_cycle$yield_GPP_ratio, 
                     "DS_max" = C_cycle$DS_max)
  training_data = cbind(outputs_df, input_pars)
  
  
  # if final wave, also return the LAI time series, otherwise just return dataframe of input & outputs
  if(is_final_wave == FALSE){
    return(training_data)
  } else{
    # call the C_cycle function once
    return(list("training_data" = training_data, "LAI_curves" = current_LAI))
  }
  
}