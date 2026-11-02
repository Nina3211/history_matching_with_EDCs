# the second set of EDCs requires us to run DALEC Crop to check the simulator output
# REVISION: only change is the bug fix below (EDC_SOM1 was never applied; fail_EDC_SOM2
# appeared twice) and a counter for DALEC-Crop evaluations (used if hm_runner.R is loaded).
check_edcs2 = function(input_pars){
  if (exists("count_dalec")) count_dalec(nrow(input_pars))
  # this function assumes that we have the drivers and CARDAMOM already set up
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
  
  # We have three EDCs to check, plus whether all fluxes are positive
  fail_DS = C_cycle$DS_max < 1.9
  fail_EDC_SOM1 = abs(C_cycle$EDC_SOM1) > log(10)
  fail_EDC_SOM2 = abs(C_cycle$EDC_SOM2) > 0.2
  cat("fail DS: ", sum(fail_DS), "\n")
  cat("fail fail_EDC_SOM1: ", sum(fail_EDC_SOM1), "\n")
  cat("fail fail_EDC_SOM2: ", sum(fail_EDC_SOM2), "\n")
  
  fail_EDCs = fail_DS | fail_EDC_SOM1 | fail_EDC_SOM2   # BUG FIX: was fail_EDC_SOM2 twice
  pass_EDCs = !fail_EDCs
  
  # ensure minimum pool values are >= 0 and /= NaN (which is automatic in checking >=0)and not infinity 
  pass_pool1 = (apply(C_cycle$labile_gCm2, 1, min) >= 0) & (!is.infinite(apply(C_cycle$labile_gCm2, 1, max)))
  pass_pool2 = (apply(C_cycle$foliage_gCm2, 1, min) >= 0) & (!is.infinite(apply(C_cycle$foliage_gCm2, 1, max)))
  pass_pool3 = (apply(C_cycle$roots_gCm2, 1, min) >= 0) & (!is.infinite(apply(C_cycle$roots_gCm2, 1, max)))
  pass_pool4 = (apply(C_cycle$wood_gCm2, 1, min) >= 0) & (!is.infinite(apply(C_cycle$wood_gCm2, 1, max)))
  pass_pool5 = (apply(C_cycle$litter_gCm2, 1, min) >= 0) & (!is.infinite(apply(C_cycle$litter_gCm2, 1, max)))
  pass_pool6 = (apply(C_cycle$som_gCm2, 1, min) >= 0) & (!is.infinite(apply(C_cycle$som_gCm2, 1, max)))
  pass_pool7 = (apply(C_cycle$autotrophic_gCm2, 1, min) >= 0) & (!is.infinite(apply(C_cycle$autotrophic_gCm2, 1, max)))
  pass_pool9 = (apply(C_cycle$StorageOrgan_gCm2, 1, min) >= 0) & (!is.infinite(apply(C_cycle$StorageOrgan_gCm2, 1, max)))
  pass_pool10 = (apply(C_cycle$DeadFoliage_gCm2, 1, min) >= 0) & (!is.infinite(apply(C_cycle$DeadFoliage_gCm2, 1, max)))
  
  pass_pools = pass_pool1 & pass_pool2 & pass_pool3 & pass_pool4 & pass_pool5 & 
    pass_pool6 & pass_pool7 & pass_pool9 & pass_pool10
  cat("pass_pools", sum(pass_pools), "\n")
  
  # ensure fluxes are not infinity or nan
  # check if all fluxes are not infinite
  flux_names = c("gpp_gCm2day", "rauto_gCm2day", "rhet_litter_gCm2day", 
                 "rhet_som_gCm2day", "harvest_gCm2day", "extracted_residue_gCm2day", "alloc_foliage_gCm2day",
                 "alloc_labile_gCm2day", "alloc_roots_gCm2day", "alloc_wood_gCm2day", "labile_to_foliage_gCm2day", 
                 "alloc_autotrophic_gCm2day", "alloc_StorageOrgan_gCm2day", "foliage_to_litter_gCm2day", 
                 "roots_to_litter_gCm2day", "wood_to_litter_gCm2day", "litter_to_som_gCm2day", 
                 "rauto_maintenance_gCm2day", "rauto_labile_to_foliage_gCm2day", "rauto_npp_to_labile_gCm2day", 
                 "rauto_foliage_to_litter_gCm2day", "rauto_wood_to_litter_gCm2day", "HARVESTextracted_foliage_gCm2day", 
                 "HARVESTextracted_wood_gCm2day", "HARVESTextracted_DeadFoliage_gCm2day", "HARVESTextracted_labile_gCm2day", 
                 "HARVESTlitter_foliage_gCm2day", "HARVESTlitter_wood_gCm2day", "HARVESTlitter_DeadFoliage_gCm2day", 
                 "HARVESTlitter_autotrophic_gCm2day", "HARVESTlitter_labile_gCm2day", "PLOUGHlitter_roots_gCm2day", 
                 "ET_kgH2Om2day", "Etrans_kgH2Om2day", "Esoil_kgH2Om2day", "Ewetcanopy_kgH2Om2day", 
                 "runoff_kgH2Om2day", "underflow_kgH2Om2day", "SurfWater_kgH2Om2")
  pass_fluxes_df = data.frame(matrix(nrow = length(pass_pools), ncol = length(flux_names))) 
  colnames(pass_fluxes_df) = flux_names
  for (i in seq_along(flux_names)) {
    current_flux = flux_names[i]
    pass_fluxes_df[,current_flux] = !is.infinite(apply(C_cycle[[current_flux]], 1, max))
  }
  pass_fluxes = apply(pass_fluxes_df, 1, all)
  
  cat("pass_fluxs", sum(pass_fluxes), "\n")
  
  # now combine all passes
  pass_all = pass_EDCs & pass_pools &  pass_fluxes
  cat("So many points pass the 2nd EDCs: ", sum(pass_all), "\n")
  return(input_pars[pass_all,])
}

assess_edcs2 = function(input_pars, parnames = vars_DALEC.A3.C3.H2.M1, 
                        fixed_pars = list("sow_day" = 274, "harvest_length" = 342)){

  if(ncol(input_pars) < length(parnames)){
    missing_pars = parnames[!parnames %in% colnames(input_pars)]
    # find out which parameters are missing
    cat("the following pars are missing", missing_pars, "\n")
    inputs_par_train = list()
    for(i in missing_pars){
      inputs_par_train[[i]] = rep(fixed_pars[[i]], fixed_parstimes = nrow(input_pars))
      
    }
    fixed_df = do.call("cbind", inputs_par_train)
    input_pars = cbind(input_pars, fixed_df)[parnames] # sort parameters in right order
  }
  
  if(nrow(input_pars) > 20000){
    # Calculate the total number of batches
    batch_size <- 20000
    print(nrow(input_pars))
    
    num_batches <- ceiling(nrow(input_pars) / batch_size)
    batches <- list()
    for (i in 1:num_batches) {
      print("splitting LHD into batches")
      cat("Batch ", i, 1, "\n")
      start_row <- ((i - 1) * batch_size) + 1
      end_row <- min(i * batch_size, nrow(input_pars))  # Ensure not to exceed the number of rows
      batches[[i]] <- check_edcs2(input_pars[start_row:end_row, ])
    }
    input_pars = do.call('rbind', batches)
  }
  else{
    input_pars = check_edcs2(input_pars)
  }
  return(input_pars)
}

