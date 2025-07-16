initial_LHD = function(n_points_requested, use_EDCs = TRUE, 
                       other_fixed_vars, pars_with_priors, fixed_vars,
                       ranges_df_with_priors){
  
  points.factor <- 4 # how many more points we are proposing
  
  nu_of_pars_all = nrow(ranges_df_with_priors)
  nu_of_pars_LHS = nu_of_pars_all - length(fixed_vars) - length(other_fixed_vars)
  # which parameters for LHD?
  parameters_active = ranges_df_with_priors$par_name[!ranges_df_with_priors$par_name %in% c(fixed_vars, other_fixed_vars)]
  print(parameters_active)
  # do we want to set any priors?
  if(length(pars_with_priors) > 0){
    ranges_df_with_priors[ranges_df_with_priors$par_name %in% pars_with_priors, c("par_min", "par_max")] = ranges_df_with_priors[ranges_df_with_priors$par_name %in% pars_with_priors, c("potential_min", "potential_max")]
  }
  
  # drop from ranges df all constant cols
  ranges_df_with_priors_active = ranges_df_with_priors[ranges_df_with_priors$par_name %in% parameters_active, c("par_min", "par_max", "par_name")]
  # produce ranges list for emulator with final ranges
  wave0_ranges = split(ranges_df_with_priors_active[, c("par_min", "par_max")], seq(nrow(ranges_df_with_priors_active)))
  names(wave0_ranges) = parameters_active

  
  final_LHD = data.frame()
  
  while(nrow(final_LHD) < n_points_requested){
    
    # first generate LHD scaled between [-1,1]
    LHD = setNames(
      data.frame(2* (lhs::randomLHS(n_points_requested * points.factor, length(parameters_active)) - 0.5)),
      names(parameters_active)
    )
    
    # Scale each column to the original parameter range
    scaled_lhs <- as.data.frame(
      mapply(
        function(values, par_min, par_max) par_min + (par_max - par_min) * ((values + 1) / 2),
        LHD,
        ranges_df_with_priors_active$par_min,
        ranges_df_with_priors_active$par_max
      )
    )
    
    # Ensure column names are preserved
    names(scaled_lhs) <- ranges_df_with_priors_active$par_name

    if(use_EDCs){
      scaled_lhs = assess_EDC1(scaled_lhs, fixed_pars = prior_vals_list) # this adds the fixed variables
      if(nrow(scaled_lhs)>0){
        print("assessing EDCs 2")
        scaled_lhs = assess_edcs2(scaled_lhs)}
    }
    
    final_LHD = rbind(final_LHD, scaled_lhs)
    
  }
  
  # choose final subset using maximin
  print("choosing final points using maximin sampling.")
  final_LHD = hmer:::maximin_sample(final_LHD, n = n_points_requested, nms = names(scaled_lhs))

  return(list("initial_LHD" = final_LHD, "wave0_ranges" = wave0_ranges))
}