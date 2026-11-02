parmin = vector(length = 38)
parmax = vector(length = 38)
parname = vector(length = 38)

# Decomposition rate (frac/day)
parname[1] = "decomp_rate"
parmin[1] = 0.0001141 # 24   years at 0oC
parmax[1] = 0.02      # 0.13 years at 0oC

# Fraction of GPP to autotrophic pool
# NOTE: this does not include Ra due to C transfers 
# after initial fixation
parname[2] = "frac_GPP_resp_auto"
parmin[2] = 0.2
parmax[2] = 0.6

# max development rate (day-1) DS (0->1)
parname[3] = "DR_pre"
parmin[3] = 0.01 
parmax[3] = 0.05

# max development rate (day-1) DS (1->2)
parname[4] = "DR_post"
parmin[4] = 0.01 
parmax[4] = 0.05

# turnover rate foliage (frac/day)
parname[5] = "TOR_fol" 
parmin[5] = 0.002737851 # 365.25 days
parmax[5] = 0.03571429  # 28 days 

# TOR stem* - 1% loss per year value (day-1)
parname[6] = "TOR_stem"
parmin[6] = 0.002737851 # 365.25 days
parmax[6] = 0.03571429  # 28 days
# maximum rate of foliar turnover (day-1) due to self-shading
parname[7] = "RDRSHMAX"
parmin[7] = 0.02380952  # 42 days
parmax[7] = 0.07142857  # 17 days

# effective vernalisation days when plants are 50 % vernalised
parname[8] = "VDh"
parmin[8] = 12 
parmax[8] = 32

# Turnover of litter (fraction; temperature adjusted)
parname[9] = "min_litter"
parmin[9] = 0.0001141 # 24   years at 0oC
parmax[9] = 0.02      # 0.13 years at 0oC

# Turnover of som to Rhet (fraction; temperature adjusted)
parname[10] = "min_som"
parmin[10] = 1.368925e-06   # 2000 years at 0oC
parmax[10] = 9.126169e-05   #   30 years at 0oC #0.0001368926 #   20 years at 0oC

# Photosynthetic nitrogen use efficiency (gC/gN/m2/day)
parname[11] = "NUE"
parmin[11] = 1
parmax[11] = 40

# sow day
parname[12] = "sow_day"
parmin[12] = 365.25 #100 
parmax[12] = 2*365.25

# phenological heat units required for emergence
parname[13] = "PHUem"
parmin[13] = 100 
parmax[13] = 150

# Growing season length, i.e. harvest date is sowing + this value
parname[14] = "harvest_length"
parmin[14] = 90
parmax[14] = 365

# Intercept value for canopy N (gN/m2) dilution model
parname[15] = "N_interc"
parmin[15] = 2.5
parmax[15] = 7.0 # 5.5

# Coefficient relating crop development stage to N dilution
parname[16] = "N_slope"
parmin[16] = -0.04
parmax[16] = -0.005

# LMA
parname[17] = "LMA"
parmin[17] = 20 
parmax[17] = 100

#
  # NOTE number order not consistent
#
  
  # minimum temperature for development (oC)
parname[26] = "tmin"
parmin[26] = (-1+273.15) 
parmax[26] = ( 8+273.15)  # -1,8
# maximum temperature for development (oC)
parname[27] = "tmax"
parmin[27] = (18+273.15) 
parmax[27] = (36+273.15)   # 18,36
# optimum temperature for development (oC)
parname[28] = "topt"
parmin[28] = ( 8+273.15) 
parmax[28] = (22+273.15)   # 8,22

# minimum temperature for vernalisation (oC)
parname[29] = "tmin_v"
parmin[29] = (-5.3+273.15)
parmax[29] = (-0.3+273.15)   # -5,0
# maximum temperature for vernalisation (oC)
parname[30] = "tmax_v"
parmin[30] = (12.0+273.15) 
parmax[30] = (19.0+273.15)    # 12,19
# optimum temperature for vernalisation (oC)
parname[31] = "topt_v"
parmin[31] = (2.0+273.15) 
parmax[31] = (7.0+273.15)   # 2,7

# critical photoperiod for development (hrs)
parname[32] = "crit_photo"
parmin[32] = 6 
parmax[32] = 12
# photoperiod sensitivity
parname[33] = "photo_sens"
parmin[33] = 0.10 
parmax[33] = 0.35

# turnover rate of labile
parname[34] = "TOR_lab"
parmin[34] = 1e-5
parmax[34] = 0.1
# turnover rate of autotrophic pool
parname[35] = "TOR_resp_auto"
parmin[35] = 0.65 #    1.5 days 
parmax[35] = 0.99  # ~ 1.0 day

# BUCKET - root biomass needed to reach 50 % of max depth
parname[36] = "root_k"
parmin[36] = 10
parmax[36] = 100

# BUCKET - maximum rooting depth
parname[37] = "max_depth"
parmin[37] = 0.5
parmax[37] = 5 #20.0

  # INITIAL VALUES DECLARED HERE
#
  
  # C labile (gC.m-2)
parname[18] = "C_lab"
parmin[18] = 1 
parmax[18] = 10
# C foliar (gC.m-2)
parname[19] = "C_fol"
parmin[19] = 1 
parmax[19] = 5
# C roots (gC.m-2)
parname[20] = "C_roots"
parmin[20] = 1 
parmax[20] = 5
# C_wood (gC.m-2)
parname[21] = "C_stem"
parmin[21] = 1 
parmax[21] = 5
# C litter (gC.m-2)
parname[22] = "C_litter"
parmin[22] = 1 
parmax[22] = 10
# C_som (gC.m-2)
parname[23] = "C_som"
parmin[23] = 100 
parmax[23] = 200000
# C autotrophic pool (gC.m-2)
parname[24] = "C_auto"
parmin[24] = 0.1 
parmax[24] = 5
# C storage organ (gC.m-2)
parname[25] = "C_stor"
parmin[25] = 0.1 
parmax[25] = 1

# Soil water fraction (m3/m3)
parname[38] = "soil_water_frac"
parmin[38] = 0.1 
parmax[38] = 0.9

ranges_df = data.frame(
  "par_name" = parname,
  "par_min" = parmin, 
  "par_max" = parmax
)

vars_DALEC.A3.C3.H2.M1 = ranges_df$par_name

original_ranges = split(ranges_df[, c("par_min", "par_max")], seq(nrow(ranges_df)))
names(original_ranges) = vars_DALEC.A3.C3.H2.M1

# we can also decide to include some priors
PARPRIORS = vector(length = length(vars_DALEC.A3.C3.H2.M1))
PARPRIORUNC = vector(length = length(vars_DALEC.A3.C3.H2.M1))
# derived from the ATEC experiment
#PARPRIORS[1] = 0.00238848     ; PARPRIORUNC[1]  = 0.016370    # Decomposition of litter to som (fraction/day at 0C)
PARPRIORS[2] = 0.43619389     ; PARPRIORUNC[2]  = 0.179778     # Fraction of GPP allocated to autotrophic pool
PARPRIORS[3] = 0.03897653     ; PARPRIORUNC[3]  = 0.026342     # Development rate coefficient DS 0-1
PARPRIORS[4] = 0.01967937     ; PARPRIORUNC[4]  = 0.028072     # Development rate coefficient DS 1-2
PARPRIORS[5] = 0.00896606     ; PARPRIORUNC[5]  = 0.027046     # Potential turnover rate of foliage due to DS (fraction/day)
PARPRIORS[6] = 0.01014908     ; PARPRIORUNC[6]  = 0.027931     # Potential turnover rate of stem due to DS (fraction/day)
PARPRIORS[7] = 0.04080012     ; PARPRIORUNC[7]  = 0.042867     # Potential turnover rate of foliage due to self-shading (fraction/day)
PARPRIORS[8] = 19.3456794     ; PARPRIORUNC[8]  = 17.70075     # No. of vernalisation days for plants to be 50 % vernalised
#PARPRIORS[9] = 0.00106709     ; PARPRIORUNC[9]  = 0.01343840   # Mineralisation rate of litter (fraction/day at 0C)
#PARPRIORS[10] = 7.7304288e-06 ; PARPRIORUNC[10] = 5.771160e-05 # Mineralisation rate of som (fraction/day at 0C)
PARPRIORS[11] = 21.1491       ; PARPRIORUNC[11] = 0.5*8.534234 # ; PARPRIORWEIGHT[11] = noyears # Ceff: derived from multiple trait values from Kattge et al., (2011)
PARPRIORS[12] = 403.87042     ; PARPRIORUNC[12] = 343.6277     # Day of year for sowing (applied as modulus 365.25)
PARPRIORS[13] = 123.09865     ; PARPRIORUNC[13] = 45.06803     # Phenological heat units for seed emergence (aka growing degree days)
PARPRIORS[14] = 230.11399     ; PARPRIORUNC[14] = 161.8359     # Growing season length sowing->harvest days
#      #PARPRIORS[12]=OBS$planting_doy       ; PARPRIORUNC[12]=OBS$planting_doy_unc # Sow day of year, applied as p12%%365.25
#      #PARPRIORS[14]=OBS$harvest_doy        ; PARPRIORUNC[14]=OBS$harvest_doy_unc  # Growing season length sowing->harvest days
# if (max(OBS$LAI) > 0) {
#   # Prior on canopy N derived from ATEC experiment assuming max LAI is related to canopy N
#   PARPRIORS[15] = min(7.0,max(OBS$LAI) * 0.488 + 2.5131)
#   PARPRIORUNC[15] = 1.2 # mean confidence interval of linear regression for LAI ranges 1-6
#   # PARPRIORWEIGHT[15] = noyears
# } else {
#   PARPRIORS[15] = 4.088        ; PARPRIORUNC[15] = 0.6052851    # Constant for canopy N dilution model (gN/m2leaf)
# }
PARPRIORS[16] = -0.0252      ; PARPRIORUNC[16] = 0.00092      # Coefficient relating foliar C to N dilution
PARPRIORS[17] = 25.89064     ; PARPRIORUNC[17] = 12.45272     # Leaf carbon per unit area (gC/m2)
#      #PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
PARPRIORS[18] = 3.23791      ; PARPRIORUNC[18] = 7.735398     # initial labile (gC/m2)
PARPRIORS[19] = 2.25043      ; PARPRIORUNC[19] = 3.543016     # initial foliage (gC/m2)
#      #PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
PARPRIORS[20] = 2.263929     ; PARPRIORUNC[20] = 3.550366     # initial fine roots (gC/m2)
#      #PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
PARPRIORS[21] = 2.233714     ; PARPRIORUNC[21] = 3.542520     # initial stem (gC/m2)
#      #PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
PARPRIORS[22] = 3.125807     ; PARPRIORUNC[22] = 7.797632     # initial litter (gC/m2)
#      #PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
# PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
# the prior on SOM comes from somewhere else
PARPRIORS[23] = 15338.74; PARPRIORUNC[23] = 3047.499
PARPRIORS[24] = 0.726655     ; PARPRIORUNC[24] = 4.105285     # initial autotrophic (gC/m2)
PARPRIORS[25] = 0.320946     ; PARPRIORUNC[25] = 0.786673     # initial storage organ (gC/m2)
PARPRIORS[26] = 275.2997     ; PARPRIORUNC[26] = 5.533301     # min temperature for development (K)
PARPRIORS[27] = 301.8037     ; PARPRIORUNC[27] = 15.51328     # max temperature for development (K)
PARPRIORS[28] = 286.6288     ; PARPRIORUNC[28] = 10.69780     # optimum temperature for development (K)
PARPRIORS[29] = 270.3326     ; PARPRIORUNC[29] = 4.524566     # min temperature for vernalisation (K)
PARPRIORS[30] = 288.5911     ; PARPRIORUNC[30] = 6.244053     # max temperature for vernalisation (K)
PARPRIORS[31] = 277.6463     ; PARPRIORUNC[31] = 4.458580     # optimum temperature for vernalisation (K)
PARPRIORS[32] = 8.172631     ; PARPRIORUNC[32] = 5.096538     # critical photoperiod for development to begin (days)
PARPRIORS[33] = 0.212208     ; PARPRIORUNC[33] = 0.218979     # sensitivity to photoperiod above critical threshold
PARPRIORS[34] = 0.028113     ; PARPRIORUNC[34] = 0.090153     # potential turnover rate of labile (fraction/day)
PARPRIORS[35] = 0.915873     ; PARPRIORUNC[35] = 0.199231     # autotrophic turnover rate (/day)
PARPRIORS[36] = 25.17841     ; PARPRIORUNC[36] = 71.154717    # fine root biomass (g/m2) needed to reach 50 % of maximum depth
PARPRIORS[37] = 2.188784     ; PARPRIORUNC[37] = 4.002502     # maximum rooting depth (m)

PARPRIORS[PARPRIORS == 0] = NA

ranges_df_with_priors = cbind(ranges_df, "prior_mean" = PARPRIORS, "prior_unc" = PARPRIORUNC)
ranges_df_with_priors = ranges_df_with_priors %>% mutate("prior_min" = prior_mean - 3*prior_unc, 
                                                         "prior_max" = prior_mean + 3*prior_unc)
ranges_df_with_priors = ranges_df_with_priors %>%   mutate(
  potential_min = if_else(par_min > prior_min, par_min, prior_min),
  potential_max = if_else(par_max < prior_max, par_max, prior_max)
)

ranges_df_with_priors[is.na(ranges_df_with_priors$potential_min), "potential_min"] = ranges_df_with_priors[is.na(ranges_df_with_priors$potential_min), "par_min"] 
ranges_df_with_priors[is.na(ranges_df_with_priors$potential_max), "potential_max"] = ranges_df_with_priors[is.na(ranges_df_with_priors$potential_max), "par_max"] 
ranges_df_with_priors$prior_ranges_differ = (ranges_df_with_priors$par_min != ranges_df_with_priors$potential_min | ranges_df_with_priors$par_max != ranges_df_with_priors$potential_max)

