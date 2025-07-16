library(hmer)
library(tidyr)
library(dplyr)
library(ggplot2)

# read in targets
trial_plot_data_2017_2018 <- readr::read_csv("synth_experiment/synth_jan/trial_plot_data_2017_2018.csv",
                                             show_col_types = FALSE)
trial_plot_data_2017_2018 = trial_plot_data_2017_2018 %>% 
  mutate("LAI_unc" = 0.12*LAI, 
         "yield_conv" = 0.48*100*YIELD) 

# read in initial ranges for LHS
# also get ranges_list
# declare before which variables are fixed
source("synth_experiment/synth_jan/utils/parameter_ranges.R") # gives us original_ranges, ranges_df and ranges_df_with_priors
fixed_vars = c("sow_day", "harvest_length")

find_ranges_fct = function(input_data, original_ranges){
  input_data = input_data[names(original_ranges)]
  ranges <- setNames(
    map(
      names(original_ranges), function(nm) {
        c(
          max(original_ranges[[nm]][1], min(input_data[,nm])-0.01*diff(range(input_data[,nm]))),
          min(original_ranges[[nm]][2], max(input_data[,nm])+0.01*diff(range(input_data[,nm])))
        )
      }
    ),
    names(original_ranges))
  ranges <- hmer:::convertRanges(ranges)
  return(ranges)
}

# get all functions
source("synth_experiment/synth_jan/utils/edc1.R")
source("synth_experiment/synth_jan/utils/edc2.R")
source("synth_experiment/synth_jan/utils/build_initial_LHD.R")
source("synth_experiment/synth_jan/utils/point_proposal.R")
source("synth_experiment/synth_jan/utils/generate_input_output_DALEC_Crop.R")

source("./R_functions/load_all_cardamom_functions.r")

#################### which plot are we running?
year = 2018
project_name = "test"
if(year == 2017){
  obs_dates_names = c("date20170420", "date20170516", "date20170530", "date20170612", "date20170626", "date20170710")
  obs_dates = as.Date(c("2017-04-20", "2017-05-16", "2017-05-30", "2017-06-12", "2017-06-26", "2017-07-10"))
  load("CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/run_ATEC_year1/infofile.RData")
} else{
  obs_dates_names =  c("date20180507", "date20180613", "date20180619", "date20180704")
  obs_dates = as.Date(c("2018-05-07", "2018-06-13", "2018-06-19", "2018-07-04"))
  load("CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/run_ATEC_year2/infofile.RData")
  
}

prior_vals_list = list("sow_day" = 273, "harvest_length" = 329)
# get drivers
load("CARDAMOM_outputs/DALEC.A3.C3.H2.M1.#_MHMCMC/ATEC_year2_jan25/RESULTS_PROCESSED/ES_plot1_y2_2016_2018_parameters.RData")
new_PROJECT = PROJECT
library(purrr)

###################

### begin with history matching
n_points_train = 2*10*36 # start with a bigger LHD for wave 0
# more_fixed_vars = NULL
# more_fixed_vars = c("C_fol", "C_roots", "C_stem", "C_auto", "C_stor")
# more_fixed_vars = c("C_fol", "C_roots", "C_stem", "C_litter", "C_auto", "C_stor")
# more_prior_vars = c("C_som")
more_fixed_vars = NULL
more_prior_vars = NULL
more_prior_vars = c("NUE", "N_slope", "LMA", "C_som")
prior_vals_list = list("sow_day" = 273, "harvest_length" = 329) # , 
                       # "C_fol" = ranges_df_with_priors[ranges_df_with_priors$par_name == "C_fol","prior_mean"], 
                       # "C_roots" = ranges_df_with_priors[ranges_df_with_priors$par_name == "C_roots","prior_mean"], 
                       # "C_stem" = ranges_df_with_priors[ranges_df_with_priors$par_name == "C_stem","prior_mean"], 
                       # "C_litter" = ranges_df_with_priors[ranges_df_with_priors$par_name == "C_litter","prior_mean"], 
                       # "C_auto" = ranges_df_with_priors[ranges_df_with_priors$par_name == "C_auto","prior_mean"], 
                       # "C_stor" = ranges_df_with_priors[ranges_df_with_priors$par_name == "C_stor","prior_mean"])

for (plot_id in c(49:50)) {
# c(29,31,36,43,48,4,9,11,16,23,28,35,40,42,47)
  project_hm_name = paste0("HM_year2_no_ranges_june_plot", plot_id)
  print(project_hm_name)

# grab the right targets
LAI_plot = trial_plot_data_2017_2018 %>% filter(PLOT_ID == plot_id, YEAR == year) %>% arrange(DOY) %>% select(LAI) %>% pull()
LAI_unc_plot = trial_plot_data_2017_2018 %>% filter(PLOT_ID == plot_id, YEAR == year) %>% arrange(DOY) %>% select(LAI_unc) %>% pull()
targets = lapply(seq_along(LAI_plot), function(i) list(val = LAI_plot[i], sigma = LAI_unc_plot[i]))
names(targets) = obs_dates_names
targets$ra_GPP_ratio = list(val = 0.54, sigma = 0.12*(1/sqrt(3)))
targets$yield_GPP_ratio = list(val = 0.38, sigma = 0.087*(1/sqrt(3)))
targets$EDC_SOM1 = list(val = 0, sigma = log(10)/3)
targets$EDC_SOM2 = list(val = 0, sigma = 0.2/3)
targets$DS_max = list(val = 1.9, sigma = 1)

discs_prop = 0.1
model_sigma_discs = discs_prop*sapply(targets, "[[", 1) 
model_discs = purrr::map(seq_along(model_sigma_discs),
                         ~list(internal = model_sigma_discs[[.]], external = NULL))
max_LAI =max(unlist(lapply(targets[1:4], "[[", 1)))
# potential_ranges = truncnorm::qtruncnorm(p = c(0.005, 0.995), a = 2.5, b = 7, mean = max_LAI, sd = 1.2)
# if(potential_ranges[1] > 2.5){
#   ranges_df_with_priors[15, "par_min"] = potential_ranges[1]
# }
# if(potential_ranges[2] < 7){
#   ranges_df_with_priors[15, "par_max"] = potential_ranges[2]
# }

# other info that will be useful later
plot_info_list = list(
  "plot_id" = plot_id,
  "true_yield" = trial_plot_data_2017_2018 %>% filter(PLOT_ID == plot_id, YEAR == year) %>% select(yield_conv) %>% unique() %>% pull(),
  "N_level" = trial_plot_data_2017_2018 %>% filter(PLOT_ID == plot_id, YEAR == year) %>% select(FERTILISER_LEVEL) %>% unique() %>% pull()
)

# generate first training data set
inputs_train_wave0 = initial_LHD(n_points_train, TRUE, other_fixed_vars = more_fixed_vars, pars_with_priors = more_prior_vars,
                                 fixed_vars = fixed_vars, ranges_df_with_priors = ranges_df_with_priors)
initial_ranges_list = inputs_train_wave0$wave0_ranges
initial_ranges_list <- lapply(initial_ranges_list, function(df) {
  colnames(df) <- c("min", "max")
  df
})

initial_ranges_list <- lapply(initial_ranges_list, function(df) {
  as.numeric(unlist(df)) # Convert data frame to numeric vector
})

# generate first testing data set
n_points_test = 20
inputs_test_wave0 = initial_LHD(n_points_test, TRUE, other_fixed_vars = more_fixed_vars, pars_with_priors = more_prior_vars,
                                fixed_vars = fixed_vars, ranges_df_with_priors = ranges_df_with_priors)

# run DALEC Crop for inputs to get training and testing data set
data_wave0 = list("training_data" = generate_input_output_DALEC_Crop(inputs_train_wave0$initial_LHD, 
                                                                     obs_dates = obs_dates, obs_dates_names = obs_dates_names, 
                                                                     fixed_pars = prior_vals_list, 
                                                                     year = 2018),
                  "testing_data" = generate_input_output_DALEC_Crop(inputs_test_wave0$initial_LHD, 
                                                                    obs_dates = obs_dates, obs_dates_names = obs_dates_names, 
                                                                    fixed_pars = prior_vals_list, year = 2018)
)

output_names = colnames(data_wave0$training_data)[!colnames(data_wave0$training_data) %in% vars_DALEC.A3.C3.H2.M1]
output_names = output_names[!output_names %in% "harvest_vector"]

# start entering loop for multiple waves
ems_list = list()
data_list = list("wave0" = data_wave0)

ranges_list = list()
ranges_list[[1]] = initial_ranges_list

n_points_train = 10*length(initial_ranges_list)
for(wave in c(1:3)){ # remember that we need to start this loop at 1, 
  
  # build emulator
  new_ems = emulator_from_data(
    input_data = data_list[[wave]]$training_data, 
    output_names = output_names, 
    specified_priors = list(hyper_p = rep(0.6,length(targets))),
    ranges = ranges_list[[wave]],
    beta.var = TRUE,
    check.ranges = FALSE,
    # corr_name = "matern", 
    targets = targets,
    add = TRUE, 
    order = 2, 
    discrepancies = model_discs
  )
  ems_list = c(list(new_ems), ems_list)
  
  # use testing data from previous wave to look at validation diagnostics
  validation_diagnostics(emulators = ems_list[[wave]], 
                         targets = targets,
                         validation = data_list[[wave]]$testing_data)
  
  # propose new training data
  new_training_inputs = generate_new_design_edcs(all_ems = ems_list, n_requested_points = n_points_train, targets = targets,  
                                                 cutoff = 3, n = 1)
  
  # propose new testing data
  new_testing_inputs = generate_new_design_edcs(all_ems = ems_list, n_requested_points = n_points_test, targets = targets,  
                                                cutoff = 3, n = 1)
  
  data_list[[wave+1]] = list("training_data" = generate_input_output_DALEC_Crop(new_training_inputs, 
                                                                       obs_dates = obs_dates, obs_dates_names = obs_dates_names, 
                                                                       fixed_pars = prior_vals_list, 
                                                                       year = 2018),
                    "testing_data" = generate_input_output_DALEC_Crop(new_testing_inputs, 
                                                                      obs_dates = obs_dates, obs_dates_names = obs_dates_names, 
                                                                      fixed_pars = prior_vals_list, year = 2018)
  )
  
  ranges_list[[wave + 1]] = find_ranges_fct(new_training_inputs, original_ranges = ranges_list[[wave]])
  
  print(mean(data_list[[wave+1]]$training_data$harvest_vector))
  print(median(data_list[[wave+1]]$training_data$harvest_vector))
  
}

wave = wave + 1

# build emulator
new_ems = emulator_from_data(
  input_data = data_list[[wave]]$training_data, 
  output_names = output_names, 
  specified_priors = list(hyper_p = rep(0.6,length(targets))),
  ranges = ranges_list[[wave]],
  beta.var = TRUE,
  check.ranges = FALSE,
  # corr_name = "matern", 
  targets = targets,
  add = TRUE, 
  order = 2, 
  discrepancies = model_discs
)
ems_list = c(list(new_ems), ems_list)

# use testing data from previous wave to look at validation diagnostics
validation_diagnostics(emulators = ems_list[[wave]], 
                       targets = targets,
                       validation = data_list[[wave]]$testing_data)

# propose new training data
new_training_inputs = generate_new_design_edcs(all_ems = ems_list, n_requested_points = 300, targets = targets,  
                                               cutoff = 3, n = 1)

data_list[[wave+1]] = list("training_data" = generate_input_output_DALEC_Crop(new_training_inputs, 
                                 obs_dates = obs_dates, obs_dates_names = obs_dates_names, is_final_wave = TRUE, 
                                 fixed_pars = prior_vals_list, 
                                 year = 2018))

# compile all data into one rds file
# compile all data into one rds file
this_plots_data = list(# "emulators" = ems_list, 
  "data" = data_list[[wave+1]], 
  "plot_info" = plot_info_list, 
  "NROY_reduced" = ranges_list[[4]][names(initial_ranges_list)[!mapply(identical, initial_ranges_list, ranges_list[[4]])]]
)

LAI_hm_long = data_list[[5]]$training_data$LAI_curves %>%
  pivot_longer(., -date, values_to = "LAI", names_to = "TS.nu")

targets_LAI = targets[1:4]

# need to combine the LAIs in a dataframe for easier plotting
LAI_df = data.frame("val" = unlist(lapply(targets_LAI, FUN = "[[", 1)),
                    "sigma" = unlist(lapply(targets_LAI, FUN = "[[", 2)),
                    "date" = obs_dates)

ggplot() +
  geom_line(data = LAI_hm_long %>% group_by(date) %>% summarise(av_LAI = mean(LAI)),
            aes(x = date, y = av_LAI, color = "HM"), size=2) +
  geom_ribbon(data = LAI_hm_long %>% group_by(date) %>%
                mutate(min_LAI = mean(LAI) - 2*sd(LAI),
                       max_LAI = mean(LAI) + 2*sd(LAI)),
              aes(ymin = min_LAI, ymax = max_LAI, x = date), colour = "#648FFF", alpha = 0.1, linetype=2) +
  geom_errorbar(data = LAI_df, aes(x = as.Date(date), ymin = val - 2*sigma,
                                   ymax = val + 2*sigma)) +
  geom_point(data = LAI_df, aes(x = as.Date(date), y = val))
# saveRDS(this_plots_data, paste0("~/scratch/ATEC_year1/", project_hm_name, ".rds"))

df = data_list[[5]]$training_data$training_data %>% 
  mutate(group = ifelse(harvest_vector < 200, "<200", ">=200")) %>% 
  pivot_longer(cols = vars_DALEC.A3.C3.H2.M1, names_to = "variable", values_to = "value")

ggplot(df, aes(x = value, fill = group)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~variable, scales = "free") +
  theme_minimal()

print(this_plots_data$data$training_data$training_data$harvest_vector %>% mean())


# check which parameters are getting reduced:
diffs <- names(initial_ranges_list)[!mapply(identical, initial_ranges_list, ranges_list[[4]])]
ranges_list[[4]][diffs]


saveRDS(this_plots_data,
        paste0("~/scratch/ATEC_year2/", project_hm_name, ".rds"))
ems_list = list()
}


