# History matching with EDCs

Repo containing the code to reproduce the results from paper on history matching with EDCs.
Additionally to the code provided below, to generate DALEC Crop runs, we require to have a local copy of the CARDAMOM repository (https://github.com/GCEL/CARDAMOM), 
and of DALEC Crop (https://github.com/GCEL/DALEC-CropN)

The majority of the history matching code is based on the R package hmer: Iskauskas A, Vernon I, Goldstein M, Scarponi D, McCreesh N, McKinley TJ, White RG
  (2024). "Emulation and History Matching Using the hmer Package." _Journal of
  Statistical Software_, *109*(10), 1-48. doi:10.18637/jss.v109.i10
  <https://doi.org/10.18637/jss.v109.i10>.

## Requirements
- R
- linux computer/ server

## Contributors
- Nina Fischer (School of Maths, University of Edinburgh)
- Luke Smallman (School of Geosciences, University of Edinburgh)

## This repository containts the following files: 
### Code  
- DALEC crop ranges and priors (for CARDAMOM) (parameter_ranges.R)
- code to run the history matching for the synthetic experiment (synth_example_hm_ranges.R)
- code to run the history matching for the ATEC experiment (ATEC_year1_hm_ranges.R and ATEC_year2_hm_ranges.R)
- code to apply the Ecological and Dynamic Constraints (EDCs) (edc1.R and edc2.R)
- code to build the initial Latin Hypercube Design to train the first wave of emulators (build_initial_LHD.R)
- code to propose new points at each wave (point_proposal.R)
- code to run DALEC Crop (generate_input_output_DALEC_Crop.R). This also requires loading of the CARDAMOM functions in the background: ./R_functions/load_all_cardamom_functions.r from the CARDAMOM repo)

### Data
- synthetic observations for synthetic experiment (synth_obs_feb.RDS)
- ATEC Trial Plot data (trial_plot_data_2017_2018.csv) from https://datashare.ed.ac.uk/handle/10283/3697 (Florence et al., 2021, https://doi.org/10.3390/agriculture11030258)

## Results 
We can provide the following results (too big to upload here)
- synthetic experiment results (in synthetic_data.zip): For the following history matching experiments: 1a) history matching with EDCs and 10\% model discrepancy (HM EDCs and md) and 1b) history matching with EDCs and 0\% model discrepancy (HM EDCs) and 2) history matching without EDCs and 10\% model discrepancy. All .RDS files provide each wave's emulators, the input-DALEC Crop output data for each wave, information on any priors used (none) and for which parameters the ranges were reduced at the final wave.
- ATEC experiment results (in ATEC_2016.zip and ATEC_2017.zip). All .RDS files provide the final wave's input-DALEC Crop output data, information on any priors used (none) and for which parameters the ranges were reduced at the final wave.
- compiled versions of DALEC Crop (written in fortran) to run DALEC Crop via generate_input_output_DALEC_Crop.R: run_ATEC_year1/infofile.RData and run_ATEC_year2/infofile.RData

