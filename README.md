# History matching with EDCs

Repo containing the code to reproduce the results from the paper *History matching with ecological constraints for efficient calibration of the DALEC-Crop simulator*.

- `revision_2026/` contains the code for the revised manuscript. All results, tables and figures in the revised paper were produced with this code.
- The files at the top level are the code of the original submission, kept for reference (also available as the release `v1-original-submission`).

Additionally to the code provided below, to generate DALEC Crop runs, we require a local copy of the CARDAMOM repository (https://github.com/GCEL/CARDAMOM, version of August 2024) and of DALEC Crop (https://github.com/GCEL/DALEC-CropN). For the revised manuscript, four CARDAMOM files were modified: `cardamom_io.f90`, `MODEL_LIKELIHOOD.f90`, `DALEC.A3.C3.H2.M1.#.f90` and `DALEC.A3.C3.H2.M1.#_R_interface.f90`.

The majority of the history matching code is based on the R package hmer: Iskauskas A, Vernon I, Goldstein M, Scarponi D, McCreesh N, McKinley TJ, White RG
  (2024). "Emulation and History Matching Using the hmer Package." _Journal of
  Statistical Software_, *109*(10), 1-48. doi:10.18637/jss.v109.i10
  <https://doi.org/10.18637/jss.v109.i10>.

## Requirements
- R (we used R 4.4.3 and hmer 1.6.0), with the packages dplyr, tidyr, purrr, tibble, readr, ggplot2, ggpubr and scales
- linux computer/ server

## Contributors
- Nina Fischer (School of Maths, University of Edinburgh)
- Luke Smallman (School of Geosciences, University of Edinburgh)

## The folder revision_2026 contains the following files:
All scripts are run from the CARDAMOM working directory. They expect the code in `synth_experiment/revision/` and write all results to `~/scratch`; both can be changed by setting `code_dir` and `scratch` before running a script.

### Code to run history matching
- DALEC crop ranges and priors (for CARDAMOM) (parameter_ranges.R)
- code to apply the Ecological and Dynamic Constraints (EDCs) (edc1.R and edc2.R)
- code to build the initial Latin Hypercube Design to train the first wave of emulators (build_initial_LHD.R)
- code to propose new points at each wave (point_proposal.R)
- code to run DALEC Crop (generate_input_output_DALEC_Crop.R). This also requires loading of the CARDAMOM functions in the background: ./R_functions/load_all_cardamom_functions.r from the CARDAMOM repo
- the history matching workflow as one function, run_history_matching() (hm_runner.R)
- code to run all synthetic experiments: Exp. 1a, 1b and 2, the sensitivity to the model discrepancy and the repeated runs (run_synthetic_all.R)
- code to run the ATEC experiment for both seasons (run_ATEC_all.R, set `year <- 2017` or `year <- 2018` first)

### Code to analyse the results
- checks and diagnostics for each synthetic run, e.g. whether the true parameters were ruled out at any wave (posthoc_synthetic_analysis.R)
- summaries of all synthetic runs (summarise_synthetic.R)
- sensitivity to the correlation length (theta_all_waves.R)
- time of one DALEC Crop run vs one emulator evaluation (timing_benchmark.R)
- emulator validation and stopping diagnostics for the synthetic experiment (validation_and_stopping.R)
- summaries of the ATEC history matching results (check_existing_results.R)
- CARDAMOM results for the ATEC plots, including R-hat, and the comparison with history matching (summarise_CARDAMOM.R)
- the "All" rows of Table 11 (table11_all_rows.R)

### Code to make the figures
- colours and line types shared by all figures (figure_scheme.R)
- synthetic experiment, Figures 2, 3, 4, 9 and 10 (make_figures_synthetic.R, which uses figure2_synthetic.R and figure4_10_parameters.R)
- ATEC experiment, Figures 5 and 6 (make_figures_ATEC.R, which uses figure5_LAI.R and figure6_yield.R)
- emulator validation, Figures 12 and 13 (validation_and_stopping.R, figure_validation_ATEC.R and redraw_validation_figures.R)
- LAI for further ATEC plots, Figure 15 (lai_panels_overview.R)

The figure scripts save intermediate results (forward_runs_cache.rds, figure2_data.rds, figure4_10_data.rds and the folder figure5_cache), so that figures can be redrawn without running DALEC Crop again. Delete these after rerunning any experiment.

### Order to run the scripts
1. run_synthetic_all.R and run_ATEC_all.R (history matching), and CARDAMOM for the synthetic experiment and the ATEC plots
2. posthoc_synthetic_analysis.R (for each synthetic run), summarise_synthetic.R, theta_all_waves.R, timing_benchmark.R, validation_and_stopping.R, check_existing_results.R, summarise_CARDAMOM.R and table11_all_rows.R
3. make_figures_synthetic.R, make_figures_ATEC.R, figure_validation_ATEC.R, redraw_validation_figures.R and lai_panels_overview.R

## Data
- synthetic observations for synthetic experiment (synth_obs_feb.RDS)
- ATEC Trial Plot data (trial_plot_data_2017_2018.csv) from Edinburgh DataShare, https://doi.org/10.7488/ds/2989 (Florence et al., 2021, https://doi.org/10.3390/agriculture11030258)

## Results
We can provide the following results (too big to upload here):
- synthetic experiment results: for each history matching experiment (1a) history matching with EDCs and 10\% model discrepancy, 1b) history matching with EDCs and 0\% model discrepancy, 2) history matching without EDCs and 10\% model discrepancy, the discrepancy levels 5\%, 20\% and 30\%, and the repeated runs of 1a), an .RDS file with each wave's emulators, the input-DALEC Crop output data for each wave, the parameter ranges at each wave, the numbers of proposed and accepted points and the run times.
- ATEC experiment results: for each plot, an .RDS file with the final wave's input-DALEC Crop output data, the parameter ranges at each wave, the emulator validation results and the numbers of proposed and accepted points.
- compiled versions of DALEC Crop (written in fortran) to run DALEC Crop via generate_input_output_DALEC_Crop.R: run_ATEC_year1/infofile.RData and run_ATEC_year2/infofile.RData

The random seeds are set in run_synthetic_all.R (one per synthetic run) and run_ATEC_all.R (the plot number for each ATEC calibration).
