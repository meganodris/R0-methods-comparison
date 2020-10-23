## R0-methods-comparison

Code and data to reproduce the analysis in "A comparative analysis of statistical methods to estimate the reproduction number in emerging epidemics with implications for the current COVID-19 pandemic", https://doi.org/10.1093/cid/ciaa1599.

***Code***
- 'run_analysis.R' is the core R script to reproduce all key analyses.
- 'methods.R' contains functions for each method that will return R0 estimates and associated 95% confidence intervals.
- 'fit_R0_seq.R' contains a function that will fit all methods to an increasing number of time points, in sequence, in the case time series.
- 'data_simulation.R' contains a function to simulate epidemic curves from an SEIR model, returning 3 versions of the simulated data to reflect varying levels of noise in the case time series:
   - no noise
   - Poisson noise
   - Negative Binomial noise
- 'assess_performance.R' contains functions to calculate and plot performance metrics from simulation results.
- 'R0_plots.R' contains a function to plot average R0 estimates obtained at sequential time points against the corresponding case time series.


***Datasets***
- 'ZikaCases_LatinAmerica&Caribbean.RDS' contains national Zika case surveillance data for 37 countries.
