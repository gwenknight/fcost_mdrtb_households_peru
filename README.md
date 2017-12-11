###***** fcost_mdrtb_households_peru *****###
Code for estimating fitness costs to MDR-TB in households in Lima, Peru

## Abbreviations:
“Bi”: burn in of 10years included
“Fp1”: model 1 (transmission cost)
“Fp2”: model 2 (progression cost)
“Fp3”: model 3 (combined cost)
“Bi30”: sensitivity analysis of burn-in of 30years
“lhs”: latin hypercube sampling code


## Main model code:
bi_gillespie_simple_peru_code_withRecovered_allM2.R
The natural history model is included in here. 

## LHS functions:
bi_fp#_mcmc_peru_lhs_c.R
The initial sampling for the standard deviation for the subsequent MCMC sampler is given here. 

## MCMC functions:
bi_mcmc_functions.R
The functions for running the MCMC sampler with burn-in are here.

bi_fp#_look_at_output_parallel.R
Checks of the MCMC sampler with plots are here. 

## Look at output
bi_fp#_plotoutput.R
Plots are generated with this function. 
