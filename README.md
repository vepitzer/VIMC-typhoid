# VIMC-typhoid
MATLAB code for VIMC typhoid transmission model. Data inputs not included.

The main file for running the model is "vax_sims_vimc_full.m". 
The differential equations for the model are defined in "fn_SIR_vaccine_vnv_strata.m" and "fn_SIR_vaccine_vnv_strata_h20.m". 
  The former is used for the burn-in period while the latter includes a decrease in the transmission rate over time associated with improvements in WASH.
The "vimc_generate_tables.m" files saves the model output to the VIMC-formatted output csv files. 
  (Note: this took a long time to run and generates >10GB files (x3) for the stochastic output (200 runs for 93 countries))
  
 NOTES: 
 Future changes to code should include:
 -Change to loop structure. Right now, different stochastic samples are generated for each country (and in each loop); this is not necessary for most parameters.
 -Limit what is saved to the "output2" structure. For the full runs, this variable becomes too large to save. Not all of the components of the structure are necesssary.
