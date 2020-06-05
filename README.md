# Pooled_Testing_HIV
 This project contains all of the code used to produce the results in the paper titled "Prediction-Driven Pooled Testing Methods: Application to HIV Treatment Monitoring in Kalisizo, Uganda".

Below are descriptions of programs, documentation and datasets used on this project.

IMPORTANT: the programs will not run until the working directories and filenames are set up in each of the programs. Many of the programs refer to/source other programs, so working directories need to be set appropriately to run. Make sure when copying programs to your system that you declare all of the working directories pointing to the correct locations. There are comments noting where to set working directories in the programs. Search for 'setwd' to find them. Also, the method_eval programs need a location to write results to. In the 'write.table' statements in those programs, put the desired file path in the quotations before the name of the result file.

The real data collected from Uganda is not included in the repository. All of the code for cleaning/formatting the Uganda data is included, but will not run as expected without the course data file.

In order to replicate the simulation results presented in the paper:

1) run data_gen.R. This will create 6 datasets needed to run the simulations. There will be errors as part of the data_gen.R program sources the data cleaning program for the real data, which is not included in this repository. There is also a shiny program we used for visualizing distributions of VLs. The shiny folder is called shinypopfit, and is included in this repository. Downloading/using this shiny program is optional though. It is not needed for reproducing results.

    a) there are write.table statements which the user needs to fill in to save the datasets to the location of their choice. They are commented out, so uncomment them when you set the working directory for their saved location.

    b) the name of the 6 datasets as R objects are:

      simdata - simulation data with SD=1.0 used for method evaluation

      simdata_train - sim data with SD=1.0 used as training set to get estimated  model betas

      simdata_rev - this is the simdata_train dataset with direction of covariate association with VL reversed

      simdata2 - same as simdata, but with SD=0

      simdata2_train - same as simdata_train, but with SD=0

      simdata2_rev - same as simdata_rev, but with SD=0

2) run program sim_data_betas.R. This program uses the training set for simulation data, simdata_train, to estimate the betas of different models. For the simulations in the paper, we used only the estimated betas from the data where SD=1.0 to predict VL. We did not set seeds for ridge regression, so betas may vary slightly. The betas we used are included in the method_eval programs (more on these below).

3) Run a method_eval program from the 'Method_eval' folder. Each program sources the program containing the needed functions, method_eval_source. This program must be downloaded and the working directory set to point to it before running. Nothing in method_eval_source needs to be changed.

  a) Each method_eval program evaluates the methods over 50,000 samples (500, 10 x 10 matrices) in a variety of scenarios.

  b) Which scenario each statement evaluates should be recognizable by the name of the result file in the write.table statement below the function call. The statement begins with a read.table command which tells us which data we are using, i.e., SD=1.0 or SD=0. Then the estimated betas are set which are the betas used to predict the VL based on the 2 covariates. (These are where the betas are recorded that we obtained from sim_data_betas program.)

  c) Within each evaluation scenario, we set the number of matrices, matrix size, SE (which is the ME on log10 scale in the paper). prec and precrd have to do with prediction fitting precision variables, and we leave them at 10 and 20, respectively. Cutoff is left at 1000 and represents the VL cutoff for 'failure'. We leave tstperd at 5, lowlimit at 50 and filltype="rnd" for the evaluations in the paper. There are 3 main evaluation functions within the method_eval programs:

    1) pool.alg.cov - this function evaluates the same matrices for all methods excluding the Hypred method. The output is a table of results where each row represents the performance of classifying 100 subjects.

    2) hypred - this function does the same as pool.alg.cov except for only the Hypred method. There are additional arguments which are top_percent and bot_percent. These are the cutoffs (between 0 and 1) for the top tier and the bottom tier. The Hypred function simulates individual testing for the top tier (we set it at the top 10% of predicted VLs). Hypred also simulates each of the pooled testing methods discussed in the paper for both the middle and bottom tier (we reported MiniPred for both tiers based on simulation results). To get results for 50,000 samples using Hypred, you must choose which method to use for each the middle and bottom tier, and compile those results with individual testing.

    3) hypred_uganda - the hypred method specifically for the uganda data. This was written to make sure that the number of subjects in the top tier/bottom tiers were divisible by 100. We rounded the top tier down, so there are less than 10% in the top tier (300 out of 3600).

4) Run a results_anal program. These programs read in the results datasets produced by the method_eval programs, and creates LaTeX code for tables of summary stats for method performance. There are 3 of these programs:

  a) results_anal.R - this program creates LaTeX code for 16 tables; 8 using SD=1 data and 8 using SD=0 data. All of these are reported in tables 1 and 2.

  b) results_anal_hypred.R - creates LaTeX code for the Hypred method. The tables created are risk tier specific, i.e., one for just the middle tier and another for just the bottom tier. It also compiles results from all 3 tiers and outputs mean performance summaries which are reported in tables 1-3.

  c) results_anal_uganda.R - creates LaTeX code for summary tables using the uganda data reported in table 3. The mean performance measures reported for HyPred in table 3 are output in resulta_anal_hypred.R

This will replicate the simulation results. Again, you will not be able to replicate the results from the Uganda data without the actual data which we cannot upload to the repository.

Below are specifics on the data and the programs included in this repository.

Data Dictionary (DD):
  The data dictionary is an Excel workbook with a separate sheet for the datasets included in this repository. The dictionary for each dataset contains a row for each variable in that dataset and the source for the variable. If it is a derived variable it also includes the derivation for that variable. The datasets include the simulation datasets used in Section 4 as well as the raw, intermediate and final datasets containing real-world data from Uganda used in Section 5. The raw dataset from Uganda is not included in the repository.

Programs:
  raw_data_clean.R: this program reads in the raw data from Uganda, cleans the variables and adds some derived variables as described in the DD.

  training_set.R: takes the cleaned Uganda data generated by raw_data_clean, transforms it into long format, filters out records by time (those prior to June 30th, 2010) and adds variables as described in the DD.

  test_set: does the same as training_set, but filters for all records after June 30th, 2010.

  uganda_eval_data.R - uses the predictive model obtained by model_pred.R to predict the VLs in the test_set for the uganda data, and puts the test_set data into format to use in the method evaluation.

  data_gen.R: generates 6 simulated datasets used in Section 4 of the paper as described in the DD.

  method_eval_source.R: this program contains code for all of the functions needed for the method evaluation, and must be sourced before running the code to evaluate the methods.

  Method_eval folder programs: there are 4 programs in this folder - each one evaluates the methods in different scenarios.

    method_eval3_run2.R - reads in the simulated data with SD=1.0 and evaluates the methods (minus Hypred) in 16 scenarios (not all presented in the paper)

    method_eval4_run2.R - same as method_eval3_run2.R, but using the data with SD=0

    method_run2_hypred.R - evaluates the HyPred method in the same scenarios as method_eval3_run2.R and method_eval4_run2.R does for the other methods. In order to make sure the same data is used, I extracted the 50,000 records used in the 500 matrices for the other methods, and applied them to the HyPred method to ensure they used the same records. That data is saved here Pooled_Testing_HIV\SimData\Records_used_in_500_matrices_hypred

    method_run2_uganda.R - evaluates the methods using the data from Uganda in a variety of measurement errors.

  model_pred.R - applies a variety of regression techniques for prediction on the training set of the uganda data. We ended up using ridge regression modelling single failure.

  sim_data_betas.R - applies regression techniques on the simulated training sets of 5,000 records each to get estimated betas to use in the method evaluation of the simulated data.

  results_anal.R -compiles the results of the method evaluation for the simulations (minus HyPred) and puts them into tables; outputs latex code for the summary tables.

  results_anal_hypred.R - does the same as results_anal.R, but for the HyPred method using the simulated dataset

  results_anal_uganda.R - does the same as results_anal.R but for all methods using the real uganda data.
