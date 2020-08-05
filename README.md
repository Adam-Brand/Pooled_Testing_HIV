# Pooled_Testing_HIV
 This project contains all of the code used to produce the results in the paper titled "Prediction-Driven Pooled Testing Methods: Application to HIV Treatment Monitoring in Kalisizo, Uganda".

Below are descriptions of programs, documentation and datasets used on this project.


The real data collected from Uganda is not included in the repository. All of the code for cleaning/formatting the Uganda data is included, but will not run as expected without the course data file.

In order to replicate the simulation results presented in the paper:

1) run data_gen.R. This will create 6 datasets needed to run the simulations. There will be errors as part of the data_gen.R program sources the data cleaning program for the real data, which is not included in this repository. There is also a shiny program we used for visualizing distributions of VLs. The shiny folder is called shinypopfit, and is included in this repository. Downloading/using this shiny program is optional though. It is not needed for reproducing results.

    a) the datasets will save in the SimData folder

    b) the name of the 6 datasets as R objects are:

      simdata - simulation data with SD=1.0 used for method evaluation

      simdata_train - sim data with SD=1.0 used as training set to get estimated  model betas

      simdata_rev - this is the simdata_train dataset with direction of covariate association with VL reversed

      simdata2 - same as simdata, but with SD=0

      simdata2_train - same as simdata_train, but with SD=0

      simdata2_rev - same as simdata_rev, but with SD=0

2) run program sim_data_betas.R. This program uses the training set for simulation data, simdata_train, to estimate the betas using different models. For the simulations in the paper, we used only the estimated betas from the data where SD=1.0 to predict VL using ridge regression. We did not set seeds for ridge regression, so betas may vary slightly. The betas we used are included in the method_eval programs (more on these below).

3) Run a method_eval program from the 'Method_eval' folder. Each program sources the program containing the needed functions, method_eval_source. This program must be downloaded and the working directory set to point to it before running.

    a) Each method_eval program evaluates the methods on 50,000 samples (500, 10 x 10 matrices) in a variety of scenarios.

    b) Which scenario each statement evaluates should be recognizable by the name of the result file. The statement begins with a readRDS command which tells us which data we are using, i.e., SD=1.0 or SD=0. Then the estimated betas are set which are the betas used to predict the VL based on the 2 covariates. (This is where the betas are recorded that we obtained from sim_data_betas program.)

    c) Within each evaluation scenario, we set the number of matrices, matrix size, SE (the ME on log10 scale in the paper). prec and precrd have to do with prediction fitting precision variables, and we leave them at 10 and 20, respectively. Cutoff is left at 1000 and represents the VL cutoff for 'failure'. We leave tstperd at 5, lowlimit at 50 and filltype="rnd" for the evaluations in the paper. There are 3 main evaluation functions within the method_eval programs:

        i) pool.alg.cov - this function evaluates the same matrices for all methods excluding the Hypred method. The output is a table of results where each row represents the performance of classifying 100 subjects.

        ii) hypred - this function does the same as pool.alg.cov except for only the Hypred method. There are additional arguments which are top_percent and bot_percent. These are the cutoffs (between 0 and 1) for the top tier and the bottom tier. The Hypred function simulates individual testing for the top tier (we set it at the top 10% of predicted VLs). Hypred also simulates each of the pooled testing methods discussed in the paper for both the middle and bottom tier (we reported MiniPred for both tiers based on simulation results). To get results for 50,000 samples using Hypred, you must choose which method to use for each the middle and bottom tier, and compile those results with individual testing.

        iii) hypred_uganda - the hypred method specifically for the uganda data. This was written to make sure that the number of subjects in the top tier/bottom tiers were divisible by 100. We rounded the top tier down, so there are less than 10% in the top tier (300 out of 3600).

4) Run a results_anal program. These programs read in the results datasets produced by the method_eval programs, and creates LaTeX code for tables of summary stats for method performance. There are 3 of these programs:

    a) results_anal.R - this program creates LaTeX code for 16 tables; 8 using SD=1 data and 8 using SD=0 data. All of these are reported in tables 1 and 2.

    b) results_anal_hypred.R - creates LaTeX code for the Hypred method. The tables created are risk tier specific, i.e., one for just the middle tier and another for just the bottom tier. It also compiles results from all 3 tiers and outputs mean performance summaries which are reported in tables 1-3. The means reported in tables 1-3 come from the function calls comb_hypred at the bottom of this program file.

    c) results_anal_uganda.R - creates LaTeX code for summary tables using the uganda data reported in table 3. The mean performance measures reported for HyPred in table 3 are output in resulta_anal_hypred.R

This will replicate the simulation results. Again, results from the Uganda data cannot be produced without the actual data which we cannot upload to the repository.

Below are specifics on the data and the programs included in this repository.

Data Dictionary (DD):
  The data dictionary is an Excel workbook with a separate sheet for the datasets included in this repository as well as the real Uganda data which is not. The dictionary for each dataset contains a row for each variable in that dataset and the source for the variable. If it is a derived variable it also includes the derivation for that variable. The data dictionary includes the simulation datasets used in Section 4 as well as the raw, intermediate and final datasets containing real-world data from Uganda used in Section 5. The raw dataset from Uganda is not included in the repository.

Programs:

  raw_data_clean.R: this program reads in the raw data from Uganda, cleans the variables and adds some derived variables as described in the DD.

  training_set.R: takes the cleaned Uganda data generated by raw_data_clean, transforms it into long format, filters out records by time (those prior to June 30th, 2010) and adds variables as described in the DD.

  test_set: does the same as training_set, but filters for all records after June 30th, 2010.

  uganda_eval_data.R - uses the predictive model obtained by model_pred.R to predict the VLs in the test_set for the uganda data, and formats the test_data for use in the method evaluation.

  data_gen.R: generates 6 simulated datasets used in Section 4 of the paper as described in the DD.

  method_eval_source.R: this program contains code for all of the functions needed for the method evaluation, and must be sourced before running the code to evaluate the methods.

  Method_eval folder programs: there are 4 programs in this folder - each one evaluates the methods in different scenarios.

    method_eval3_run2.R - reads in the simulated data with SD=1.0 and evaluates the methods (except Hypred) in 8 scenarios

    method_eval4_run2.R - same as method_eval3_run2.R, but using the data with SD=0

    method_run2_hypred.R - evaluates the HyPred method in the same scenarios as method_eval3_run2.R and method_eval4_run2.R. In order to make sure the same data is used, I extracted the 50,000 records used in the 500 matrices for the other methods, and applied them to the HyPred method to ensure they used the same records. That datasets are called SD1.0_data_500_run2.R and SD0_data_500_run2.R, and are saved here Pooled_Testing_HIV\SimData\Records_used_in_500_matrices_hypred

    method_run2_uganda.R - evaluates the methods using the data from Uganda in a variety of measurement errors.

  model_pred.R - applies a variety of regression techniques for prediction on the training set of the uganda data. We ended up using ridge regression modelling single failure.

  sim_data_betas.R - applies regression techniques on the simulated training sets of 5,000 records each to get estimated betas to use in the method evaluation of the simulated data.

  results_anal.R -compiles the results of the method evaluation for the simulations (minus HyPred) and puts them into tables; outputs latex code for the summary tables.

  results_anal_hypred.R - does the same as results_anal.R, but for the HyPred method using the simulated dataset

  results_anal_uganda.R - does the same as results_anal.R but for all methods using the real uganda data.
