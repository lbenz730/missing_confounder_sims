# Missing Confounder Simulations

Repository for simulations exploring causal inference methods with partially missing confounders. Comparison of complete case missing at random (CCMAR) based estimators from (Levis et al., CJS, 2024) as well as ad-hoc but reasonable methods that rely on imputation of missing data + standard causal inference methods. For more information, please see our paper 

Benz, L., Levis, A., and Haneuse, S. "Comparing Causal Inference Methods for Point Exposures with Missing Confounders: A Simulation Study." _Under Review_, 2024. [[Pre-Print](https://arxiv.org/abs/2407.06038)]

---

## R Scripts: 
* __causal_sim.R__: Wrapper to run simulation and save data
* __generate_data__: 
    * Function to generate data for simulations
    * Script also contains estimations of models for $\mu, \eta, \lambda, \pi$ along with hard-coded covariate shifts/amplifications
* __specify_inputs.R__: Script to specify and save out simulation parameters
* __estimate_treatment_effects.R__: Script to implement various causal inference estimation methods
* __helpers.R__: File containing useful helper functions used in various scripts
* __analysis.R__: Script for analyzing simulation results
* __levis_estimators.R__: Script to compute CCMAR estimators from Levis et al. paper
* __levis_helpers.R__: Helper functions in script to compute  CCMAR estimators from Levis et al. paper
* __compute_truth.R__: Function to compute true ATE for various simulation settings.


## Directories
* __inputs/__: Simulation inputs, in .rds file format
* __outputs/__: Simulation results, which contain simulation id in name of file.
* __eda/__: Directory of exploratory work better understanding skew in CCMAR IWOR estimator.
* __paper/__: Directory of scripts to create LaTeX tables and figures for paper.


## Data Driven Simulations
`sim_id` corresponds to the order in which we ran data driven simulations. Below we provide a table that shows how simulations referenced in the paper correspond to our simulations in the __input/__ and __output/__ directories. Note that `sim_id` = 19/20/21 correspond to simulation 4 presented in the paper, where `sim_id = 19` runs the Levis estimators with parametric models for the nuisance functions while `sim_id = 20`/`21` runs the Levis estimators with semi-parametric models for the nuisance functions (all other estimators are identical). The fully flexible simulation with 

| Paper Simulation | sim_id |
|:----------------:|:------:|
|        1         |   4    |
|        2         |   6    |
|        3         |   16   |
|        4         |   19   |
|        4         |   20   |
|        4         |   21   |
|        5         |   1    |
|        6         |   2    |
|        7         |   3    | 
|        8         |   5    |
|        9         |   7    |
|        10        |   8    |
|        11        |   9    |
|        12        |   10   |
|        13        |   11   |
|        14        |   12   |
|        15        |   13   |
|        16        |   14   |
|        17        |   15   |
|        18        |   17   |
|        19        |   18   |


## Fully Flexible Simulation
* __fully_flexible/__: A different simulation study of Levis estimators under fully non-parametric modeling
    * __create_datasets.R__: Save out `n_sims` $\times$ `n_folds` datasets (each a fold of a dataset)
    * __estimate_nuiscance_fx.R__: Use flexible methods to estimate all component nuisance functions (analagous for parametric)
    * __subject_IF_contributions.R__: Compute all the subject contributions to the IF based estimator (analagous for parametric)
    * __consolidate_IF_contributions.R__: Consolidate subject IF contributions to a single estimator
