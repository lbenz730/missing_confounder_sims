library(readr)

### Simulation 1
### Original parameters used in Levis et. al paper
shifts <- 
  list('eta' =  list('GENDER' = 0.35, 
                     'raceeth_hispanic' = 0.25,
                     'GENDER:raceeth_hispanic' = 0.4),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325,
                     'A:raceeth_hispanic' = 0.09,
                     'raceeth_hispanic' = 0.02,
                     'A:bbl' = -0.003),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =   list('A' = -0.6,
                         'ptwc' = -1.3,
                         'A:ptwc' = -0.4)
  )

model_types <- 
  list('mu' = 'levis_paper',
       'eta' = 'levis_paper',
       'pi' = 'levis_paper',
       'lambda' = 'levis_paper',
       'lambda_shape' = 'levis_paper')

params <- 
  list('sim_id' = 1,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Levis',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = F,
       'description' = 'Original Levis Dissertation Simulation')

write_rds(params, 'inputs/sim_params_1.rds')  


### Simulation 2 
### Better Positivity Properties
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = 0.2,
                     'GENDER:raceeth_hispanic' = 0.1),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325,
                     'A:raceeth_hispanic' = 0.09,
                     'raceeth_hispanic' = 0.02,
                     'A:bbl' = -0.003),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =   list('A' = .1,
                         'ptwc' = -.6,
                         'A:ptwc' = -0.05)
  )

model_types <- 
  list('mu' = 'levis_paper',
       'eta' = 'levis_paper',
       'pi' = 'levis_paper',
       'lambda' = 'levis_paper',
       'lambda_shape' = 'levis_paper')

params <- 
  list('sim_id' = 2,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Levis',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = F,
       'description' = 'First Scenario that Didn\'t Violate Positivity')

write_rds(params, 'inputs/sim_params_2.rds')  

### Simulation 3
### Better Positivity Properties than original
### Still trying to induce bias under misspecification
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = 0.2,
                     'GENDER:raceeth_hispanic' = 0.4),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325,
                     'A:raceeth_hispanic' = 0.09,
                     'raceeth_hispanic' = 0.02,
                     'A:bbl' = -0.003),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =   list('A' = .1,
                         'ptwc' = -.6,
                         'A:ptwc' = -0.05)
  )

model_types <- 
  list('mu' = 'levis_paper',
       'eta' = 'levis_paper',
       'pi' = 'levis_paper',
       'lambda' = 'levis_paper',
       'lambda_shape' = 'levis_paper')

params <- 
  list('sim_id' = 3,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Levis',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = F,
       'description' = 'Don\'t Violate Positivity | Try to Induce Bias')

write_rds(params, 'inputs/sim_params_3.rds')  


### Simulation 4
### Better Positivity Properties than original
### Still trying to induce bias under misspecification
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = -0.2,
                     'GENDER:raceeth_hispanic' = 0.4),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325 + 0.1,
                     'A:raceeth_hispanic' = 0.09 + 0.05,
                     'raceeth_hispanic' = 0.02,
                     'A:bbl' = -0.003-0.002),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =   list('A' = 0.3,
                         'ptwc' = -.6,
                         'A:ptwc' = -0.5)
  )

model_types <- 
  list('mu' = 'levis_paper',
       'eta' = 'levis_paper',
       'pi' = 'levis_paper',
       'lambda' = 'levis_paper',
       'lambda_shape' = 'levis_paper')

params <- 
  list('sim_id' = 4,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Levis',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = F,
       'description' = 'Don\'t Violate Positivity | Try to Induce Bias')

write_rds(params, 'inputs/sim_params_4.rds')  


### Simulation 5
### Non-Linear Interactions in Lambda
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = 0.2,
                     'GENDER:raceeth_hispanic' = 0.4),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325,
                     'A:raceeth_hispanic' = 0.09,
                     'raceeth_hispanic' = 0.02,
                     'A:bbl' = -0.003),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =   list('A' = .1,
                         'ptwc' = -.6,
                         'A:ptwc' = -0.05,
                         'bbl' = -0.01,
                         'I(bbl^2)' = 0.0001)
  )

model_types <- 
  list('mu' = 'levis_paper',
       'eta' = 'levis_paper',
       'pi' = 'levis_paper',
       'lambda' = 'interact_nonlinear',
       'lambda_shape' = 'levis_paper')

params <- 
  list('sim_id' = 5,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Levis',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = F,
       'description' = 'Additional Interactions/Non-Linearities in "Imputation" Model')

write_rds(params, 'inputs/sim_params_5.rds')  


### Simulation 6
### Non-Linear Interactions in Lambda
### Change Lambda Shape to be 1 
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = 0.2,
                     'GENDER:raceeth_hispanic' = 0.4),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325,
                     'A:raceeth_hispanic' = 0.09,
                     'raceeth_hispanic' = 0.02,
                     'A:bbl' = -0.003),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =   list('A' = .1,
                         'ptwc' = -.6,
                         'A:ptwc' = -0.05,
                         'bbl' = -0.01,
                         'I(bbl^2)' = 0.0001)
  )

model_types <- 
  list('mu' = 'levis_paper',
       'eta' = 'levis_paper',
       'pi' = 'levis_paper',
       'lambda' = 'interact_nonlinear',
       'lambda_shape' = 1)

params <- 
  list('sim_id' = 6,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Levis',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = F,
       'description' = 'Additional Interactions/Non-Linearities in "Imputation" Model + Lower Gamma Shape')

write_rds(params, 'inputs/sim_params_6.rds')  

### Simulation 7
### Simultion 2 redux w/ multiple confounders
### Better Positivity Properties
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = 0.2,
                     'GENDER:raceeth_hispanic' = 0.1),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325,
                     'A:raceeth_hispanic' = 0.09,
                     'raceeth_hispanic' = 0.02,
                     'A:bbl' = -0.003),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =   list('A' = .1,
                         'ptwc' = -.6,
                         'A:ptwc' = -0.05)
  )

model_types <- 
  list('mu' = 'levis_paper',
       'eta' = 'levis_paper',
       'pi' = 'levis_paper',
       'lambda' = 'levis_paper',
       'lambda_shape' = 'levis_paper',
       'rho' = 'durable')

params <- 
  list('sim_id' = 7,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Levis',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = T,
       'description' = 'First Scenario that Didn\'t Violate Positivity + Multiple Missing Confounders')

write_rds(params, 'inputs/sim_params_7.rds')  

### Simulation 8
### Simulation 3 Redux w/ Multiple Confounders
### Better Positivity Properties than original
### Still trying to induce bias under misspecification
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = 0.2,
                     'GENDER:raceeth_hispanic' = 0.4),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325,
                     'A:raceeth_hispanic' = 0.09,
                     'raceeth_hispanic' = 0.02,
                     'A:bbl' = -0.003),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =   list('A' = .1,
                         'ptwc' = -.6,
                         'A:ptwc' = -0.05)
  )

model_types <- 
  list('mu' = 'levis_paper',
       'eta' = 'levis_paper',
       'pi' = 'levis_paper',
       'lambda' = 'levis_paper',
       'lambda_shape' = 'levis_paper',
       'rho' = 'durable')

params <- 
  list('sim_id' = 8,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Levis',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = T,
       'description' = 'Don\'t Violate Positivity | Try to Induce Bias | Multiple Missing Confounders')

write_rds(params, 'inputs/sim_params_8.rds')  


### Simulation 9
### Simulation 4 Redux w/ Multiple Missing Confounders
### Better Positivity Properties than original
### Still trying to induce bias under misspecification
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = -0.2,
                     'GENDER:raceeth_hispanic' = 0.4),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325 + 0.1,
                     'A:raceeth_hispanic' = 0.09 + 0.05,
                     'raceeth_hispanic' = 0.02,
                     'A:bbl' = -0.003-0.002),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =   list('A' = 0.3,
                         'ptwc' = -.6,
                         'A:ptwc' = -0.5)
  )

model_types <- 
  list('mu' = 'levis_paper',
       'eta' = 'levis_paper',
       'pi' = 'levis_paper',
       'lambda' = 'levis_paper',
       'lambda_shape' = 'levis_paper',
       'rho' = 'durable')

params <- 
  list('sim_id' = 9,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Levis',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = T,
       'description' = 'Don\'t Violate Positivity | Try to Induce Bias | Multiple Missing Confounders')

write_rds(params, 'inputs/sim_params_9.rds')  


### Simulation 10
### Simulation 5 Redux w/ LP = 2 
### Non-Linear Interactions in Lambda
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = 0.2,
                     'GENDER:raceeth_hispanic' = 0.4),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325,
                     'A:raceeth_hispanic' = 0.09,
                     'raceeth_hispanic' = 0.02,
                     'A:bbl' = -0.003),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =   list('A' = .1,
                         'ptwc' = -.6,
                         'A:ptwc' = -0.05,
                         'bbl' = -0.01,
                         'I(bbl^2)' = 0.0001)
  )

model_types <- 
  list('mu' = 'levis_paper',
       'eta' = 'levis_paper',
       'pi' = 'levis_paper',
       'lambda' = 'interact_nonlinear',
       'lambda_shape' = 'levis_paper',
       'rho' = 'durable')

params <- 
  list('sim_id' = 10,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Levis',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = T,
       'description' = 'Additional Interactions/Non-Linearities in "Imputation" Model/Multiple Missing Confounders')

write_rds(params, 'inputs/sim_params_10.rds')  


### Simulation 11
### Simulation 6 Redux w/ Multiple Missing Confounders
### Non-Linear Interactions in Lambda
### Change Lambda Shape to be 1 
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = 0.2,
                     'GENDER:raceeth_hispanic' = 0.4),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325,
                     'A:raceeth_hispanic' = 0.09,
                     'raceeth_hispanic' = 0.02,
                     'A:bbl' = -0.003),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =   list('A' = .1,
                         'ptwc' = -.6,
                         'A:ptwc' = -0.05,
                         'bbl' = -0.01,
                         'I(bbl^2)' = 0.0001)
  )

model_types <- 
  list('mu' = 'levis_paper',
       'eta' = 'levis_paper',
       'pi' = 'levis_paper',
       'lambda' = 'interact_nonlinear',
       'lambda_shape' = 1,
       'rho' = 'durable')

params <- 
  list('sim_id' = 11,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Levis',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = T,
       'description' = 'Additional Interactions/Non-Linearities in "Imputation" Model + Lower Gamma Shape + Multiple Missing Confounders')

write_rds(params, 'inputs/sim_params_11.rds')  


### Simulation 12
### Non-Linear Interactions in Lambda
### Change Lambda Shape to be pretty skewed 
### lp2 model just depends on lp1
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = 0.2,
                     'GENDER:raceeth_hispanic' = 0.4),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325,
                     'A:raceeth_hispanic' = 0.09,
                     'raceeth_hispanic' = 0.02,
                     'A:bbl' = -0.003),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =   list('A' = .1,
                         'ptwc' = -.6,
                         'A:ptwc' = -0.05,
                         'bbl' = -0.01,
                         'I(bbl^2)' = 0.0001)
  )

model_types <- 
  list('mu' = 'levis_paper',
       'eta' = 'levis_paper',
       'pi' = 'levis_paper',
       'lambda' = 'interact_nonlinear',
       'lambda_shape' = 1,
       'rho' = 'cbl_only')

params <- 
  list('sim_id' = 12,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Levis',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = T,
       'description' = 'Additional Interactions/Non-Linearities in "Imputation" Model + Lower Gamma Shape + Multiple Missing Confounders +  LP2 ~ LP1 Directly')

write_rds(params, 'inputs/sim_params_12.rds')  

### Simulation 13
### Non-Linear Interactions in Lambda
### Change Lambda Shape to be pretty skewed 
### Rho Shifts
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = 0.2,
                     'GENDER:raceeth_hispanic' = 0.4),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325,
                     'A:raceeth_hispanic' = 0.09,
                     'raceeth_hispanic' = 0.02,
                     'A:bbl' = -0.003),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =  list('A' = 0.3,
                        'ptwc' = -.6,
                        'A:ptwc' = -0.5,
                        'bbl' = -0.01,
                        'I(bbl^2)' = 0.0001),
       
       'rho' = list('ptwc' = -4)
  )

model_types <- 
  list('mu' = 'levis_paper',
       'eta' = 'levis_paper',
       'pi' = 'levis_paper',
       'lambda' = 'interact_nonlinear',
       'lambda_shape' = 1,
       'rho' = 'cbl_ptwc')

params <- 
  list('sim_id' = 13,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Levis',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = T,
       'description' = 'Additional Interactions/Non-Linearities in "Imputation" Model + Lower Gamma Shape + Multiple Missing Confounders + LP2 ~ LP1 + Y Directly')

write_rds(params, 'inputs/sim_params_13.rds')  


### Simulation 14
### Non-Linear Interactions in Lambda
### Change Lambda Shape to be pretty skewed 
### Rho Shifts to amplyify lp2 ~ y relation in more complex durable model
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = 0.2,
                     'GENDER:raceeth_hispanic' = 0.4),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325,
                     'A:raceeth_hispanic' = 0.09,
                     'raceeth_hispanic' = 0.02,
                     'A:bbl' = -0.003),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =  list('A' = 0.3,
                        'ptwc' = -.6,
                        'A:ptwc' = -0.5,
                        'bbl' = -0.01,
                        'I(bbl^2)' = 0.0001),
       
       'rho' = list('ptwc' = -4)
  )

model_types <- 
  list('mu' = 'levis_paper',
       'eta' = 'levis_paper',
       'pi' = 'levis_paper',
       'lambda' = 'interact_nonlinear',
       'lambda_shape' = 1,
       'rho' = 'durable')

params <- 
  list('sim_id' = 14,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Levis',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = T,
       'description' = 'Additional Interactions/Non-Linearities in "Imputation" Model + Lower Gamma Shape + Multiple Missing Confounders + Amplify Y')

write_rds(params, 'inputs/sim_params_14.rds')  


### Simulation 15
### Non-Linear Interactions in Lambda
### Change Lambda Shape to be pretty skewed 
### A/Y Interaction in LP2
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = 0.2,
                     'GENDER:raceeth_hispanic' = 0.4),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325,
                     'A:raceeth_hispanic' = 0.09,
                     'raceeth_hispanic' = 0.02,
                     'A:bbl' = -0.003),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =  list('A' = 0.3,
                        'ptwc' = -.6,
                        'A:ptwc' = -0.5,
                        'bbl' = -0.01,
                        'I(bbl^2)' = 0.0001)
  )

model_types <- 
  list('mu' = 'levis_paper',
       'eta' = 'levis_paper',
       'pi' = 'levis_paper',
       'lambda' = 'interact_nonlinear',
       'lambda_shape' = 1,
       'rho' = 'durable_interact')

params <- 
  list('sim_id' = 15,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Levis',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = T,
       'description' = 'Additional Interactions/Non-Linearities in "Imputation" Model + Lower Gamma Shape + A/Y/LP2 Interact')

write_rds(params, 'inputs/sim_params_15.rds')  



### Simulation 16
### Non-Linear Interactions in Lambda
### Change Lambda Shape to be pretty skewed 
### A/Y, A/LP1, LP1/Y Interaction in LP2
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = 0.2,
                     'GENDER:raceeth_hispanic' = 0.4),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325,
                     'A:raceeth_hispanic' = 0.09,
                     'raceeth_hispanic' = 0.02,
                     'A:bbl' = -0.003),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =  list('A' = 0.3,
                        'ptwc' = -.6,
                        'A:ptwc' = -0.5,
                        'bbl' = -0.01,
                        'I(bbl^2)' = 0.0001)
  )

model_types <- 
  list('mu' = 'levis_paper',
       'eta' = 'levis_paper',
       'pi' = 'levis_paper',
       'lambda' = 'interact_nonlinear',
       'lambda_shape' = 1,
       'rho' = 'aylp1_interact')

params <- 
  list('sim_id' = 16,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Levis',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = T,
       'description' = 'Additional Interactions/Non-Linearities in "Imputation" Model + Lower Gamma Shape + A/Y/LP1/LP2 Interact')

write_rds(params, 'inputs/sim_params_16.rds')  

### Simulation 17
### First Alternative Factorization Simulation
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = 0.2,
                     'GENDER:raceeth_hispanic' = 0.1),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325,
                     'A:raceeth_hispanic' = 0.09,
                     'raceeth_hispanic' = 0.02,
                     'cbl' = 0.001,
                     'A:cbl' = 0.001,
                     'A:bbl' = -0.003),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =   list('GENDER' = 0.2,
                         'GENDER:raceeth_hispanic' = 0.1)
                         
  )

model_types <- 
  list('mu' = 'alternative',
       'eta' = 'alternative',
       'pi' = 'alternative',
       'lambda' = 'alternative',
       'lambda_shape' = 'levis_paper')

params <- 
  list('sim_id' = 17,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Alternative',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = F,
       'description' = 'First Scenario')

write_rds(params, 'inputs/sim_params_17.rds')  


### Simulation 18
### Alternative Factorization Simulation
### More Complex
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = 0.2,
                     'GENDER:raceeth_hispanic' = 0.1,
                     'cbl' = -0.02,
                     'cbl:raceeth_hispanic' = -0.02),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325,
                     'A:raceeth_hispanic' = 0.09,
                     'raceeth_hispanic' = 0.02,
                     'cbl' = -0.001,
                     'A:cbl' = -0.001,
                     'A:bbl' = 0.003),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =   list('GENDER' = 0.2,
                         'GENDER:raceeth_hispanic' = 0.1,
                         'bbl:GENDER' = 0.001)
       
  )

model_types <- 
  list('mu' = 'alternative',
       'eta' = 'alternative',
       'pi' = 'alternative',
       'lambda' = 'alternative',
       'lambda_shape' = 'levis_paper')

params <- 
  list('sim_id' = 18,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Alternative',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = F,
       'description' = 'Lp Interactions in Treatment Model + Other Complexities')

write_rds(params, 'inputs/sim_params_18.rds')  


### Simulation 19
### Alternative Factorization Simulation
### More Complex
### 2 missing confounders
###
### Note have to encode Lp2 Y relationship through shifts only 
### [To-do truth in this setting]
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = 0.2,
                     'GENDER:raceeth_hispanic' = 0.1,
                     'cbl' = -0.02,
                     'cbl:raceeth_hispanic' = -0.02),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325,
                     'A:raceeth_hispanic' = 0.09,
                     'raceeth_hispanic' = 0.02,
                     'cbl' = -0.001,
                     'A:cbl' = -0.001,
                     'A:bbl' = 0.003,
                     'smoker' = 0.2,
                     'A:smoker' = 0.1,
                     'smoker:cbl' = -0.001),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =   list('GENDER' = 0.2,
                         'GENDER:raceeth_hispanic' = 0.1,
                         'bbl:GENDER' = 0.001)
       
  )

model_types <- 
  list('mu' = 'alternative',
       'eta' = 'alternative',
       'pi' = 'alternative',
       'lambda' = 'alternative',
       'rho' = 'alternative',
       'lambda_shape' = 'levis_paper')

params <- 
  list('sim_id' = 19,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Alternative',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = T,
       'description' = 'Lp Interactions in Treatment Model + Other Complexities')
write_rds(params, 'inputs/sim_params_19.rds') 

### Simulation 20
### Alternative Factorization Simulation
### More Complex
### 2 missing confounders
###
## Same as 19 but make Levis Flexible
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = 0.2,
                     'GENDER:raceeth_hispanic' = 0.1,
                     'cbl' = -0.02,
                     'cbl:raceeth_hispanic' = -0.02),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325,
                     'A:raceeth_hispanic' = 0.09,
                     'raceeth_hispanic' = 0.02,
                     'cbl' = -0.001,
                     'A:cbl' = -0.001,
                     'A:bbl' = 0.003,
                     'smoker' = 0.2,
                     'A:smoker' = 0.1,
                     'smoker:cbl' = -0.001),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =   list('GENDER' = 0.2,
                         'GENDER:raceeth_hispanic' = 0.1,
                         'bbl:GENDER' = 0.001)
       
  )

model_types <- 
  list('mu' = 'levis_flexible',
       'eta' = 'levis_flexible',
       'pi' = 'levis_flexible',
       'lambda' = 'levis_flexible',
       'rho' = 'levis_flexible',
       'lambda_shape' = 'levis_paper')

params <- 
  list('sim_id' = 20,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Alternative',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = T,
       'description' = 'Flexible Levis')
write_rds(params, 'inputs/sim_params_20.rds')  



### Simulation 21
### Alternative Factorization Simulation
### More Complex
### 2 missing confounders
###
## Same as 19/20 but make Levis Flexible
### Linear model for lambda (lp1)
shifts <- 
  list('eta' =  list('GENDER' = 0.15, 
                     'raceeth_hispanic' = 0.2,
                     'GENDER:raceeth_hispanic' = 0.1,
                     'cbl' = -0.02,
                     'cbl:raceeth_hispanic' = -0.02),
       
       'mu' =   list('GENDER:raceeth_hispanic' = -0.3,
                     'A:GENDER' = 0.325,
                     'A:raceeth_hispanic' = 0.09,
                     'raceeth_hispanic' = 0.02,
                     'cbl' = -0.001,
                     'A:cbl' = -0.001,
                     'A:bbl' = 0.003,
                     'smoker' = 0.2,
                     'A:smoker' = 0.1,
                     'smoker:cbl' = -0.001),
       
       'pi' =   list('(Intercept)' = 1,
                     'raceeth_hispanic' = -3,
                     'ptwc' = 2.7,
                     'A' = 2.4,
                     'A:ptwc' = 2.29,
                     'GENDER:ptwc' = 1.99,
                     'bbl:raceeth_hispanic' = -0.14,
                     'A:bbl' = 0.18),
       
       'lambda' =   list('GENDER' = 0.2,
                         'GENDER:raceeth_hispanic' = 0.1,
                         'bbl:GENDER' = 0.001)
       
  )

model_types <- 
  list('mu' = 'levis_flexible',
       'eta' = 'levis_flexible',
       'pi' = 'levis_flexible',
       'lambda' = 'levis_flexible_linear',
       'rho' = 'levis_flexible',
       'lambda_shape' = 'levis_paper')

params <- 
  list('sim_id' = 20,
       'n_patients' = 4344,
       'n_sims' = 5000,
       'data_generating_process' = 'Alternative',
       'shifts' = shifts,
       'models' = model_types,
       'multiple_lp' = T,
       'description' = 'Flexible Levis')
write_rds(params, 'inputs/sim_params_21.rds')  
