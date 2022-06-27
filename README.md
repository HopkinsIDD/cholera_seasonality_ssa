# The seasonality of cholera in sub-Saharan Africa

This is a repository containing minimal datasets, code, and results for the manuscript:


Perez-Saez, J., Lessler, J., Lee, E.C., Luquero, F.J., Malembaka, E.B., Finger, F., Langa, J.P., Yennan, S., Zaitchik, B. and Azman, A.S., 2022. *The seasonality of cholera in sub-Saharan Africa: a statistical modelling study*, The Lancet Global Health. https://doi.org/10.1016/S2214-109X(22)00007-9.

**-- THIS IS THE BRANCH WITH ALL DATA FOR REPLICATION--**

### Content
The repository is structured as follows in its current state in the main branch:

```
1  .                                                                                                             
2   ¦--analysis                                                                                                  
3   ¦   ¦--01_data_preparation                           -- (not on this repo).                                  
4   ¦   ¦--02_modeling                                                                                           
5   ¦   ¦   ¦--01_run_model_prior.R                      --- Runs the seasonality model with priors only         
6   ¦   ¦   ¦--02_run_model.R                            --- Runs the seasonality model                          
7   ¦   ¦   ¦--03_generate_simulations.R                 --- Generates simulated outputs from model fits         
8   ¦   ¦   °--configs                                   --- Configs used for model running                      
9   ¦   ¦       ¦--cases                                 ---- Cases as thershold for cholera presence            
10  ¦   ¦       ¦--default_config.yml                    ---- Annotated example config                           
11  ¦   ¦       ¦--mean_annual_incidence                 ---- Mean annual incidence as threshold                 
12  ¦   ¦       °--occurrence                            ---- Occurrence as thershold                            
13  ¦   ¦--03_postprocessing                                                                                     
14  ¦   ¦   ¦--01_extract_prior_results_config.R         --- Extracts outputs from prior model fits              
15  ¦   ¦   ¦--02_extract_fitting_results_config.R       --- Extracts outputs from model fits                    
16  ¦   ¦   ¦--03_compute_seasonality_index.R            --- Computes the seasonality index                      
17  ¦   ¦   °--04_postprocess_outputs.R                  --- Postprocesses outputs for reporting                 
18  ¦   ¦--04_clustering                                                                                         
19  ¦   ¦   ¦--01_compute_grouping.R                     --- Runs clustering model                               
20  ¦   ¦   °--02_extract_clustering_results.R           --- Extracts clustering outputs                         
21  ¦   ¦--05_outputs                                                                                            
22  ¦   ¦   ¦--01_final_stats.R                          --- Makes statistics statements used in manuscript      
23  ¦   ¦   °--02_final_plots.R                          --- Makes all figures in manuscript                     
24  ¦   ¦--distribution_utils.R                          -- Helper functions for preparing this repo             
25  ¦   ¦--make_output_metadata.R                        -- Script to make output metadata                       
26  ¦   ¦--make_repo_tree.R                              -- Script to make the repo tree                         
27  ¦   ¦--modeling_utils.R                              -- Helper functions for seasonality models              
28  ¦   ¦--nbgraph_utils.R                               -- Helper functions for spatial adjacency               
29  ¦   ¦--plotting_utils.R                              -- Helper functions for plotting                        
30  ¦   ¦--postprocessing_utils.R                        -- Helper functions for postprocessing                  
31  ¦   ¦--simulation_utils.R                            -- Helper functions for simulations                     
32  ¦   ¦--stan                                          -- All stan code                                        
33  ¦   °--utils.R                                       -- Misc helper functions                                
34  ¦--figures                                                                                                   
35  ¦   °--pub_figures                                   -- Figures in main text and supplement                  
36  ¦--generated_data                                                                                            
37  ¦   ¦--best_seasonality_estimates_schema.json        -- Final seasonality estimates metadata                 
38  ¦   ¦--best_seasonality_estimates.csv                -- Final seasonality estimates                          
39  ¦   ¦--data_mapping                                  -- Formatted data for stan                              
40  ¦   ¦--geodata                                       -- Geographic data used in the analysis                 
41  ¦   ¦--model_outputs                                 -- Modeling outputs used to reproduce figures and tables
42  ¦   ¦--postprocessing                                -- Postprocessed data to reproduce figures and tables   
43  ¦   °--run_data                                      -- Data used to run the models                          
44  ¦       ¦--all_run_data_mean_annual_incidence.csv    --- Compilation of all run data for the main analysis   
45  ¦       °--all_run_data_schema.json                  --- Metadata of the run data                            
46  °--manuscript                                                                                                
47      °--supplement                                                                                            
48          ¦--refs.bib                                  --- References used in the supplement                   
49          ¦--supplement_cholera_seasonality.pdf        --- Supplementary material pdf                          
50          °--supplement_cholera_seasonality.Rmd        --- RMD to produce supplementary material pdf   
```

The content of the main branch consists of:

- all code used in the analysis
- raw figures found in the publication (`figures/pub_figures`)
- the analysis dataset (`generated_data/run_data/all_run_data_mean_annual_incidence.csv`)
- main outputs (seasonality estimats in `generated_data/best_seasonality_estiamtes.csv`)
- the code and pdf of the supplementary material (`manuscript/supplement`).

### Details on project implementation

Details on project implementation are given in the notebook in `notebooks/project_details.html`. This file is currently empty and will be updated soon.

### Replicate the results
All data and outputs necessary to replicate the results are available in this branch. Some files are large and require git lfs to download (https://git-lfs.github.com/). These can be pulled running:

```
git lfs pull
```