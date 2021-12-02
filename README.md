# The seasonality of cholera in sub-Saharan Africa

This is a repository containing minimal datasets, code, and results for the manuscript:

Perez-Saez, J., Lessler, J., Lee, E.C., Luquero, F.J., Malembaka, E.B., Finger, F., Langa, J.P., Yennan, S., Zaitchik, B. and Azman, A.S., 2021. *The seasonality of cholera in sub-Saharan Africa.* medRxiv. https://doi.org/10.1101/2021.11.23.21265642.

## Content
The repository is structured as follows in its current state:

```
1 .                                                                                        
2  ¦--analysis                                      # Code for analysis                    
3  ¦   ¦--make_output_metadata.R                    ## Script to make output metadata      
4  ¦   °--make_repo_tree.R                          ## Script to make the repo tree        
5  °--generated_data                                # All analysis outputs                 
6      ¦--best_seasonality_estimates_schema.json    ## Final seasonality estimates metadata
7      °--best_seasonality_estimates.csv            ## Final seasonality estimates
```

## Upcoming

- Codes for models in Stan
- Minimal datasets for testing the models
- Code for producing final figures
