## Inferring wildlife poaching in Southeast Asia with multispecies dynamic occupancy models and hidden Markov models. 
### Lucile Marescot and Olivier Gimenez, December 31, 2018.

### 1. Introduction. 

This repo contains data and codes associated with a submitted paper on inferring wildlife poaching in Southeast Asia with multispecies dynamic occupancy models and hidden Markov models. Specifically, we present the two-species dynamic occupancy model introduced by MacKenzie et al. (2006) (see Table 8.7 page 245) which we formulate as a hidden Markov model (HMM) and implement in R. 

Note that the data (in particular the shapefiles for spatial analyses) are in zip files that you will need to unzip to reproduce the results. 

### 2. Case study of poacher and wildlife interactions in Cambodia

The code "Smart data occupancy analysis.Rmd" loads the SMART data for the two study areas Serepok Wildlife Sancturary (SWS) and Phnom Prich Wildlife Sanctuary (PPWS) and builds the occupancy table of co-occurence data on the presence of wildlife (visual observation, tracks, scats) and poachers (snares). It is structured as follows:

* Loading the smart data, the shapefiles and the libraries 
* Data description and selection
   * focus on animals mainly ungulates
   * focus on illegal activities mainly snares
* Build and map the occupancy grid (Figure 1 is displayed)
* A dynamic two-species patch-occupancy model, formulated as a HMM 
   * the model: vector of initial occupancy, transition matrix and observation matrix
   * load the function showing all possible effects on the parameters and returning the associated vector of regression coefficients
   * load the likelihood function of each of these models 
   * selecting the best models regarding effects on detection, on occupancy and on the transition
   * estimating parameters
  

Eventually, we present the results through several figures: 

  * Plot the mean estimates and their associated confidence intervals obtained from the best model (Figure 2) 
  * Plot the detection probabilities as a function of patrolling effort (Figure 4)
  * Plot the state probabilities during the 4 years of the study period (Figure 5)
  * Spatial variation in detection and patrolling effort in the two study areas


### 3. Simulation study 

The aim of this study is to check for potential bias in the parameter estimates obtained from the best model introduced in the previous section. 

The codes "Simulation_Bias_OccupancyDesign.R" and "Simulation_Bias_Parameters.R" focus respectively on the occupancy design and on the parameters related to the detectability and occupancy of the two species. 

The R scripts to perform the simulation are structured as follows:

* Generate the simulated data sets 
   * for different study design "Simulation_Bias_OccupancyDesign.R" which calls source('Multispecies_occupancy_dynamic_simul_fitting.R')
   * for species associated to different detection and occupancy probabilities "Simulation_Bias_Parameters.R" which calls source('Multispecies_occupancy_dynamic_simul_fitting.R')
* Model fitting 
   * calculate the probability of obtainting the simulated data given the two-species dynamic occupancy model via source('Multispecies_occupancy_dynamic_simul_fitting.R')
* Calculate the mean bias and minimum square error 

This code was developped in R version 3.5.0 (2018-04-23). This code is licensed under the terms of the GNU General Public License v3.0.
