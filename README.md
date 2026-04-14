# Manaus COVID-19 transmission
This repository contains code and models to reproduce the analysis of COVID-19 transmission dynamics in Manaus using Bayesian inference implemented in Stan (2.33.1) using cmdstan r package (0.8.0). Model inference was done with R (4.3.1) using Rstudio (2023.06.2+561).

Analyses were published as
> van Meenen JP, van Dorp CH, de Boer RJ, van Boven M, [title]. [publisher]. [year]. [doi] 

### Installation guide
R Version 3.6.0 https://www.r-project.org/<br />
R Studio Version 1.3.959 (Interface to R) https://rstudio.com/<br />
cmdstanr R package Version 0.8.0 https://mc-stan.org/cmdstanr/<br />

### Repository structure
#### 0. Data processing steps 
Raw data are not provided in this repository. Data sources and proper citation are described in the Materials and Methods section of the associated paper. 
- `00`: Load the raw data
- `01_a`: Prepare data for Stan
- `01_b`: Convert 5-year binned contact matrix to 10-year bins
- `01_c`: Format all data for Stan, generating:
	- `job/input/data_stan_Manaus_2w_strainfreq_survey.RDS` (model input data)
	- `job/input/inits_list_2w_updated_strainfreq` (initial values)
- `01_d`: Create waning immunity prior
- `02`: Generate Stan models (stored in `job/stan/`)

#### 1. Running the model
All processed data required to run the model is found under `job/input/`. The required models are found under `job/stan/`
Run 
```
job/Manaus_stan_job_2w_updated_strainfreq_RK_NB_survey_forcedbetaifr.R
```

For the SIRS model: 
```
job/Manaus_stan_job_2w_updated_strainfreq_RK_NB_survey_forcedbetaifr_SIRS.R
```
Outputs are stored in `job/output/`

#### 2. Analysis
Our best fit and SIRS fit required for analyses are stored in `job/output`.
- `03_a`: Full model analysis (fit, parameter estimates, calculations, plots)
- `03_b`: SIRS model analysis
- `04`: Counterfactual analysis 

#### 3. Published figures
- `figures/PaperImages_updated`: Generates final figures used in the paper

### How to Cite

If you use this code or model in your work, please cite:

**[Authors]** (Year). *XXXXX*. GitHub repository.  
Available at: https://github.com/[username]/[repository-name]

You may also cite the associated publication:

**[Authors]** (Year). *[Paper title]*. *[Journal]*. https://doi.org/[DOI]