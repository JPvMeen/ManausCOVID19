options(mc.cores = parallel::detectCores())

path <- "/linuxhome/tmp/JP"

# Check if directory exists; if not, create it (including parent dirs)
if (!dir.exists(path)) {
  dir.create(path, recursive = TRUE)
  message("Directory created: ", path)
} else {
  message("Directory already exists: ", path)
}

# Set working directory
setwd(path)

### Packages
data_management_packages_to_install <- c("data.table", "tidyverse", "readxl", "yarrr", "ggplot2", "zoo", "cmdstanr")

packages_to_install <- c(data_management_packages_to_install)

# Check and install packages using a for loop
for (package in packages_to_install) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
  library(package, character.only = TRUE)
}

### file locations
job_location <- "job"

### Loading in the data
data_stan_Manaus_2w_updated_strainfreq <- readRDS(paste0(job_location,"/input/data_stan_Manaus_2w_strainfreq_survey.RDS"))

inits_list_2w_updated_strainfreq <- readRDS(paste0(job_location,"/input/inits_list_2w_strainfreq_survey"))
inits_list_2w_updated_strainfreq

manaus_model_2w_updated_strainfreq <- cmdstan_model(paste0(job_location,"/stan/manaus_2w_updated_s1re_strainfreq_RK_NB_survey.stan"))

### Fitting procedure
fit_manaus <- manaus_model_2w_updated_strainfreq$sample(data = data_stan_Manaus_2w_updated_strainfreq,
                                          init = inits_list_2w_updated_strainfreq,
                                          refresh = 50,
                                          chains = 4,
                                          parallel_chains = 4,
                                          iter_warmup=1000,
                                          iter_sampling=1500-1000,
                                          seed = 0,
                                          #control = list(adapt_delta = 0.9)
)

current_time <- format(Sys.time(), "%Y_%m_%d_t%H_%M_%S")

filename_df <- paste0(job_location,"/output/",current_time,"_manaus_2w_updated_s1re_strainfreq_RK_NB_survey")

print(filename_df)

### Generating output
manaus_draws_df_2w_updated_strainfreq <- fit_manaus$draws(format = "df")
saveRDS(manaus_draws_df_2w_updated_strainfreq, file = filename_df)