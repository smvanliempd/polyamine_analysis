library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# global vars/functions
mets <- c("Putrescine","Spermidine","Spermine")
qtl_vec <- c(0.1,0.25,0.5,0.75,.90)
seed <- 2001
calc_conc <- function(y, a, b)  (y/exp(a))^(1/b)
get_conc <- function(dt ,dt_pars = d.qnt.pars, signal_type,
                     merge_cols  =  c("Assay.Date", "Analyte")) {
  cols <- c("b_low","a_low","b_high","a_high","signal_cut",merge_cols)
  dt <- merge(dt, 
              dt_pars[,..cols],
              by = merge_cols, 
              all.x = TRUE)
  dt[ , Conc_raw_alt :=  ifelse(Signal < signal_cut, 
                                calc_conc(get(signal_type),a_low,b_low), 
                                calc_conc(get(signal_type),a_high,b_high))]
  return(dt)
}

# create dir to store graph and table output
dir.create(paste0(getwd(), "/graphs_tables"), showWarnings = FALSE)

# run LC analysis
source("chromatograms.R")

# run assay performance analysis
# get data
d.all <- read.csv(file = paste0(getwd(),"/data/assay_data_v2.csv"),sep = ",") |> 
  data.table()
d.all[d.all==""] <- NA

# run analyses
source("extraction_efficiency.R")
source("recovery.R")
source("quantification.R")
source("assay_variability.R")
source("analyte_attenuation.R")
source("combine_data.R")
source("use_case_meds.R")
source("use_case_inhibs.R")
