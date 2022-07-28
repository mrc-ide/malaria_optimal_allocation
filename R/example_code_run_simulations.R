# Run P. falciparum and P. vivax simulations -----------------------------------
# Code was run on HPC

# Setup ------------------------------------------------------------------------
library(tidyverse)
library(ICDMM)
library(here)

datadir <- here("files", "/")

# Functions --------------------------------------------------------------------

get_model_values_bednets_yearly_pf <- function(eir, itn_cov, num_int) {
  
  year <- 365
  age_vector <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80)
  
  out <-  run_model(model = 'odin_model',
                    age = age_vector,
                    time = 75*year,    # run for 75 years (to equilibrium)
                    het_brackets = 5,
                    init_ft = 0.4, 
                    init_EIR = eir,
                    itn_cov = itn_cov, 
                    num_int = num_int, 
                    ITN_on = 1, 
                    ITN_interval = year*1) %>%  
    as_tibble(rownames = NA) 
  
  m <- out %>% 
    mutate(years = ceiling(t/365)) %>% 
    group_by(years) %>%
    summarize(clinical_incidence = mean(inc) * year, 
              clinical_incidence_05 = mean(inc05) * year,
              init_EIR = eir,
              itn_cov = itn_cov)
  
  write.csv(m, paste0(datadir, "model_runs_pf/foi","_",eir,"_",itn_cov,".csv"))
  
} 

get_model_values_bednets_yearly_pv <- function(eir, itn_cov, num_int){
  
  year <- 365
  age_vector <- c(0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,3.5,5,7.5,10,15,20,30,40,50,60,70,80)
  maxyear <- 175
  
  itn_cov <- ifelse(itn_cov==1, 0.9999, itn_cov)
  
  out <-  run_model(model = 'vivax_model',
                    age = age_vector,
                    time = maxyear*year,       # run for 175 years (to equilibrium)
                    het_brackets = 5,
                    K_max=2,
                    init_ft = 0.4, 
                    init_EIR = eir,
                    itn_cov = itn_cov, 
                    irs_cov = 0, 
                    num_int = num_int, 
                    ITN_IRS_on = 1, 
                    ITN_interval = year*1,  
                    stop_rebound_switch = 1,   # enable option to prevent late rebounds due to relapse
                    hypnozoite_prev_threshold = 1e-5) %>%   # threshold at which stop_rebound_switch is triggered
    as_tibble(rownames = NA) 
  
  m <- out %>% 
    mutate(years = ceiling(t/365)) %>% 
    group_by(years) %>%
    summarize(prevLM = mean(prevLM, na.rm=T), # mean prev
              clinical_incidence = mean(inc) * year, # DMM outputs daily CI
              clinical_incidence_05 = mean(inc05) * year,
              init_EIR = eir,
              itn_cov = itn_cov)
  
  m$itn_cov[m$itn_cov==0.9999] <- 1

  write.csv(m, paste0(datadir, "model_runs_pv/foi","_",eir,"_",itn_cov,".csv"))
  
} 

# Run P. falciparum simulations ------------------------------------------------
eir <- c(0.001,0.01,0.05,seq(0.1,0.9,0.1),seq(1,10,1),seq(15,140,5)) # EIRs for model 
itn_cov <- seq(0,1,0.01)
year <- 365

# creating every combination of ITN coverage % and EIR
combo <- crossing(eir,itn_cov) %>% 
  as_tibble() %>%
  mutate(num_int = ifelse(itn_cov==0,1,2)) %>%
  dplyr::select(eir,itn_cov,num_int)

# run model for each combo
out <- apply(combo, 1, function(x) get_model_values_bednets_yearly_pf(eir = as.list(x)$eir,
                                                                          itn_cov = as.list(x)$itn_cov,
                                                                          num_int= as.list(x)$num_int))

# read in outputs
data_all <- list.files(path = paste0(datadir, "model_runs_pf"), # Identify all csv files in folder
                       pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>% # Store all files in list
  bind_rows # Combine data sets into one data set 
head(data_all)  
# Model returns case incidence detected using daily ACD

# Run P. vivax simulations ----------------------------------------------------

# This requires installing the "feat/vivax" branch of ICDMM

# Pre-selected EIRs for model based on match with current prevalence:
eir <- c(0.001, 0.002, 0.003, 0.004, 0.006, 0.010, 
         0.020, 0.030, 0.040, 0.060, 0.070, 0.080, 
         0.090, 0.100, 0.200, 0.300, 0.500, 1.400, 1.600, 1.700)

itn_cov <- seq(0,1,0.01)

# creating every combination of % and EIR
combo <- crossing(eir,itn_cov) %>% as_tibble() %>%
  mutate(num_int = ifelse(itn_cov==0,1,2)) %>%
  dplyr::select(eir,itn_cov,num_int) %>%
  mutate(itn_cov = ifelse(itn_cov==1, 0.9999, itn_cov))

# run model for each combo
out <- apply(combo, 1, function(x) 
  get_model_values_bednets_yearly_pv(eir = as.list(x)$eir,
                                     itn_cov = as.list(x)$itn_cov,
                                     num_int= as.list(x)$num_int))

# read in outputs
data_all <- list.files(path = paste0(datadir, "model_runs_pv"), # Identify all csv files in folder
                       pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>% # Store all files in list
  bind_rows # Combine data sets into one data set 
head(data_all)  
# Model returns case incidence detected using daily ACD
