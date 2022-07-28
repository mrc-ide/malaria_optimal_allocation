# Run optimisation -------------------------------------------------------------
# Code was run on HPC

# Setup ------------------------------------------------------------------------
library(tidyverse)
library(akima)
library(GenSA)
library(here)

datadir <- here("files", "/")
outdir <- here("output", "/")

# Functions --------------------------------------------------------------------

SAfun_pf_pv <- function(data, eirs, species = "pf_pv",   
                        ITNuse,    # max ITN use, 80 or 90
                        usedist, single_budget_proportion){
  
  # single_budget_proportion is a fraction of the max_budget
  
  file <- data %>%
    semi_join(eirs,by = c("init_EIR_pf", "init_EIR_pv"))
  
  # set max ITN coverage
  file <- file %>% filter(itn_cov <= ITNuse/100)
  
  if (ITNuse == 80) {
    netz <- readRDS(paste0(datadir, "netz_median_use_rate.rds"))
    conversion <- readRDS(paste0(datadir, "conversion_usage_pcnets_median_use_rate.rds")) %>%
      as.data.frame() %>%
      filter(usage<=ITNuse/100)
    conversion$pc_nets_annual <- round(conversion$pc_nets_annual,3)
  }
  
  if (ITNuse == 90) {
    netz <- readRDS(paste0(datadir, "netz_max_use_rate.rds")) 
    conversion <- readRDS(paste0(datadir, "conversion_usage_pcnets_max_use_rate.rds")) %>% 
      as.data.frame() %>% 
      filter(usage<=ITNuse/100)
    conversion$pc_nets_annual <- round(conversion$pc_nets_annual,3)
  }
  
  file <- file %>% 
    left_join(netz, by=c('itn_cov'='usage')) %>%
    mutate(cases = cases_pf+cases_pv,
           pop = pmax(pop_pf, pop_pv))
  
  # Find maximum budget to spend to reach minimum cases
  if(usedist=='linear'){
    max_budget <- file %>%
      ungroup() %>%
      mutate(threshold = ifelse(cases >= 1,1,0)) %>%  # elimination of all malaria 
      group_by(init_EIR_pf, init_EIR_pv, threshold) %>%                                
      mutate(min_itn=ifelse(threshold == 0,min(itn_cov),max(file$itn_cov))) %>% # itn_cov at which all malaria is eliminated         
      ungroup() %>%
      group_by(init_EIR_pf, init_EIR_pv) %>%
      filter(itn_cov == ifelse(min_itn<0.8,min(min_itn), 0.8)) %>%
      mutate(cost=itn_cov*pop) %>%
      ungroup() %>%
      summarise(sum(cost)) %>% as.numeric() %>% ceiling()
  }
  
  if(usedist=='loess'){
    max_budget <- file %>%
      ungroup() %>%
      mutate(threshold = ifelse(cases >= 1,1,0)) %>%  # elimination of all malaria 
      group_by(init_EIR_pf, init_EIR_pv, threshold) %>%                                
      mutate(min_itn=ifelse(threshold == 0,min(itn_cov),max(file$itn_cov))) %>% # itn_cov at which all malaria is eliminated         
      ungroup() %>%
      group_by(init_EIR_pf, init_EIR_pv) %>%
      filter(itn_cov == ifelse(min_itn<0.8,min(min_itn), 0.8)) %>%
      mutate(cost = annual_percapita_nets_distributed*pop) %>%
      ungroup() %>%
      summarise(sum(cost)) %>% as.numeric() %>% ceiling()
    
  }
  
  # Create 2 surfaces for Pf and Pv separately
  # Since EIRs have not been summed, we can work with the incidence rate instead of cases
  
  x <- file$itn_cov
  y <- file$init_EIR_pf
  z <- file$clinical_incidence_pf
  
  dimx <- 9000
  dimy <- 9000
  
  spline_pf <- akima::interp(x,y,z, yo=seq(min(y), max(y),length.out=dimy), # define spline between a min and max value
                             xo=seq(0,1,length.out=dimx), duplicate = 'mean', linear = TRUE)
  
  x2 <- file$itn_cov
  y2 <- file$init_EIR_pv
  z2 <- file$clinical_incidence_pv
  
  spline_pv <- akima::interp(x2,y2,z2, yo=seq(min(y2), max(y2),length.out=dimy), # define spline between a min and max value
                             xo=seq(0,1,length.out=dimx), duplicate = 'mean', linear = TRUE)
  
  # Find population for each combination of EIRs
  pop <- unique(file$pop)
  
  # Check that this is in the correct order:
  if(all.equal(select(file, init_EIR_pf,init_EIR_pv,pop) %>%
               distinct() %>% as.data.frame(), 
               cbind(eirs, pop))==TRUE) {
    print("Check complete")
  } else {
    stop("EIRs and population size not matching")
  }
  
  # Objective function 
  
  if (usedist == 'linear') {  
    
    # assign constraints
    upper <- pop*ITNuse/100 %>% as.vector
    lower <- rep(0,nrow(eirs)) %>% as.vector
    
    fn <- function(spending, pop, eirs_pf, eirs_pv, budget){
      
      r <- 10000000000 # set a high penalty when the constraint is not respected
      ITN_percentage <- ifelse(spending/pop > ITNuse/100, ITNuse/100, spending/pop)
      
      spline_results <- spline_pf[["z"]][round(ITN_percentage*(dimx-1)+1),
                                         round(((eirs_pf - min(eirs_pf))/(max(eirs_pf) - min(eirs_pf)))*(dimy-1)+1)] +
        spline_pv[["z"]][round(ITN_percentage*(dimx-1)+1),
                         round(((eirs_pv - min(eirs_pv))/(max(eirs_pv) - min(eirs_pv)))*(dimy-1)+1)]
      spline_results <- diag(spline_results) * pop
      
      return(sum(spline_results, na.rm='always')  + r * max(sum(spending) - budget, 0))
    }
    
  }
  
  if (usedist == 'loess') { 
    
    # assign constraints
    upper <- pop * max(conversion[conversion$usage==ITNuse/100,]$pc_nets_annual) %>% as.vector
    lower <- rep(0,nrow(eirs)) %>% as.vector
    
    fn <- function(spending, pop, eirs_pf, eirs_pv, budget){
      
      r <- 10000000000 # set a high penalty when the constraint is not respected
      ITN_percentage1 <- spending/pop
      # Convert proportion to nets distributed
      ITN_percentage1 <- sapply(ITN_percentage1, function(x){
        ifelse(x > max(conversion$pc_nets_annual), max(conversion$pc_nets_annual), x)
      })
      ITN_percentage1 <- round(ITN_percentage1,3)
      # Convert back to matching usage to relate back to spline
      ITN_percentage <- sapply(ITN_percentage1, function(x){
        conversion[conversion$pc_nets_annual==x,]$usage
      })
      
      spline_results <- spline_pf[["z"]][round(ITN_percentage*(dimx-1)+1),
                                         round(((eirs_pf - min(eirs_pf))/(max(eirs_pf) - min(eirs_pf)))*(dimy-1)+1)] +
        spline_pv[["z"]][round(ITN_percentage*(dimx-1)+1),
                         round(((eirs_pv - min(eirs_pv))/(max(eirs_pv) - min(eirs_pv)))*(dimy-1)+1)]
      
      spline_results <- diag(spline_results) * pop
      return(sum(spline_results, na.rm='always')  + r * max(sum(spending) - budget, 0))
      
    }
  }
  
  # read in manual strategies to set starting values below to speed up optimization
  file2 <- readRDS(paste0(outdir, "bednets_data_equilibrium_",species,'_',ITNuse,'_',usedist,'.rds')) %>% 
    filter(scenario == 'high burden to high impact')
  
  # set function to run SA process
  solver <- function(budget) {
    
    par <- file2 %>% group_by(init_EIR) %>% filter(tdollar < budget) %>% 
      filter(tdollar == max(tdollar)) %>% 
      filter(itn_cov == max(itn_cov)) %>%
      ungroup() %>% 
      arrange(init_EIR) %>% mutate(dollar = usedist*pop) %>% select(init_EIR, dollar,pop) 
    
    sum_eirs <- eirs %>%
      mutate(init_EIR = round(init_EIR_pf+init_EIR_pv,6)) 
    freq <- data.frame(init_EIR = sort(unique(sum_eirs$init_EIR)),
                       freq = as.numeric(table(sum_eirs$init_EIR)))
    par <- sum_eirs %>%
      full_join(freq, by = c("init_EIR")) %>%
      left_join(par, by=c("init_EIR")) %>%
      mutate(dollar = ifelse(is.na(dollar),0,dollar/freq)) %>%  # divide by freq because EIRs are not summed here
      select(dollar) %>% unlist %>% as.numeric()
    
    result <- GenSA(par = par,
                    fn = fn,  
                    lower = lower, 
                    upper = upper, 
                    pop = pop, 
                    eirs_pf = eirs$init_EIR_pf, 
                    eirs_pv = eirs$init_EIR_pv, 
                    control = list(maxit = 5e6, temperature = 1e6),
                    budget = budget)
    
    data.frame(b=result$par, eirs=eirs, B=budget)
    
  }
  
  # set value for each step increase in $
  
  values <- floor(single_budget_proportion * max_budget)
  
  set.seed(123)
  opt_data <- map_dfr(.x=values, .f=solver, .id=NULL) # iterate
  write.csv(opt_data, paste0(datadir, "optim_output/SA_", nrow(eirs),'_', ITNuse, '_', usedist, '_', values, ".csv"))
  
  
}

pprocessloess_pf_pv <- function(file, n=63, ITNuse, usedist, PYRcost) {
  # Read in output
  files <- list.files(path = paste0(datadir, "optim_output/"), # Identify all csv files in folder
                      pattern = paste0("SA_", n, '_', ITNuse, '_', usedist, '_*'), full.names = TRUE)
  
  dat_list <- lapply(files, function (x) read_csv(x))
  dat <- bind_rows(dat_list) %>% # Combine data sets into one data set 
    rename(init_EIR_pf = eirs.init_EIR_pf,
           init_EIR_pv = eirs.init_EIR_pv) %>%
    select(-1)
  
  dmm_output_foi_bednets_wpop <- readRDS(paste0(datadir, file, '.rds')) %>% filter(itn_cov <= ITNuse/100)
  pop <- dmm_output_foi_bednets_wpop %>% select(init_EIR_pf, init_EIR_pv, pop_pf, pop_pv) %>% distinct()
  cinc <- dmm_output_foi_bednets_wpop  %>% select(init_EIR_pf, init_EIR_pv, 
                                                  itn_cov, clinical_incidence_pf,
                                                  clinical_incidence_pv)
  
  if (ITNuse/100 == .80) {
    conversion <- readRDS(paste0(datadir, "conversion_usage_pcnets_median_use_rate.rds")) %>% 
      as.data.frame() %>%
      filter(usage<=ITNuse/100)
    conversion$pc_nets_annual <- round(conversion$pc_nets_annual,3)
  }
  
  if (ITNuse/100 == .90) {
    conversion <- readRDS(paste0(datadir, "conversion_usage_pcnets_max_use_rate.rds")) %>% 
      as.data.frame() %>% 
      filter(usage<=ITNuse/100)
    conversion$pc_nets_annual <- round(conversion$pc_nets_annual,3)
  }
  
  # Elimination in the combined approach:
  # threshold remains the same - if total cases <1, malaria is eliminated
  # but the population at risk may vary depending on whether only P. falciparum or
  # P. vivax are eliminated
  
  if(usedist == "loess") {
    
    output <- dat %>% 
      group_by(B) %>% 
      mutate(b_sum=sum(b, na.rm=T)) %>%
      ungroup() %>%
      left_join(pop, by=c("init_EIR_pf", "init_EIR_pv")) %>%
      mutate(pop = pmax(pop_pf, pop_pv),
             dollar = round((b)/(pop), 3)) %>%
      left_join(conversion, by=c('dollar'='pc_nets_annual')) %>%
      mutate(itn_cov=round(usage,2)) %>%
      left_join(cinc, by=c("init_EIR_pf", "init_EIR_pv", "itn_cov")) %>%
      mutate(cases_pf = clinical_incidence_pf * pop_pf,
             cases_pv = clinical_incidence_pv * pop_pv,
             cases = cases_pf+cases_pv,
             popatrisk_pf = ifelse(cases_pf<1,0, pop_pf),
             popatrisk_pv = ifelse(cases_pv<1,0, pop_pv),
             popatrisk = pmax(popatrisk_pf, popatrisk_pv),
             B = B * PYRcost,
             b = b * PYRcost,
             b_sum = b_sum * PYRcost)
    # pop is for calculation of cost and popatrisk the resulting population at risk 
    # for a given budget being invested
    
  } else if (usedist == "linear") {
    
    output <- dat %>% 
      group_by(B) %>% 
      mutate(b_sum=sum(b, na.rm=T)) %>%
      ungroup() %>%
      left_join(pop, by=c("init_EIR_pf", "init_EIR_pv")) %>%
      mutate(pop = pmax(pop_pf, pop_pv),
             itn_cov = round(b/(pop), 2)) %>%
      left_join(cinc, by=c("init_EIR_pf", "init_EIR_pv", "itn_cov")) %>%
      mutate(cases_pf = clinical_incidence_pf * pop_pf,
             cases_pv = clinical_incidence_pv * pop_pv,
             cases = cases_pf+cases_pv,
             popatrisk_pf = ifelse(cases_pf<1,0, pop_pf),
             popatrisk_pv = ifelse(cases_pv<1,0, pop_pv),
             popatrisk = pmax(popatrisk_pf, popatrisk_pv),
             B = B * PYRcost,
             b = b * PYRcost,
             b_sum = b_sum * PYRcost)
  }
  
  # Note all outputs are summarised on the global level:
  # the population at risk of malaria represents the maximum
  # of the global sum of the population at risk of Pf and the  
  # global sum of the population at risk of Pv (NOT the sum of the global 
  # population at risk of any malaria across all settings)

}


# Run optimiser ----------------------------------------------------------------

data <- readRDS(paste0(datadir, 'dmm_output_foi_bednets_wpop_pf_pv.rds'))
eirs <- select(data, init_EIR_pf, init_EIR_pv) %>% distinct()
budgets <- c(seq(0.05,1.1,0.05))    # a proportion of max_budget is used in the function

out <- SAfun_pf_pv(data=data, eirs=eirs,
                   species='pf_pv', ITNuse=80,
                   usedist='loess', single_budget_proportion = budgets[10])

# Max. 80% coverage loess
for (i in 1:length(budgets)) {
  out <- SAfun_pf_pv(data=data, eirs=eirs,
                          species='pf_pv', ITNuse=80,
                          usedist='loess', single_budget_proportion = budgets[i])
}

# Max. 80% coverage linear
for (i in 1:length(budgets)) {
  out <- SAfun_pf_pv(data=data, eirs=eirs,
                          species='pf_pv', ITNuse=80,
                          usedist='linear', single_budget_proportion = budgets[i])
}

# Post-processing --------------------------------------------------------------
output_loess <- pprocessloess_pf_pv(file="dmm_output_foi_bednets_wpop_pf_pv",
                    ITNuse=80, usedist="loess", PYRcost = 3.5)

output_linear <- pprocessloess_pf_pv(file="dmm_output_foi_bednets_wpop_pf_pv",
                    ITNuse=80, usedist="linear", PYRcost = 3.5)