## Setup ------------------------------------------------------------------------

# packages
library(tidyverse)
library(fuzzyjoin)
library(here)

# working directory
datadir <- here("files", "/")
outdir <- here("output", "/")

# set threshold to define elimination
elimination_threshold <- 1 # 1 case

# set cost per bednet 
PYRcost <- 3.5

## Prepare data: group Pf+Pv dataset by sum of EIRs ----------------------------

# Load files with population size incorporated
sims_pf_pv <- readRDS(file.path(datadir, "dmm_output_foi_bednets_wpop_pf_pv.rds"))

# Assume EIR is the sum of Pf and Pv EIRs
# Cases always need to be calculated with the respective population at risk (not total)
options(pillar.sigfigs=100)
sims_pf_pv <- sims_pf_pv  %>% 
  mutate(init_EIR = round(init_EIR_pf+init_EIR_pv,6))

# Checks
sum(filter(sims_pf_pv, itn_cov==0)$n)==105    # test # of countries
# total cases
sum(filter(sims_pf_pv, itn_cov==0)$cases_pf+filter(sims_pf_pv, itn_cov==0)$cases_pv)
# 321327854
# global population at risk
pmax(sum(filter(sims_pf_pv, itn_cov==0)$pop_pf),sum(filter(sims_pf_pv, itn_cov==0)$pop_pv))
# 4147518844

# Combine those that have the same sum of EIRs
sims_pf_pv <- select(sims_pf_pv, -init_EIR_pf, -init_EIR_pv) %>%
  group_by(init_EIR, itn_cov) %>%
  summarise(pop_pf=sum(pop_pf),    
            pop_pv=sum(pop_pv),
            cases_pf = sum(cases_pf),
            cases_pv = sum(cases_pv),
            n = sum(n)) %>%
  mutate(pop = pmax(pop_pf, pop_pv),
         cases = cases_pf+cases_pv)
# Need to drop clinical_incidence here and work with previously calculated cases_pf and _pv

# Checks
sum(filter(sims_pf_pv, itn_cov==0)$n)==105    # test # of countries
all.equal(sum(filter(sims_pf_pv, itn_cov==0)$cases_pf+filter(sims_pf_pv, itn_cov==0)$cases_pv),
          321327854)
all.equal(pmax(sum(filter(sims_pf_pv, itn_cov==0)$pop_pf),sum(filter(sims_pf_pv, itn_cov==0)$pop_pv)),
          4147518844)
length(unique(sims_pf_pv$init_EIR)) # 55 unique EIRs

## Calculate weights for proportional allocation -------------------------------

# Weights for proportional allocation are based on Global Fund burden metric
# scroll down to 'all annexes in excel format': https://www.who.int/teams/global-malaria-programme/reports/world-malaria-report-2020
calculate_weights_gf <- function(match_dataset, species) {
  # Prepare dataset with country-specific burden metrics
  who_data <- read.csv(file.path(datadir, "WMR_2020_annex3F.csv")) %>% 
    select(ISO, country, year, population, cases, deaths) %>%
    mutate(cases = as.numeric(cases),
           deaths = as.numeric(deaths),
           incidence_rate= cases/population,
           mortality_rate = cases/population) 
  
  par_latest <- group_by(who_data, ISO) %>% filter(year==max(year)) %>%
    mutate(par_latest=population) %>%
    select(ISO, par_latest)
  
  normalize <- function(x, na.rm = TRUE) {
    return((x- min(x)) /(max(x)-min(x)))
  }
  
  # Burden metric is based on formula in Global Fund document (Description of the 2020-2022 Allocation Methodology)
  # Based on mean values from 2000-2004 for each country
  country_burden_metrics <- filter(who_data, year %in% c(2000:2004)) %>%
    group_by(ISO) %>%
    summarise(cases=mean(cases),
              deaths=mean(deaths),
              incidence_rate=mean(incidence_rate),
              mortality_rate=mean(mortality_rate),
              par=mean(population)) %>%
    left_join(par_latest, by="ISO") %>%
    mutate(burden_metric = (normalize(cases)+normalize(deaths)+0.05*normalize(incidence_rate)+0.05*normalize(mortality_rate))*(par_latest/par))
  # Weights need to be calculated after grouping countries by EIR to sum to 1
  
  if (species == "pf_pv") {
    # Combine match_dataset into a single EIR (sum of Pf and Pv)
    options(pillar.sigfigs=100)
    match_dataset <- match_dataset %>% 
      mutate(init_EIR = round(init_EIR_pf+init_EIR_pv,6),  # avoid dplyr rounding issues
             cases_pf = clinical_incidence_pf*pop_pf,
             cases_pv = clinical_incidence_pv*pop_pv) %>%
      select(-init_EIR_pf, -init_EIR_pv) %>%
      group_by(ISO,init_EIR) %>%
      summarise(pop_pf=sum(pop_pf),    
                pop_pv=sum(pop_pv),
                cases_pf = sum(cases_pf),
                cases_pv = sum(cases_pv)) %>%
      mutate(pop_admins0 = pmax(pop_pf, pop_pv),
             cases = cases_pf+cases_pv)
  }
  
  out <-select(country_burden_metrics, ISO, burden_metric) %>%
    full_join(match_dataset %>% dplyr::select(ISO,init_EIR,pop_admins0), by="ISO") %>%
    group_by(init_EIR) %>% filter(!is.na(init_EIR) & !is.na(burden_metric)) %>%
    # Ignore those countries that had no data to calculate a burden metric
    # (most have little malaria)
    summarize(pop = sum(pop_admins0, na.rm = T),
              burden_metric= sum(burden_metric, na.rm = T)) %>%
    mutate(w = burden_metric/sum(burden_metric,na.rm = T)) %>% 
    dplyr::select(init_EIR, w)
  
}

# read in data-sets with one entry per country per individual species
match_both <- readRDS(paste0(datadir, "country_match_both.rds"))

# calculate weights and make sure that totals sum to 1
weights_gf_both <- calculate_weights_gf(match_both, species = "pf_pv")
sum(weights_gf_both$w)   # check weights sum to 1 

## Manual allocation strategies for P. falciparum and P. vivax combined --------
# Allocation function  for cases and population at risk output -----------------
casefunction_pf_pv <- function(bednet, weights, netz, netz_longer, usedist){
  
  # number of total cases in the world (should be close to 300 million)
  totalcase <- bednet %>% filter(itn_cov==0) %>% ungroup() %>% summarize(totalcase=sum(cases)) %>% as.numeric()
  totalcase_pf <- bednet %>% filter(itn_cov==0) %>% ungroup() %>% summarize(totalcase=sum(cases_pf)) %>% as.numeric()
  totalcase_pv <- bednet %>% filter(itn_cov==0) %>% ungroup() %>% summarize(totalcase=sum(cases_pv)) %>% as.numeric()
  
  # number of total population in the world (should still be around 4 billion)
  totalpop_pf <- bednet %>% filter(itn_cov==0) %>% ungroup() %>% summarize(totalpop=sum(pop_pf)) %>% as.numeric()
  totalpop_pv <- bednet %>% filter(itn_cov==0) %>% ungroup() %>% summarize(totalpop=sum(pop_pv)) %>% as.numeric()
  totalpop <- max(totalpop_pf, totalpop_pv)
  
  # matching ITN distribution to use rate via netz
  netz <- netz %>% select(usage, annual_percapita_nets_distributed)
  
  # merging netz values to main dataset
  bednet <- bednet %>% left_join(netz, by=c('itn_cov'='usage'))
  
  # decide if a loess or linear relationship between use and distribution should be used
  bednet <- bednet %>% mutate(usedist = case_when(usedist == 'linear' ~ itn_cov, 
                                                  usedist == 'loess' ~ annual_percapita_nets_distributed))
  
  # Elimination for Pf and Pv combined:
  # threshold remains the same - if total cases <1, malaria is eliminated
  # but the population at risk may vary depending on whether only P. falciparum or
  # P. vivax are eliminated
  
  # HIGH BURDEN TO HIGH IMPACT
  # Funding is allocated starting from highest to lowest EIR, with increasing coverage
  # in each until eradication is reached
  output1bn <- bednet %>% 
    group_by(init_EIR) %>%
    arrange(-init_EIR, itn_cov) %>%   # sort by  decreasing EIR and increasing coverage
    # add cumulative case count
    mutate(rank = row_number(),
           changecase = cases - lag(cases),
           changecase = ifelse(is.na(changecase),0,changecase),
           cumulative = cumsum(replace_na(changecase, 0)),
           changecase_pf = cases_pf - lag(cases_pf),
           changecase_pf = ifelse(is.na(changecase_pf),0,changecase_pf),
           cumulative_pf = cumsum(replace_na(changecase_pf, 0)),
           changecase_pv = cases_pv - lag(cases_pv),
           changecase_pv = ifelse(is.na(changecase_pv),0,changecase_pv),
           cumulative_pv = cumsum(replace_na(changecase_pv, 0))) %>%
    ungroup() %>%
    mutate(tcumulative = cumsum(replace_na(changecase, 0)),  # cumulative cases averted
           tcumulative = totalcase + tcumulative,            # remaining case count by dollars spent
           tcumulative_pf = cumsum(replace_na(changecase_pf, 0)),  
           tcumulative_pf = totalcase_pf + tcumulative_pf,   
           tcumulative_pv = cumsum(replace_na(changecase_pv, 0)),  # cumulative cases averted
           tcumulative_pv = totalcase_pv + tcumulative_pv,   
           threshold_pf = ifelse(cases_pf >= elimination_threshold,1,0), 
           threshold_pv = ifelse(cases_pv >= elimination_threshold,1,0),
           threshold = ifelse(cases >= elimination_threshold,1,0)) %>%  # elimination of all malaria 
    group_by(init_EIR, threshold) %>%                                
    mutate(min_itn=ifelse(threshold == 0,min(itn_cov),max(bednet$itn_cov))) %>% # itn_cov at which all malaria is eliminated         
    ungroup() %>% 
    # population at risk and cumulative dollars spent
    mutate(popatrisk_pf=ifelse(cases_pf >= elimination_threshold, pop_pf, 0),   # pop. at risk is 0 if there is less than 1 case
           popatrisk_pv=ifelse(cases_pv >= elimination_threshold, pop_pv, 0),
           popatrisk=pmax(popatrisk_pf, popatrisk_pv)) %>%  
    group_by(init_EIR) %>%
    mutate(dollar = ifelse(itn_cov<=min_itn, usedist*pop*PYRcost, 0),
           diff_popatrisk=c(0,diff(popatrisk)),
           diff_popatrisk_pf=c(0,diff(popatrisk_pf)),
           diff_popatrisk_pv=c(0,diff(popatrisk_pv)),
           diff_dollar=ifelse(itn_cov<=min_itn,c(0,diff(dollar)),0)) %>%
    ungroup() %>%
    mutate(cum_diff_popatrisk=cumsum(diff_popatrisk),  
           cum_diff_popatrisk_pf=cumsum(diff_popatrisk_pf), 
           cum_diff_popatrisk_pv=cumsum(diff_popatrisk_pv), 
           tdollar=cumsum(diff_dollar),
           pcumulative_pf = totalpop_pf+cum_diff_popatrisk_pf,
           pcumulative_pv = totalpop_pv+cum_diff_popatrisk_pv,
           pcumulative = pmax(pcumulative_pf, pcumulative_pv)) %>%  
    select(init_EIR, itn_cov, pop_pf, pop_pv, pop, cases_pf, cases_pv, cases, rank, 
           tdollar, tcumulative, popatrisk, pcumulative, usedist, 
           tcumulative_pf, tcumulative_pv, pcumulative_pf, pcumulative_pv)
  
  # SHRINK THE MAP
  # Funding is allocated starting from lowest to highest EIR, with increasing coverage
  # in each until eradication is reached
  output2bn <- bednet %>% 
    group_by(init_EIR) %>%
    arrange(init_EIR, itn_cov) %>%   # sort by increasing EIR and increasing coverage
    # add cumulative case count
    mutate(rank = row_number(),
           changecase = cases - lag(cases),
           changecase = ifelse(is.na(changecase),0,changecase),
           cumulative = cumsum(replace_na(changecase, 0)),
           changecase_pf = cases_pf - lag(cases_pf),
           changecase_pf = ifelse(is.na(changecase_pf),0,changecase_pf),
           cumulative_pf = cumsum(replace_na(changecase_pf, 0)),
           changecase_pv = cases_pv - lag(cases_pv),
           changecase_pv = ifelse(is.na(changecase_pv),0,changecase_pv),
           cumulative_pv = cumsum(replace_na(changecase_pv, 0))) %>%
    ungroup() %>%
    mutate(tcumulative = cumsum(replace_na(changecase, 0)),  # cumulative cases averted
           tcumulative = totalcase + tcumulative,           # remaining case count by dollars spent
           tcumulative_pf = cumsum(replace_na(changecase_pf, 0)),  
           tcumulative_pf = totalcase_pf + tcumulative_pf,   
           tcumulative_pv = cumsum(replace_na(changecase_pv, 0)),  # cumulative cases averted
           tcumulative_pv = totalcase_pv + tcumulative_pv,  
           threshold_pf = ifelse(cases_pf >= elimination_threshold,1,0), 
           threshold_pv = ifelse(cases_pv >= elimination_threshold,1,0),
           threshold = ifelse(cases >= elimination_threshold,1,0)) %>% 
    group_by(init_EIR, threshold) %>%                                
    mutate(min_itn=ifelse(threshold == 0,min(itn_cov),max(bednet$itn_cov))) %>%  # itn_cov at which all malaria is eliminated  
    ungroup() %>% 
    # population at risk and cumulative dollars spent
    mutate(popatrisk_pf=ifelse(cases_pf >= elimination_threshold, pop_pf, 0),   # pop. at risk is 0 if there is less than 1 case
           popatrisk_pv=ifelse(cases_pv >= elimination_threshold, pop_pv, 0),
           popatrisk=pmax(popatrisk_pf, popatrisk_pv)) %>%  
    group_by(init_EIR) %>%
    mutate(dollar = ifelse(itn_cov<=min_itn, usedist*pop*PYRcost, 0),
           diff_popatrisk=c(0,diff(popatrisk)),
           diff_popatrisk_pf=c(0,diff(popatrisk_pf)),
           diff_popatrisk_pv=c(0,diff(popatrisk_pv)),
           diff_dollar=ifelse(itn_cov<=min_itn,c(0,diff(dollar)),0)) %>%
    ungroup() %>%
    mutate(cum_diff_popatrisk=cumsum(diff_popatrisk),  
           cum_diff_popatrisk_pf=cumsum(diff_popatrisk_pf), 
           cum_diff_popatrisk_pv=cumsum(diff_popatrisk_pv), 
           tdollar=cumsum(diff_dollar),
           pcumulative_pf = totalpop_pf+cum_diff_popatrisk_pf,
           pcumulative_pv = totalpop_pv+cum_diff_popatrisk_pv,
           pcumulative = pmax(pcumulative_pf, pcumulative_pv)) %>%  # remaining pop. at risk by dollars spent
    select(init_EIR, itn_cov, pop_pf, pop_pv, pop, cases_pf, cases_pv, cases, rank, 
           tdollar, tcumulative, popatrisk, pcumulative, usedist,
           tcumulative_pf, tcumulative_pv, pcumulative_pf, pcumulative_pv)
  
  # PROPORTIONAL ALLOCATION
  # Iterate through a sequence of budgets and calculate the $ amount allocated
  # to each setting for each based on the weight
  maxx <- output2bn %>% summarize(maxx = max(tdollar)) %>% as.numeric()
  budgets <- c(seq(0, maxx, 50000000),maxx)
  
  if (usedist == "linear") {
    output3bn <- list()
    for (i in 1:length(budgets)) {
      output3bn[[i]] <- bednet %>% 
        left_join(weights, by='init_EIR') %>%
        mutate(b=w*budgets[i],
               funded_cov =ifelse(round((b/PYRcost)/pop, 2)<=max(bednet$itn_cov), 
                                  round((b/PYRcost)/pop, 2), 
                                  max(bednet$itn_cov))) %>%
        group_by(init_EIR) %>%
        filter(itn_cov==funded_cov) %>%
        mutate(popatrisk_pf=ifelse(cases_pf >= elimination_threshold, pop_pf, 0),   # pop. at risk is 0 if there is less than 1 case
               popatrisk_pv=ifelse(cases_pv >= elimination_threshold, pop_pv, 0),
               popatrisk=pmax(popatrisk_pf, popatrisk_pv)) %>%
        ungroup() %>%
        summarise(tdollar= budgets[i],
                  tcumulative=sum(cases,na.rm=T), 
                  tcumulative_pf=sum(cases_pf,na.rm=T), 
                  tcumulative_pv=sum(cases_pv,na.rm=T), 
                  #pcumulative=sum(popatrisk,na.rm=T),
                  pcumulative_pf=sum(popatrisk_pf,na.rm=T),
                  pcumulative_pv=sum(popatrisk_pv,na.rm=T),
                  pcumulative = pmax(pcumulative_pf, pcumulative_pv))
      
    }
    
    output3bn <- do.call("rbind", output3bn)
    
  } else if (usedist == "loess") {
    
    output3bn <- list()
    for (i in 1:length(budgets)) {
      output3bn[[i]] <- bednet %>% 
        left_join(weights, by='init_EIR') %>%
        mutate(b=w*budgets[i],
               funded_nets = ifelse(round((b/PYRcost)/pop, 3)<=max(bednet$annual_percapita_nets_distributed), 
                                    round((b/PYRcost)/pop, 3), 
                                    max(bednet$annual_percapita_nets_distributed))) %>%
        fuzzyjoin::difference_left_join(as_tibble(netz_longer),
                                        by=c("funded_nets"="pc_nets_annual"), 
                                        max_dist=0.01, distance_col="dist") %>%
        group_by(init_EIR, funded_nets) %>% 
        slice_min(dist) %>%
        ungroup() %>%
        group_by(init_EIR) %>%
        filter(round(itn_cov,2)==round(usage,2)) %>%
        mutate(popatrisk_pf=ifelse(cases_pf >= elimination_threshold, pop_pf, 0),   # pop. at risk is 0 if there is less than 1 case
               popatrisk_pv=ifelse(cases_pv >= elimination_threshold, pop_pv, 0),
               popatrisk=pmax(popatrisk_pf, popatrisk_pv)) %>%
        ungroup() %>%
        summarise(tdollar= budgets[i],
                  tcumulative=sum(cases,na.rm=T), 
                  tcumulative_pf=sum(cases_pf,na.rm=T), 
                  tcumulative_pv=sum(cases_pv,na.rm=T), 
                  #pcumulative=sum(popatrisk,na.rm=T),
                  pcumulative_pf=sum(popatrisk_pf,na.rm=T),
                  pcumulative_pv=sum(popatrisk_pv,na.rm=T),
                  pcumulative = pmax(pcumulative_pf, pcumulative_pv))
    }
    
    output3bn <- do.call("rbind", output3bn)
  }
  
  
  join_all <- function(data, label){
    data %>% 
      mutate(scenario=label)
  }
  
  # merge absolute difference outputs together
  output <- join_all(output1bn, "high burden to high impact") %>% 
    full_join(join_all(output2bn, "shrink the map")) %>%
    full_join(join_all(output3bn, "proportional allocation"))
  
  return(output)
  
}

# Apply allocation function -----------------------------------------------------------

# up to 80% ITN coverage - falciparum+vivax - linear costing assumption
alloc_pf_pv_linear <- casefunction_pf_pv(bednet = sims_pf_pv  %>% filter(itn_cov <= 0.80),
                                     weights =  weights_gf_both,
                                     netz = readRDS(paste0(datadir, "netz_median_use_rate.rds")), 
                                     netz_longer = readRDS(paste0(datadir, "conversion_usage_pcnets_median_use_rate.rds")),
                                     usedist = "linear")

saveRDS(alloc_pf_pv_linear, paste0(outdir, "bednets_data_equilibrium_pf_pv_80_linear.rds"))

# up to 80% ITN coverage - falciparum+vivax - loess costing assumption
alloc_pf_pv_loess <- casefunction_pf_pv(bednet = sims_pf_pv  %>% filter(itn_cov <= 0.80),
                                     weights =  weights_gf_both,
                                     netz = readRDS(paste0(datadir, "netz_median_use_rate.rds")), 
                                     netz_longer = readRDS(paste0(datadir, "conversion_usage_pcnets_median_use_rate.rds")),
                                     usedist = "loess")

saveRDS(alloc_pf_pv_loess, paste0(outdir, "bednets_data_equilibrium_pf_pv_80_loess.rds"))

# Check output  -----------------------------------------------------------------
# Impact on global cases
ggplot(alloc_pf_pv_loess) +
  geom_line(aes(x=tdollar, y = tcumulative, group = scenario, colour=scenario)) +
  labs(x="Budget", y="Total malaria cases", colour = "Scenario") +
  theme_classic() 

# Impact on global population at risk
ggplot(alloc_pf_pv_loess) +
  geom_line(aes(x=tdollar, y = pcumulative, group = scenario, colour=scenario)) + 
  labs(x="Budget", y="Total population at risk of malaria", colour = "Scenario") +
  theme_classic()

# Impact on global cases by species
ggplot(alloc_pf_pv_loess) +
  geom_line(aes(x=tdollar/1000000, y = tcumulative_pf, colour = "Pf")) +
  geom_line(aes(x=tdollar/1000000, y = tcumulative_pv, colour = "Pv")) +
  facet_wrap(~scenario) +
  labs(x="Budget (millions)", y="Species-specific clinical cases", colour = "Species") +
  theme_classic()

# Impact on global population at risk by species
ggplot(alloc_pf_pv_loess) +
  geom_line(aes(x=tdollar, y = pcumulative_pf, colour = "Pf")) +
  geom_line(aes(x=tdollar, y = pcumulative_pv, colour = "Pv")) +
  geom_line(aes(x=tdollar, y = pcumulative, colour = "Total"), col = "black") +
  facet_wrap(~scenario) +
  labs(x="Budget", y="Species-specific population at risk", colour = "Species") + 
  theme_classic()

# Allocation function for proportion of budget allocated -----------------------

# This function takes as input the output of casefunction() (manual) and the original inputs
# (bednets_data_equilibrium_pf_pv, weights_gf_both) for the proportional allocation scenario
# Takes inputs for only a single year
calculate_allocation_prop <- function(manual, bednet, weights, netz, netz_longer, 
                                      usedist_var) {
  
  # matching ITN distribution to use rate via netz
  netz <- netz %>% select(usage, annual_percapita_nets_distributed)
  
  # merging netz values to main dataset
  bednet <- bednet %>% left_join(netz, by=c('itn_cov'='usage'))
  
  # decide if a loess or linear relationship between use and distribution should be used
  bednet <- bednet %>% mutate(usedist = case_when(usedist_var == 'linear' ~ itn_cov, 
                                                  usedist_var == 'loess' ~ annual_percapita_nets_distributed))
  
  manual <- manual %>% mutate(usedist = case_when(usedist_var == 'linear' ~ itn_cov, 
                                                  usedist_var == 'loess' ~ usedist))
  
  
  # HIGH BURDEN TO HIGH IMPACT
  scenario1 <- filter(manual, scenario == "high burden to high impact") %>%
    mutate(dollar=usedist*pop*PYRcost) %>%
    group_by(init_EIR) %>%
    mutate(changedollar=c(0,diff(tdollar))) %>%
    filter(changedollar != 0)  # need to exclude rows where no further investment is made
  
  prop_allocated1 <- list()
  for (i in seq(1,length(scenario1$tdollar),10)) {
    prop_allocated1[[i]] <- scenario1 %>% filter(tdollar<=scenario1$tdollar[i]) %>%
      ungroup() %>%
      mutate(B=max(tdollar)) %>%
      group_by(init_EIR) %>%
      filter(itn_cov==max(itn_cov)) %>%
      summarise(scenario=scenario,
                B=B, 
                b=sum(dollar),
                prop=b/B,
                final_cov=itn_cov)
  }
  prop_allocated1 <- do.call("rbind", prop_allocated1)
  # Fill in 0s:
  prop_allocated1 <-full_join(prop_allocated1,
                              crossing(init_EIR=unique(prop_allocated1$init_EIR),
                                       scenario="high burden to high impact",
                                       B=c(0,unique(prop_allocated1$B)))) %>%
    replace_na(list(b=0, prop=0, final_cov=0))
  
  
  # Test
  if(
    (any(between((filter(prop_allocated1, B!=0) %>% group_by(B) %>% summarise(prop=sum(prop))     %>% pull(prop)),0.99, 1.01) == FALSE)) == TRUE) {
    stop("Proportions not adding up to 1")
  }
  
  # SHRINK THE MAP
  scenario2 <- filter(manual, scenario == "shrink the map") %>%
    mutate(dollar=usedist*pop*PYRcost) %>%
    group_by(init_EIR) %>%
    mutate(changedollar=c(0,diff(tdollar))) %>%
    filter(changedollar != 0)  # need to exclude rows where no further investment is made
  
  prop_allocated2 <- list()
  for (i in seq(1,length(scenario2$tdollar),10)) {
    prop_allocated2[[i]] <- scenario2 %>% filter(tdollar<=scenario2$tdollar[i]) %>%
      ungroup() %>%
      mutate(B=max(tdollar)) %>%
      group_by(init_EIR) %>%
      filter(itn_cov==max(itn_cov)) %>%
      summarise(scenario=scenario,
                B=B, 
                b=sum(dollar),
                prop=b/B,
                final_cov=itn_cov) 
  }
  prop_allocated2 <- do.call("rbind", prop_allocated2)
  # Fill in 0s:
  prop_allocated2 <-full_join(prop_allocated2,
                              crossing(init_EIR=unique(prop_allocated2$init_EIR),
                                       scenario="shrink the map",
                                       B=c(0,unique(prop_allocated2$B)))) %>%
    replace_na(list(b=0, prop=0, final_cov=0))
  
  # Test
  if(
    (any(between((filter(prop_allocated2, B!=0) %>% group_by(B) %>% summarise(prop=sum(prop))     %>% pull(prop)),0.99, 1.01) == FALSE)) == TRUE) {
    stop("Proportions not adding up to 1")
  }
  
  # PROPORTIONAL ALLOCATION (different method - same as in casefunction()) 
  maxx <- manual %>% summarize(maxx = max(tdollar)) %>% as.numeric()
  budgets <- c(seq(0, maxx, 10000000),maxx)
  
  if (usedist_var == "linear") {
    prop_allocated3 <- list()
    for (i in 1:length(budgets)) {
      prop_allocated3[[i]] <- bednet %>% 
        left_join(weights, by='init_EIR') %>%
        mutate(b=w*budgets[i],
               funded_cov =ifelse(round((b/PYRcost)/pop, 2)<=max(bednet$itn_cov), 
                                  round((b/PYRcost)/pop, 2), 
                                  max(bednet$itn_cov))) %>%
        group_by(init_EIR) %>%
        filter(itn_cov==funded_cov) %>%
        mutate(B= budgets[i],
               prop = b/B,
               scenario="proportional allocation",
               final_cov=funded_cov) %>%
        select(init_EIR, scenario, B, b, prop, final_cov) %>%
        replace_na(list(prop=0))
    }
    prop_allocated3 <- do.call("rbind", prop_allocated3)
    # final_cov here can go up to maximum even if this wasn't needed to achieve the
    # respective cases
    
  } else if (usedist_var == "loess") {
    
    prop_allocated3 <- list()
    for (i in 1:length(budgets)) {
      prop_allocated3[[i]] <- bednet %>% 
        left_join(weights, by='init_EIR') %>%
        mutate(b=w*budgets[i],
               funded_nets = ifelse(round((b/PYRcost)/pop, 3)<=max(bednet$annual_percapita_nets_distributed), 
                                    round((b/PYRcost)/pop, 3), 
                                    max(bednet$annual_percapita_nets_distributed))) %>%
        fuzzyjoin::difference_left_join(as_tibble(netz_longer),
                                        by=c("funded_nets"="pc_nets_annual"), 
                                        max_dist=0.01, distance_col="dist") %>%
        group_by(init_EIR, funded_nets) %>% 
        slice_min(dist) %>%
        ungroup() %>%
        group_by(init_EIR) %>%
        filter(round(itn_cov,2)==round(usage,2)) %>%
        mutate(B= budgets[i],
               prop = b/B,
               scenario="proportional allocation",
               final_cov=usage) %>%
        select(init_EIR, scenario, B, b, prop, final_cov) %>%
        replace_na(list(prop=0))
    }
    prop_allocated3 <- do.call("rbind", prop_allocated3)
    
  }
  
  # Test
  if(
    (any(between((filter(prop_allocated3, B!=0) %>% group_by(B) %>% summarise(prop=sum(prop))     %>% pull(prop)),0.99, 1.01) == FALSE)) == TRUE) {
    stop("Proportions not adding up to 1")
  }
  
  out <- rbind(prop_allocated1,prop_allocated2, prop_allocated3)
  
  return(out)
  
}

# Apply allocation pattern function --------------------------------------------

# up to 80% ITN coverage - falciparum+vivax - linear
allocation_pf_pv_linear80 <- calculate_allocation_prop(
  manual=readRDS(paste0(outdir, "bednets_data_equilibrium_pf_pv_80_linear.rds")), 
  bednet=sims_pf_pv  %>% filter(itn_cov <= 0.80),   
  weights=weights_gf_both,
  netz = readRDS(paste0(datadir, "netz_median_use_rate.rds")),
  netz_longer = readRDS(paste0(datadir, "conversion_usage_pcnets_median_use_rate.rds")) ,
  usedist_var = 'linear')

# up to 80% ITN coverage - falciparum+vivax - loess
allocation_pf_pv_loess80 <- calculate_allocation_prop(
  manual=readRDS(paste0(outdir, "bednets_data_equilibrium_pf_pv_80_loess.rds")), 
  bednet=sims_pf_pv  %>% filter(itn_cov <= 0.80),   
  weights=weights_gf_both,
  netz = readRDS(paste0(datadir, "netz_median_use_rate.rds")),
  netz_longer = readRDS(paste0(datadir, "conversion_usage_pcnets_median_use_rate.rds")) ,
  usedist_var = 'loess')

# Check output ----

# Function for plotting the proportion of budget allocated to each EIR or
# the funded itn_cov
plot_allocation_pattern <- function(dataset, outcome="prop", n_eir_groups = "WHO") {
  # specify outcome as "prop" or "itn_cov" 
  # specify the number of EIR groups to differentiate by colour on plots 
  # or select WHO classification
  
  if (n_eir_groups == "WHO") {
    allocation_df <- dataset %>%
      mutate(eir_range = case_when(init_EIR < 0.1 ~ 'Very low',
                                   init_EIR >= 0.1 & init_EIR < 1 ~ 'Low',
                                   init_EIR >= 1 & init_EIR < 7 ~ 'Moderate',
                                   init_EIR >= 7 ~ 'High'),
             eir_range = factor(eir_range, levels = 
                                  c("Very low", "Low", "Moderate", "High"))) 
  } else {
    eir_groups <- data.frame(init_EIR=unique(dataset$init_EIR)) 
    eir_groups <- eir_groups %>%
      mutate(quartile_eir = ntile(init_EIR, 4))
    
    allocation_df <- dataset %>%
      ungroup() %>%
      left_join(eir_groups, by = "init_EIR") %>%
      group_by(quartile_eir) %>%
      mutate(min_eir = min(init_EIR),
             eir_range = paste0(min(init_EIR),"-",max(init_EIR)))
  }
  
  allocation_df$scenario <- factor(allocation_df$scenario)
  levels(allocation_df$scenario) <- list("Prioritize high-transmission settings" = 
                                           "high burden to high impact",
                                         "Proportional allocation" =
                                           "proportional allocation",
                                         "Prioritize low-transmission settings" = 
                                           "shrink the map")
  
  spectral_colours <- c("#D53E4F", "#FC8D59", "#FEE08B", "#FFFFBF", "#E6F598", 
                        "#99D594", "#3288BD")
  
  if(outcome=="prop") {
    
    ggplot(allocation_df %>% group_by(B, scenario, eir_range) %>% summarise(b=sum(b))) +
      geom_area(aes(x=B/1000000000, y = b/B*100, fill=eir_range, 
                    group = eir_range), size=1) +
      facet_wrap(~scenario) +
      scale_fill_manual(values= spectral_colours[c(7,6,3,1)]) + 
      labs(fill="Transmission intensity") +
      ylab("Proportion of total budget (%)") +
      xlab("Budget (billion $)") +
      xlim(c(0,max(dataset$B)/1000000000)) +
      theme_classic() +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
    
  } else if (outcome=="itn_cov") {
    
    ggplot(allocation_df %>% group_by(B, scenario, eir_range) %>% summarise(usage=mean(final_cov))) +
      geom_line(aes(x=B/1000000000, y = usage*100, colour=eir_range, 
                    group = eir_range), size=1) +
      facet_wrap(~scenario) +
      scale_colour_manual(values= spectral_colours[c(7,6,3,1)]) + 
      labs(colour="Transmission intensity") +
      ylab("Funded mean ITN usage (%)") +
      xlab("Budget (million $)") +
      xlim(c(0,max(dataset$B)/1000000000)) +
      theme_classic() +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
    
  } else if (outcome=="both") {
    
    allocation_df <- allocation_df %>% 
      group_by(B, scenario, eir_range) %>% 
      summarise(b=sum(b),
                usage=mean(final_cov)) %>%
      mutate(prop = b/B) %>%
      replace_na(list(prop=0)) 
    
    p1 <- ggplot(allocation_df) +
      geom_area(aes(x=B/1000000000, y = prop*100, fill=eir_range, 
                    group = eir_range), size=1) +
      facet_wrap(~scenario) +
      scale_fill_manual(values= spectral_colours[c(7,6,3,1)]) + 
      labs(fill="Transmission intensity") +
      ylab("Proportion of total budget (%)") +
      xlab("Budget (billion $)") +
      xlim(c(0,max(dataset$B)/1000000000)) +
      theme_classic() +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 
    
    p2 <- ggplot(allocation_df) +
      geom_line(aes(x=B/1000000000, y = usage*100, colour=eir_range, 
                    group = eir_range), size=1) +
      facet_wrap(~scenario) +
      scale_colour_manual(values= spectral_colours[c(7,6,3,1)]) + 
      labs(colour="Transmission intensity") +
      ylab("Funded mean ITN usage (%)") +
      xlab("Budget (billion $)") +
      xlim(c(0,max(dataset$B)/1000000000)) +
      theme_classic() +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            legend.position = "none") 
    
    (p1 + p2) + plot_layout(guides = "collect", nrow=2) + 
      plot_annotation(tag_levels = c("A", "B")) 
    
  } else {
    print("outcome has to be prop or itn_cov or both")
  }
  
} 

plot_allocation_pattern(allocation_pf_pv_loess80, outcome="both", n_eir_groups ="WHO")
plot_allocation_pattern(allocation_pf_pv_linear80, outcome="both", n_eir_groups ="WHO")

