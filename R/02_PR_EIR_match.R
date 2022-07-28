# Match PfPR and PvPR to EIR values --------------------------------------------
library(tidyverse)
library(here)

# working directory
datadir <- here("files", "/")

# One-off scaling of clinical incidence to represent weekly active case detection 
# instead of daily (DO NOT RE-RUN)

# Applied to simulations of clinical incidence detected using daily active case detection
# generated in "example_code_run_simulations.R": dmm_output_foi_bednets_pf.rds and
# dmm_output_foi_bednets_pv.rds

# Scaling factor for P. falciparum from Griffin et al 2014 (Nat Comms)
# For P. vivax from Battle et al 2016 (Defining the relationship between Plasmodium vivax parasite rate and clinical disease, Malaria Journal)
# after White et al 2018 (Nat Comms)

# dmm_output_foi_bednets_pf <- readRDS(file.path(datadir, "dmm_output_foi_bednets_pf.rds"))
# saveRDS(dmm_output_foi_bednets_pf, file.path(datadir, "dmm_output_foi_bednets_pf_daily_ACD.rds"))
# dmm_output_foi_bednets_pf_weekly_ACD <- dmm_output_foi_bednets_pf %>%
#    mutate(clinical_incidence = clinical_incidence*0.723,
#            clinical_incidence_05 = clinical_incidence_05*0.723)
# saveRDS(dmm_output_foi_bednets_pf_weekly_ACD, file.path(datadir, "dmm_output_foi_bednets_pf.rds"))
# 
# dmm_output_foi_bednets_yearly_pf <- readRDS(file.path(datadir, "dmm_output_foi_bednets_yearly_pf.rds"))
# saveRDS(dmm_output_foi_bednets_yearly_pf, file.path(datadir, "dmm_output_foi_bednets_yearly_pf_daily_ACD.rds"))
# dmm_output_foi_bednets_yearly_pf_weekly_ACD <- dmm_output_foi_bednets_yearly_pf %>%
#    mutate(clinical_incidence = clinical_incidence*0.723,
#          clinical_incidence_05 = clinical_incidence_05*0.723)
# saveRDS(dmm_output_foi_bednets_yearly_pf_weekly_ACD, file.path(datadir, "dmm_output_foi_bednets_yearly_pf.rds"))

# dmm_output_foi_bednets_pv <- readRDS(file.path(datadir, "dmm_output_foi_bednets_pv.rds"))
# saveRDS(dmm_output_foi_bednets_pv, file.path(datadir, "dmm_output_foi_bednets_pv_daily_ACD.rds"))
# dmm_output_foi_bednets_pv_weekly_ACD <- dmm_output_foi_bednets_pv %>%
#    mutate(clinical_incidence = clinical_incidence*0.1336,
#           clinical_incidence_05 = clinical_incidence_05*0.1336)
# saveRDS(dmm_output_foi_bednets_pv_weekly_ACD, file.path(datadir, "dmm_output_foi_bednets_pv.rds"))
# 
# dmm_output_foi_bednets_yearly_pv <- readRDS(file.path(datadir, "dmm_output_foi_bednets_yearly_pv.rds"))
# saveRDS(dmm_output_foi_bednets_yearly_pv, file.path(datadir, "dmm_output_foi_bednets_yearly_pv_daily_ACD.rds"))
# dmm_output_foi_bednets_yearly_pv_weekly_ACD <-dmm_output_foi_bednets_yearly_pv %>%
#    mutate(clinical_incidence = clinical_incidence*0.1336,
#         clinical_incidence_05 = clinical_incidence_05*0.1336)
# saveRDS(dmm_output_foi_bednets_yearly_pv_weekly_ACD, file.path(datadir, "dmm_output_foi_bednets_yearly_pv.rds"))

# Matching ---------------------------------------------------------------------
# read in simulation results: ITN usage and EIR 
sims_pf <- readRDS(file.path(datadir, "dmm_output_foi_bednets_pf.rds")) %>%  
  select(init_EIR, prev210, clinical_incidence, clinical_incidence_05, itn_cov) %>%  
  rename(init_EIR_pf=init_EIR, prev_pf=prev210, clinical_incidence_pf=clinical_incidence, clinical_incidence_05_pf=clinical_incidence_05)

sims_pv <- readRDS(file.path(datadir, "dmm_output_foi_bednets_pv.rds")) %>%
  select(init_EIR, prevLM, clinical_incidence, clinical_incidence_05, itn_cov) %>%  
  rename(init_EIR_pv=init_EIR, prev_pv=prevLM, clinical_incidence_pv=clinical_incidence, clinical_incidence_05_pv=clinical_incidence_05)

# save outputs when ITN usage == 0
outputs_pf <- sims_pf %>%
  filter(itn_cov==0) %>% select(-itn_cov)

outputs_pv <- sims_pv %>% 
  filter(itn_cov==0) %>% select(-itn_cov)

# read in geoboundary country outlines and prevalence values
admin0_pf <- readRDS(file.path(datadir, "admins0_all_2000_pf.rds"))
admin0_pv <- readRDS(file.path(datadir, "admins0_all_2000_pv.rds"))


# To each country assign one Pf EIR and one Pv EIR:
admin0_combined <- full_join(admin0_pf, 
                             admin0_pv, 
                             by = c("ISO", "Country","region")) %>%
  replace_na(list(PR_admins0_pf=0, PR_admins0_pv=0, pop_pf=0, pop_pv=0)) %>%
# remove inconsistencies between population at risk and prevalence
  mutate(pop_pf = ifelse(PR_admins0_pf==0,0,pop_pf),
         PR_admins0_pf = ifelse(pop_pf==0,0,PR_admins0_pf),
         pop_pv = ifelse(PR_admins0_pv==0,0,pop_pv),
         PR_admins0_pv = ifelse(pop_pv==0,0,PR_admins0_pv))

match <- admin0_combined %>% 
  fuzzyjoin::difference_left_join(outputs_pf,
                                  by=c("PR_admins0_pf"="prev_pf"), 
                                  max_dist=10, distance_col="dist_pf") %>%
  group_by(Country) %>% slice_min(dist_pf) %>%
  
  fuzzyjoin::difference_left_join(outputs_pv,
                                  by=c("PR_admins0_pv"="prev_pv"), 
                                  max_dist=10, distance_col="dist_pv") %>%
  group_by(Country) %>% slice_min(dist_pv) %>% ungroup() %>%
  
  # correct wrongly matched EIRS where a prev of 0 matched to an EIR of 0.001
  mutate(init_EIR_pf=ifelse(PR_admins0_pf==0, 0, init_EIR_pf),
         init_EIR_pv=ifelse(PR_admins0_pv==0, 0, init_EIR_pv),
         prev_pf=ifelse(PR_admins0_pf==0, 0, prev_pf),
         prev_pv=ifelse(PR_admins0_pv==0, 0, prev_pv),
         clinical_incidence_pf=ifelse(PR_admins0_pf==0, 0, clinical_incidence_pf),
         clinical_incidence_pv=ifelse(PR_admins0_pv==0, 0, clinical_incidence_pv))

# Now we have a Pf and Pv EIR and a population for each country
# From this can calculate, for each country, the total clinical incidence (Pf+Pv) at each EIR

# Add the clinical incidence values for P. vivax and P. falciparum for all itn_cov levels
# Do this separately depending on whether countries have Pf, Pv or both

# first make versions of datasets with ALL encompassing Pf, Pv, and both species
country_match_pf <- match %>% filter(init_EIR_pf!=0) %>% 
  group_by(ISO, Country, region) %>%
  mutate(init_EIR = init_EIR_pf,
         pop_admins0 = pop_pf,
         clinical_incidence = clinical_incidence_pf,
         cases = clinical_incidence_pf*pop_pf)

country_match_pv <- match %>% filter(init_EIR_pv!=0) %>%
  group_by(ISO, Country, region) %>%
  mutate(init_EIR = init_EIR_pv,
         pop_admins0 = pop_pv,
         clinical_incidence = clinical_incidence_pv,
         cases = clinical_incidence_pv*pop_pv)

country_match_both <- match 

saveRDS(country_match_pf, file= paste0(datadir, "country_match_pf.rds"))
saveRDS(country_match_pv, file= paste0(datadir, "country_match_pv.rds"))
saveRDS(country_match_both, file= paste0(datadir, "country_match_both.rds"))

# Create separate Pf and Pf datasets

# matches grouped by EIR
prevmatch_pf <- country_match_pf %>%
  group_by(init_EIR) %>%
  summarize(n=n(), pop=sum(pop_pf)) %>%
  select(init_EIR, n, pop)

prevmatch_pv <- country_match_pv %>%
  group_by(init_EIR) %>%
  summarize(n=n(), pop=sum(pop_pv)) %>%
  select(init_EIR, n, pop)

# link to clinical incidence
datapf <- sims_pf %>% 
  rename(init_EIR = init_EIR_pf,
         clinical_incidence_bednet = clinical_incidence_pf,
         prev = prev_pf) %>%
  left_join(prevmatch_pf, by='init_EIR') %>%
  filter(!is.na(pop)) %>%
  mutate(cases = clinical_incidence_bednet * pop) %>%
  select(init_EIR, pop, itn_cov, clinical_incidence_bednet, cases, prev, n)

saveRDS(datapf, file= paste0(datadir, "dmm_output_foi_bednets_wpop_pf.rds"))

datapv <- sims_pv %>% 
  rename(init_EIR = init_EIR_pv,
         clinical_incidence_bednet = clinical_incidence_pv,
         prev = prev_pv) %>%
  left_join(prevmatch_pv, by='init_EIR') %>%
  filter(!is.na(pop)) %>%
  mutate(cases = clinical_incidence_bednet * pop) %>%
  select(init_EIR, pop, itn_cov, clinical_incidence_bednet, cases, prev, n)

saveRDS(datapv, file= paste0(datadir, "dmm_output_foi_bednets_wpop_pv.rds"))

# Create 3 year output dataset for Pf
sims_pf_3 <- readRDS(file.path(datadir, "dmm_output_foi_bednets_yearly_pf.rds")) %>% 
  filter(years == 3) %>%
  select(years, init_EIR, prev210, clinical_incidence, clinical_incidence_05, itn_cov) 

datapf_3 <- sims_pf_3 %>% 
  rename(clinical_incidence_bednet = clinical_incidence,
         prev = prev210) %>%
  left_join(prevmatch_pf, by='init_EIR') %>%
  filter(!is.na(pop)) %>%
  mutate(cases = clinical_incidence_bednet * pop) %>%
  select(init_EIR, pop, itn_cov, clinical_incidence_bednet, cases, prev, n)
#saveRDS(datapf_3, file= paste0(datadir, "/bednets_data_year3_pf.rds"))

# Create 39 year output dataset for Pf
sims_pf_39 <- readRDS(file.path(datadir, "dmm_output_foi_bednets_yearly_pf.rds")) %>% 
  filter(years == 39) %>%
  select(years, init_EIR, prev210, clinical_incidence, clinical_incidence_05, itn_cov) 

datapf_39 <- sims_pf_39 %>% 
  rename(clinical_incidence_bednet = clinical_incidence,
         prev = prev210) %>%
  left_join(prevmatch_pf, by='init_EIR') %>%
  filter(!is.na(pop)) %>%
  mutate(cases = clinical_incidence_bednet * pop) %>%
  select(init_EIR, pop, itn_cov, clinical_incidence_bednet, cases, prev, n)
#saveRDS(datapf_39, file= paste0(datadir, "/bednets_data_year39_pf.rds"))

# Create a combined Pf and Pv dataset

# now make versions of datasets grouped by EIR instead of by country
countries_pf <- filter(match, init_EIR_pf!=0 & init_EIR_pv==0)
countries_pv <- filter(match, init_EIR_pf==0 & init_EIR_pv!=0)
countries_both <- filter(match, init_EIR_pf!=0 & init_EIR_pv!=0)

# now join all simulation runs for ITN usage 0 to 1
join_countries_pf <- countries_pf %>%
  select(ISO, Country, pop_pf, region, init_EIR_pf, init_EIR_pv) %>%
  full_join(sims_pf, by = "init_EIR_pf") %>%
  mutate(clinical_incidence_pv=0, clinical_incidence_05_pv=0, pop_pv=0, prev_pv=0) %>%
  drop_na(Country)

join_countries_pv <- countries_pv %>%
  select(ISO, Country, pop_pv, region, init_EIR_pf, init_EIR_pv) %>%
  full_join(sims_pv, by = "init_EIR_pv") %>%
  mutate(clinical_incidence_pf=0, clinical_incidence_05_pf=0, pop_pf=0, prev_pf=0) %>%
  drop_na(Country)

join_countries_both <- countries_both %>%
  select(ISO, Country, pop_pf, pop_pv, region, init_EIR_pf, init_EIR_pv) %>%
  full_join(sims_pf, by = "init_EIR_pf") %>%
  full_join(sims_pv, by = c("init_EIR_pv", "itn_cov")) %>%
  drop_na(Country)

# Combine all
all_countries_df <- rbind(join_countries_pf, join_countries_pv, join_countries_both) %>%
  mutate(cases_pf=clinical_incidence_pf*pop_pf,
         cases_pv=clinical_incidence_pv*pop_pv,
         cases_both=cases_pf+cases_pv) 

# Total populations at risk and # of EIRs are in line with individual datasets
sum(filter(all_countries_df, itn_cov==0)$pop_pf)==sum(filter(datapf, itn_cov==0)$pop)
sum(filter(all_countries_df, itn_cov==0)$pop_pv)==sum(filter(datapv, itn_cov==0)$pop)
# PAR for Pf = 4,147,518,844
# PAR for Pv = 3,954,391,575

# For the optimisation, reduce this to a dataset characterised by distinct combinations
# of Pf and Pv EIR. Total population in each needs to then combine several countries.

# matches grouped by EIR
group_countries <- all_countries_df %>%
  select(init_EIR_pf, init_EIR_pv, pop_pf, pop_pv, itn_cov) %>%
  filter(itn_cov==0) %>%
  group_by(init_EIR_pf, init_EIR_pv) %>%
  summarize(n=n(), pop_pf=sum(pop_pf), pop_pv=sum(pop_pv)) %>%
  ungroup()

data <- select(all_countries_df, init_EIR_pf, init_EIR_pv, itn_cov, clinical_incidence_pf,
               clinical_incidence_pv, prev_pf, prev_pv) %>%
  left_join(group_countries, by=c("init_EIR_pf", "init_EIR_pv")) %>%
  distinct() %>%
  mutate(cases_pf = clinical_incidence_pf * pop_pf,
         cases_pv = clinical_incidence_pv * pop_pv)

saveRDS(data, file= paste0(datadir, "dmm_output_foi_bednets_wpop_pf_pv.rds"))

# Checks ------------------------------------------------------------------------

# check total population and total case load
datapf %>% filter(itn_cov==0) %>% ungroup() %>% 
  summarize(tpop = sum(pop), tcase = sum(cases), n = sum(n))
# For weekly ACD, global cases = 251,989,681 (with treatment)

datapv %>% filter(itn_cov==0) %>% ungroup() %>% 
  summarize(tpop = sum(pop), tcase = sum(cases), n = sum(n))
# For weekly ACD, global cases = 69,338,173 (with treatment)

# Check this is the same in combined dataset:
sum(filter(data, itn_cov==0)$cases_pf)   
sum(filter(data, itn_cov==0)$cases_pv)   

# Group by country and check that the ITN range is from 0-1
summary(join_countries_pf$itn_cov)
summary(join_countries_pv$itn_cov)
summary(join_countries_both$itn_cov)

# Check for 105 individual countries
(length(unique(join_countries_pf$Country)) +  length(unique(join_countries_pv$Country)) +
  length(unique(join_countries_both$Country))) == 105

length(unique(c(unique(join_countries_pf$Country), unique(join_countries_pv$Country), 
                unique(join_countries_both$Country)))) == 105

# Check all EIR combinations have 101 rows for itn_cov (should return FALSE)
any((group_by(data, init_EIR_pf, init_EIR_pv) %>% summarise(rows=n()) %>% 
       ungroup() %>% select(rows)) != 101)

