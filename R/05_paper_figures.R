# Set-up and functions ---------------------------------------------------------

# packages
library(tidyverse)
library(fuzzyjoin)
library(ICDMM)
library(akima)
library(here)
library(patchwork)
library(ggrepel)
library(gridExtra)
library(plotly)
library(ggpubr)

# working directory
datadir <- here("files", "/")
outdir <- here("output", "/")

transmission_intensity_group_colours <- 
  c("#3288BD","#99D594", "#FEE08B", "#D53E4F")
# Note WHO classification of transmission intensities was derived by 
# matching sum of Pf+Pv prev to sum of Pf+Pv EIR (and minor tweak for low/moderate grouping)
# Prevalence grouping in https://www.who.int/publications/i/item/guidelines-for-malaria 

pf_pv_group_colours <- c("#FDAE6B", "#E6550D", "#FCCC9C")  # For Pf only,both, Pv only

# Functions to load
combine_strategies_impact <- function(optimizecases, manual) {
  # totals
  totalcases <-  manual %>% filter(scenario == 'high burden to high impact') %>% filter(itn_cov==0) %>% summarize(cases=sum(cases, na.rm=T)) %>% as.numeric()
  totalcases_pf <-  manual %>% filter(scenario == 'high burden to high impact') %>% filter(itn_cov==0) %>% summarize(cases=sum(cases_pf, na.rm=T)) %>% as.numeric()
  totalcases_pv <-  manual %>% filter(scenario == 'high burden to high impact') %>% filter(itn_cov==0) %>% summarize(cases=sum(cases_pv, na.rm=T)) %>% as.numeric()
  
  totalpop_pf <-  manual %>% filter(scenario == 'high burden to high impact') %>% filter(itn_cov==0) %>% summarize(tpop=sum(pop_pf, na.rm=T)) %>% as.numeric()
  totalpop_pv <-  manual %>% filter(scenario == 'high burden to high impact') %>% filter(itn_cov==0) %>% summarize(tpop=sum(pop_pv, na.rm=T)) %>% as.numeric()
  totalpop <- pmax(totalpop_pf, totalpop_pv)
  
  totaldollar <- manual %>% filter(scenario == 'high burden to high impact') %>% summarize(t=max(tdollar)) %>% as.numeric()
  
  optimizecases <- optimizecases %>% 
    group_by(B, b_sum) %>% 
    summarise(tcumulative=sum(cases, na.rm=T),
              tcumulative_pf=sum(cases_pf, na.rm=T),
              tcumulative_pv=sum(cases_pv, na.rm=T),
              pcumulative_pf=sum(popatrisk_pf, na.rm=T),
              pcumulative_pv=sum(popatrisk_pv, na.rm=T),
              pcumulative=pmax(pcumulative_pf, pcumulative_pv)) %>% 
    mutate(tdollar=b_sum,
           scenario = 'optimization: total cases') %>% distinct() %>%
    filter(round(tdollar,0)<=round(max(manual$tdollar),0))
  
  # Manually extend line to maximum dollar amount and 0
  optimizecases <- rbind(optimizecases,
                         data.frame(B = 0,
                                    tcumulative=totalcases, 
                                    pcumulative=totalpop, 
                                    tcumulative_pf=totalcases_pf, 
                                    tcumulative_pv=totalcases_pv, 
                                    pcumulative_pf=totalpop_pf,
                                    pcumulative_pv=totalpop_pv,
                                    tdollar=0, 
                                    scenario="optimization: total cases"),
                         data.frame(B = totaldollar,
                                    tcumulative=min(manual$tcumulative), 
                                    pcumulative=min(manual$pcumulative), 
                                    tcumulative_pf=min(manual$tcumulative_pf), 
                                    tcumulative_pv=min(manual$tcumulative_pv), 
                                    pcumulative_pf=min(manual$pcumulative_pf), 
                                    pcumulative_pv=min(manual$pcumulative_pv),
                                    tdollar=totaldollar, 
                                    scenario="optimization: total cases"))
  
  # factoring the scenario variable so that legend is in the correct order
  output <- manual %>%  full_join(optimizecases) %>% 
    mutate(scenario = factor(scenario, levels=c(
      'optimization: total cases',
      'high burden to high impact',
      'proportional allocation',
      'shrink the map'),
      labels = c('Optimized for case reduction',
                 'Prioritize high-transmission\nsettings',
                 'Proportional allocation',
                 'Prioritize low-transmission\nsettings'
      )))
  
  return(output)
  
}


## Figure 1. Global distribution of cases and population at risk ----------------
options(pillar.sigfigs=100)
sims_pf_pv <- readRDS(file.path(datadir, "dmm_output_foi_bednets_wpop_pf_pv.rds")) %>%
  filter(itn_cov==0) %>%
  mutate(init_EIR = round(init_EIR_pf+init_EIR_pv,6),
         cases = cases_pf+cases_pv,
         popatrisk_pf = ifelse(cases_pf<1,0,pop_pf),
         popatrisk_pv = ifelse(cases_pv<1,0,pop_pv),
         popatrisk=pmax(popatrisk_pf, popatrisk_pv)) %>%
  mutate(eir_groups_who = case_when(init_EIR < 0.1 ~ 'Very low',
                                    init_EIR >= 0.1 & init_EIR < 1 ~ 'Low',
                                    init_EIR >= 1 & init_EIR < 7 ~ 'Moderate',
                                    init_EIR >= 7 ~ 'High'),
         eir_groups_who = factor(eir_groups_who, levels = 
                                   c("Very low", "Low", "Moderate", "High"))) %>%
  # Note WHO classification of transmission intensities was derived by 
  # matching sum of Pf+Pv prev to sum of Pf+Pv EIR (and minor tweak for low/moderate grouping)
  # Prevalence grouping in https://www.who.int/publications/i/item/guidelines-for-malaria 
  mutate(pf_pv_group = case_when(init_EIR_pf == 0 & init_EIR_pv != 0 ~ "Pv",
                                 init_EIR_pf != 0 & init_EIR_pv == 0 ~ "Pf",
                                 init_EIR_pf != 0 & init_EIR_pv != 0 ~ "Both"),
         pf_pv_group = factor(pf_pv_group, levels = 
                                c("Pf", "Both", "Pv"))) 

sims_pf_pv <- group_by(sims_pf_pv, eir_groups_who, pf_pv_group) %>%
  summarise(cases=sum(cases),
            popatrisk_pf = sum(popatrisk_pf),
            popatrisk_pv = sum(popatrisk_pv),
            popatrisk=pmax(popatrisk_pf, popatrisk_pv),
            n = sum(n)) 

pf_pv_group_names <- c(
  "Pf" = "Settings endemic for\nP. falciparum only",
  "Both" = "Settings co-endemic for\nP. falciparum and P. vivax",
  "Pv" = "Settings endemic for\nP. vivax only"
)

plot_cases <- ggplot(data.frame(sims_pf_pv)) +
  geom_col(aes(x=eir_groups_who, y = cases/1000000), 
           colour="black", size = 1) +
  geom_col(aes(x=eir_groups_who, y = cases/1000000, fill = pf_pv_group)) +
  facet_grid(~pf_pv_group, labeller = as_labeller(pf_pv_group_names)) +
  labs(x = 'Transmission intensity', y = 'Clinical malaria cases\nper year (millions)') + 
  scale_fill_manual(values = pf_pv_group_colours) +
  theme_classic() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(-0.5,1, 0, 0), "lines"))

plot_pop <- ggplot(sims_pf_pv) +
  geom_col(aes(x=eir_groups_who, y = popatrisk/1000000000), colour="black", size =1) +
  geom_col(aes(x=eir_groups_who, y = popatrisk/1000000000, fill = pf_pv_group)) +
  facet_grid(~pf_pv_group, labeller = as_labeller(pf_pv_group_names)) +
  labs(x = 'Transmission intensity', y = 'Population at risk\nof malaria (billions)') + 
  scale_fill_manual(values = pf_pv_group_colours) +
  theme_classic() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = unit(c(-0.5,1, 0, 0), "lines"))

# Need to add 0s for number of countries
sims_pf_pv2 <- select(sims_pf_pv, eir_groups_who, pf_pv_group, n) %>%
  ungroup() %>%
  add_row(eir_groups_who = c("Moderate", "High", "High"), 
          pf_pv_group = c("Pv", "Both", "Pv"),
          n = c(0,0,0)) %>%
  mutate(eir_groups_who = factor(eir_groups_who, levels = 
                                   c("Very low", "Low", "Moderate", "High")),
         pf_pv_group = factor(pf_pv_group, levels = 
                                c("Pf", "Both", "Pv")))

plot_n <- ggplot(sims_pf_pv2, aes(x = eir_groups_who, y = 1, colour=pf_pv_group,
                                  label = format(n))) +
  geom_text(size = 3) +
  facet_grid(~pf_pv_group, labeller = as_labeller(pf_pv_group_names),
             scales='free', space='free') +
  scale_colour_manual(values = pf_pv_group_colours) +
  labs(x=NULL, y="Countries", color=NULL) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.ticks =  element_blank(),
        strip.background = NULL, 
        strip.text = element_blank(),
        plot.margin = unit(c(-0.5,1, 0, 0), "lines")) 

(plot_cases + plot_pop + plot_n + plot_layout(nrow=3, heights=c(3,3, 1)) + 
    plot_annotation(tag_levels = list(c("A","B", ""))))

#ggsave(paste0(outdir,'case_pop_distribution.png'), width = 20, height = 12.5, 
#       dpi = 300, units = 'cm')

# Proportion of pop at risk in very low+low coendemic settings
pmax(sum(filter(sims_pf_pv, pf_pv_group =="Both" & (eir_groups_who == "Low" | 
                                                      eir_groups_who == "Very low"))$popatrisk_pf),
     sum(filter(sims_pf_pv, pf_pv_group =="Both" & (eir_groups_who == "Low" | 
                                                      eir_groups_who == "Very low"))$popatrisk_pv))/
  pmax(sum(sims_pf_pv$popatrisk_pf),sum(sims_pf_pv$popatrisk_pf))

# Global distribution by country -----------------------------------------------
country_match_both <- readRDS(paste0(datadir, "country_match_both.rds"))

options(pillar.sigfigs=100)
country_match_both <- country_match_both %>%
  mutate(init_EIR = round(init_EIR_pf+init_EIR_pv,6)) %>%
  mutate(eir_groups_who = case_when(init_EIR < 0.1 ~ 'Very low',
                                    init_EIR >= 0.1 & init_EIR < 1 ~ 'Low',
                                    init_EIR >= 1 & init_EIR < 7 ~ 'Moderate',
                                    init_EIR >= 7 ~ 'High'),
         eir_groups_who = factor(eir_groups_who, levels = 
                                   c("Very low", "Low", "Moderate", "High"))) %>%
  mutate(pf_pv_group = case_when(init_EIR_pf == 0 & init_EIR_pv != 0 ~ "Pv",
                                 init_EIR_pf != 0 & init_EIR_pv == 0 ~ "Pf",
                                 init_EIR_pf != 0 & init_EIR_pv != 0 ~ "Both"),
         pf_pv_group = factor(pf_pv_group, levels = 
                                c("Pf", "Both", "Pv"))) 

# How many are African countries in each setting?
country_match_both %>% 
  mutate(region2 = ifelse(region=="Africa", "Africa", "Not Africa")) %>%
  group_by(eir_groups_who, pf_pv_group, region2) %>%
  summarise(count = n_distinct(ISO)) 

## Figure 2. Impact of ITN usage by EIR -----------------------------------------
# Changing A) annual clinical incidence (all ages), B) annual prevalence (2-10 years), C) annual clinical incidence (0-5 years), and D) number of adult mosquitoes per person, across various entomological inoculation rates (EIR). Lines are colored to represent percent reduction in mosquito to human force of infection (FOI). 

# plotting function
eir_plot2 <- function(var,text, tag, breaks){
  ggplot(output) + 
    geom_line(aes(x=itn_cov*100, y={{var}}*100, group=init_EIR, color=init_EIR)) +  # use x=diff for EIR diff
    labs(x="ITN usage (%)",y=text, color="Entomological\ninoculation\nrate",
         tag = tag) +                  # use x="% of initial EIR" for EIR diff
    scale_color_distiller(palette="Spectral", direction = -1, breaks = breaks) +
    theme_bw() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          plot.tag = element_text(), legend.title = element_text(size=10))
}

# Pf 
EIRdist <- readRDS(file.path(datadir, "dmm_output_foi_bednets_wpop_pf_pv.rds")) %>% 
  filter(init_EIR_pf !=0) %>%
  ungroup() %>% select(init_EIR_pf) %>% distinct() %>% unlist(use.names = F)

output <- readRDS(file.path(datadir, "dmm_output_foi_bednets_wpop_pf_pv.rds")) %>%
  filter(init_EIR_pf !=0) %>%
  rename(init_EIR = init_EIR_pf,
         clinical_incidence = clinical_incidence_pf,
         pop = pop_pf,
         prev = prev_pf) %>%
  filter(init_EIR %in% EIRdist[])

# plotting outcomes in patchwork
A <- eir_plot2(clinical_incidence, 
               text = "Clinical incidence\nper 100 person-years", tag = "A",
               breaks = waiver())
B <- eir_plot2(prev, "\nPrevalence 2-10 years (%)", tag = "B",
               breaks = waiver())

# Pv
EIRdist <- readRDS(file.path(datadir, "dmm_output_foi_bednets_wpop_pf_pv.rds")) %>% 
  filter(init_EIR_pv !=0) %>%
  ungroup() %>% select(init_EIR_pv) %>% distinct() %>% unlist(use.names = F)

output <- readRDS(file.path(datadir, "dmm_output_foi_bednets_wpop_pf_pv.rds")) %>%
  filter(init_EIR_pv !=0) %>%
  rename(init_EIR = init_EIR_pv,
         clinical_incidence = clinical_incidence_pv,
         pop = pop_pv,
         prev = prev_pv) %>%
  filter(init_EIR %in% EIRdist[])

# plotting outcomes in patchwork
D <- eir_plot2(clinical_incidence, "Clinical incidence\nper 100 person-years",
               tag="C", breaks = waiver())
E <- eir_plot2(prev, "\nPrevalence 0-99 years (%)", tag = "D", breaks = waiver())

# Create final plot with proper alignment of legend
pf_title <- text_grob(expression(italic('P. falciparum')),size = 12)
pf_title <- as_ggplot(pf_title) + theme(plot.margin = margin(0,0,0,-6, "cm"))   #0.25,11
pv_title <- text_grob(expression(italic('P. vivax')),size = 12)
pv_title <- as_ggplot(pv_title) + theme(plot.margin = margin(0,0,0,-7, "cm"))

X <- ggarrange(pf_title, NULL, A,B, 
               ncol=2, nrow=2, common.legend = TRUE, legend="right",
               heights = c(1,5))

Y <- ggarrange(pv_title, NULL, D, E, 
               ncol=2, nrow=2, common.legend = TRUE, legend="right",
               heights = c(1,5))
ggarrange(X,Y, nrow=2)

#ggsave(paste0(outdir,'bednet_outcomes_all.pdf'), width = 20, height = 12, 
#       dpi = 300, units = 'cm')

# Summary statistics -----------------------------------------------------------
sims_pf_pv <- readRDS(file.path(datadir, "dmm_output_foi_bednets_wpop_pf_pv.rds")) %>%
  filter(itn_cov <=0.8) %>%
  select(init_EIR_pf, init_EIR_pv, itn_cov, clinical_incidence_pf, clinical_incidence_pv, 
         cases_pf, cases_pv)
sims_pf_pv_baseline <- filter(sims_pf_pv, itn_cov==0) %>%
  rename(binc_pf = clinical_incidence_pf,
         binc_pv = clinical_incidence_pv,
         bcases_pf = cases_pf, 
         bcases_pv = cases_pv)

sims_pf_pv <- left_join(sims_pf_pv, sims_pf_pv_baseline, by = c("init_EIR_pf",
                                                                "init_EIR_pv")) %>%
  mutate(rel_red_inc_pf = (binc_pf-clinical_incidence_pf)/binc_pf,
         rel_red_inc_pv = (binc_pv-clinical_incidence_pv)/binc_pv)

# Plot relative reduction in clinical incidence
ggplot(filter(sims_pf_pv, itn_cov.x %in% c(0.25,0.5,0.8))) +
  geom_point(aes(x=init_EIR_pf, y = rel_red_inc_pf)) +
  geom_point(aes(x=init_EIR_pv, y = rel_red_inc_pv), col="red") +
  facet_wrap(~itn_cov.x)

# Cases at baseline
sum(filter(sims_pf_pv, itn_cov.x == 0)$cases_pf)
sum(filter(sims_pf_pv, itn_cov.x == 0)$cases_pv)

# Case reductions with 80% usage
(sum(filter(sims_pf_pv, itn_cov.x == 0)$cases_pf)-
    sum(filter(sims_pf_pv, itn_cov.x == 0.8)$cases_pf))/
  sum(filter(sims_pf_pv, itn_cov.x == 0)$cases_pf)

(sum(filter(sims_pf_pv, itn_cov.x == 0)$cases_pv)-
    sum(filter(sims_pf_pv, itn_cov.x == 0.8)$cases_pv))/
  sum(filter(sims_pf_pv, itn_cov.x == 0)$cases_pv)

## Figure 3. Impact of strategies on cases and population at risk (Pf+Pv)---------

manual <- readRDS(paste0(outdir, "bednets_data_equilibrium_pf_pv_80_loess.rds"))
optimizecases <- readRDS(paste0(outdir, 'SA_63_pf_pv_80_loess.rds'))

output <- combine_strategies_impact(manual = manual,
                                    optimizecases = optimizecases)

maxx <- max(output$tdollar)/1000000000

colors <- c('#F8766D', '#00B9E3','#619CFF','#00BA38','#D39200')

CP_A <- ggplot(output, aes(x=tdollar/1000000000,y=tcumulative/1000000,color=scenario)) + 
  geom_line(size=1.1, alpha=0.5, show.legend = F) +  
  labs(x="Budget (billion $)",
       y="Total malaria cases\n(millions)",
       color='Allocation strategy', tag = "A") +
  scale_color_manual(guide = guide_legend(reverse = F), values = colors) +
  scale_x_continuous(labels = scales::comma, limits=c(0, maxx), breaks=seq(0,maxx,1)) + 
  theme_classic() 

CP_B <- ggplot(output, aes(x=tdollar/1000000000,y=pcumulative/1000000000,color=scenario)) + 
  geom_step(size=1.1, alpha=0.5) +
  labs(x="Budget (billion $)",
       y="Total population at risk of malaria\n(billions)",
       color='Allocation strategy', tag = "B") +
  scale_color_manual(guide = guide_legend(reverse = F), values = colors) +
  scale_x_continuous(labels = scales::comma, limits=c(0, maxx), breaks=seq(0,maxx,1)) +    
  theme_classic() 

x <- (CP_A + CP_B) + plot_layout(nrow=1) 
x
# ggsave(paste0(outdir,'model_compare_pf_pv_80_loess.png'), 
#        x, width = 25, height = 10, dpi=300, 
#        units = 'cm')

# Plot cases by species
CP_C <- ggplot(output) + 
  geom_line(aes(x=tdollar/1000000000,y=tcumulative_pf/1000000,color="P. falciparum"),
            size=1.1) +
  geom_line(aes(x=tdollar/1000000000,y=tcumulative_pv/1000000,color="P. vivax"),
            size=1.1) +
  facet_wrap(~scenario, nrow=1) +
  labs(x="Budget (billion $)",
       y="Species-specific\ncases (millions)",
       color='Parasite species', tag = "C") +
  scale_colour_manual(values = pf_pv_group_colours[c(1,3)],
                      labels = c("P. falciparum" = expression(italic('P. falciparum')),
                                 "P. vivax" = expression(italic('P. vivax')))) +
  scale_x_continuous(labels = scales::comma, limits=c(0, maxx), breaks=seq(0,maxx,1)) + 
  theme_classic() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        legend.title = element_blank(),
        legend.position = c(0.18, 0.76))

CP_D <- ggplot(output) + 
  geom_step(aes(x=tdollar/1000000000,y=pcumulative_pv/1000000000,color="P. vivax"),
            size=1.1) +
  geom_step(aes(x=tdollar/1000000000,y=pcumulative_pf/1000000000,color="P. falciparum"),
            size=1.1) +
  facet_wrap(~scenario, nrow=1) +
  labs(x="Budget (billion $)",
       y="Species-specific\npopulation at risk (billions)",
       color='Parasite species', tag = "D") +
  scale_colour_manual(values = pf_pv_group_colours[c(1,3)],                      
                      labels = c("P. falciparum" = expression(italic('P. falciparum')),
                                 "P. vivax" = expression(italic('P. vivax')))) +
  scale_y_continuous(breaks = c(0:4), labels = c("0", "1", "2", "3",
                                                 "    4")) +
  scale_x_continuous(labels = scales::comma, limits=c(0, maxx), breaks=seq(0,maxx,1)) + 
  theme_classic() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(), 
        legend.position = "none")

(CP_C + CP_D) + plot_layout(nrow=2, guides = "collect") 

#ggsave(paste0(outdir,'model_compare_pf_pv_80_loess_by_species.png'), width = 25, height = 15, units = 'cm')

# Combine both
CP_X <- ggarrange(CP_A,CP_B, 
                  ncol=2, nrow=1, common.legend = TRUE, legend="right")

CP_Y <- ggarrange(CP_C, CP_D,
                  ncol=1, nrow=2)
ggarrange(CP_X,CP_Y, nrow=2, heights=c(3,4))

# ggsave(paste0(outdir,'model_compare_pf_pv_80_loess_full.pdf'), width = 25, 
#        height = 20, dpi=300, units = 'cm')

# Summary statistics on impact of different strategies -------------------------

# Plot proportion of cases attributable to P. vivax
ggplot(output) + 
  geom_line(aes(x=tdollar/1000000000,y=tcumulative_pv/tcumulative),
            size=1.1) +
  facet_wrap(~scenario, nrow=1) +
  labs(x="Budget (billion $)",
       y="Proportion of cases attributable to vivax") +
  scale_x_continuous(labels = scales::comma, limits=c(0, maxx), breaks=seq(0,maxx,250)) + 
  theme_classic() 

# Table of these results

# Extract statistics (chose 1 budget_levels)
budget_levels <- c(0.25,0.5,0.75,1)*max(output$tdollar)
budget_levels <- seq(0,1,0.05)*max(output$tdollar)

# Interpolate values for required budget levels
extractoutput <- list()
for (i in 1:4) {
  extractoutput[[i]] <-
    data.frame(scenario = unique(output$scenario)[i],
               budget=budget_levels,
               cases = 
                 approx(x=filter(output, scenario == unique(output$scenario)[i])$tdollar,
                        y=filter(output, scenario == unique(output$scenario)[i])$tcumulative,
                        xout =  budget_levels)$y,
               popatrisk = 
                 approx(x=filter(output, scenario == unique(output$scenario)[i])$tdollar,
                        y=filter(output, scenario == unique(output$scenario)[i])$pcumulative,
                        xout =  budget_levels)$y)
}
extractoutput <- do.call("rbind", extractoutput)    # checked it aligns with plot

totalcases <-  manual %>% filter(scenario == 'high burden to high impact') %>% 
  filter(itn_cov==0) %>% 
  summarize(cases=sum(cases, na.rm=T)) %>% as.numeric()

totalpop_pf <-  manual %>% filter(scenario == 'high burden to high impact') %>% 
  filter(itn_cov==0) %>% summarize(tpop=sum(pop_pf, na.rm=T)) %>% as.numeric()
totalpop_pv <-  manual %>% filter(scenario == 'high burden to high impact') %>% 
  filter(itn_cov==0) %>% summarize(tpop=sum(pop_pv, na.rm=T)) %>% as.numeric()
totalpop <- pmax(totalpop_pf, totalpop_pv)

extractoutput_table <- mutate(extractoutput,
                              cases_abs_red = totalcases-cases,
                              cases_millions = round(cases/1000000,1), 
                              cases_rel_red = round(cases_abs_red/totalcases*100,0),
                              popatrisk_billions = round(popatrisk/1000000000,1),
                              pop_abs_red = totalpop-popatrisk,
                              pop_rel_red = round(pop_abs_red/totalpop*100,0)) %>%
  arrange(scenario, budget) %>%
  select(-cases_abs_red, - pop_abs_red, - cases, -popatrisk)
#write.csv(extractoutput_table, paste0(outdir, "impact_table.csv"))

# Compare cases to optimiser
optim_comp <- select(extractoutput, budget, scenario, cases) %>%
  pivot_wider(id_cols = "budget",
              names_from = "scenario", values_from ="cases") 
colnames(optim_comp) <- c("budget", "high_transmission", "low_transmission", "proportional",
                          "optim")
optim_comp <- optim_comp %>%
  mutate(high_transmission_diff = round(high_transmission-optim,0),
         high_transmission_ratio = high_transmission/optim,
         low_transmission_diff = round(low_transmission-optim,0),
         low_transmission_ratio = low_transmission/optim,
         proportional_diff = round(proportional-optim,0),
         proportional_ratio =proportional/optim,
         percent_lower_than_high_transmission = 1-1/high_transmission_ratio,
         percent_lower_than_low_transmission = 1-1/low_transmission_ratio,
         percent_lower_than_proportional = 1-1/proportional_ratio)

# Maximum difference between optimiser and each strategy
optim_comp %>%
  filter(percent_lower_than_high_transmission==max(percent_lower_than_high_transmission)) %>%
  select(budget, percent_lower_than_high_transmission,high_transmission_diff) %>%
  mutate(budget_prop = budget/max(output$tdollar))

optim_comp %>%
  filter(percent_lower_than_proportional==max(percent_lower_than_proportional)) %>%
  select(budget, percent_lower_than_proportional, proportional_diff) %>%
  mutate(budget_prop = budget/max(output$tdollar))

optim_comp %>%
  filter(percent_lower_than_low_transmission==max(percent_lower_than_low_transmission)) %>%
  select(budget, percent_lower_than_low_transmission, low_transmission_diff) %>%
  mutate(budget_prop = budget/max(output$tdollar))

# Compare relative reductions to optimiser
optim_comp2 <- select(extractoutput_table, budget, scenario, cases_rel_red) %>%
  pivot_wider(id_cols = "budget",
              names_from = "scenario", values_from ="cases_rel_red") 
colnames(optim_comp2) <- c("budget", "optim", "high_transmission", 
                           "proportional", "low_transmission")
optim_comp2 <- optim_comp2 %>%
  mutate(high_transmission_rel_diff = (optim-high_transmission)/optim,
         low_transmission_rel_diff = (optim-low_transmission)/optim,
         proportional_rel_diff = (optim-proportional)/optim,
         budget_prop = budget/max(output$tdollar)) 
View(filter(optim_comp2, budget_prop %in% c(0.25,0.5,0.75)))

# At 50% budget, how many settings have eliminated in shrink the map?
closest_budget_50 <- manual$tdollar[manual$scenario == "shrink the map"]
closest_budget_50 <- closest_budget_50[which(abs(closest_budget_50-0.5*max(manual$tdollar))==
                                               min(abs(closest_budget_50-0.5*max(manual$tdollar))))]

eir_selection <- manual %>% filter(scenario == "shrink the map" & tdollar <= closest_budget_50) %>%
  ungroup() %>%
  mutate(B=max(tdollar)) %>%
  group_by(init_EIR) %>%
  filter(itn_cov==max(itn_cov)) %>% 
  full_join(filter(manual, scenario == "shrink the map" & itn_cov==0)) %>%
  group_by(init_EIR) %>%
  filter(itn_cov==max(itn_cov)) %>%
  mutate(eliminated = ifelse(popatrisk==0, "yes", "no"))

group_by(eir_selection, eliminated) %>% summarise(n())

## Figure 4. Compare allocation pattern for Pf Pv case optimiser -----------------
case_optimiser <- readRDS(paste0(outdir, "SA_63_pf_pv_80_loess.rds"))
n_eir_groups <- 5

options(pillar.sigfigs=100)
case_optimiser <- case_optimiser %>%
  mutate(scenario = "optimization: total cases",
         prop=b/B,
         final_cov=usage,
         init_EIR = round(init_EIR_pf+init_EIR_pv,6)) %>%
  replace_na(list(prop=0)) %>%
  select(init_EIR_pf, init_EIR_pv,init_EIR, scenario, B, b, prop, final_cov, cases_pf, cases_pv, cases)

allocation_df <- case_optimiser %>%
  mutate(quartile_eir = ntile(init_EIR, n_eir_groups)) %>%
  group_by(quartile_eir) %>%
  mutate(min_eir = min(init_EIR),
         eir_range = paste0(min(init_EIR),"-",max(init_EIR))) %>%
  ungroup() %>%
  mutate(eir_groups_who = case_when(init_EIR < 0.1 ~ 'Very low',
                                    init_EIR >= 0.1 & init_EIR < 1 ~ 'Low',
                                    init_EIR >= 1 & init_EIR < 7 ~ 'Moderate',
                                    init_EIR >= 7 ~ 'High'),
         eir_groups_who = factor(eir_groups_who, levels = 
                                   c("Very low", "Low", "Moderate", "High"))) %>%
  mutate(pf_pv_group = case_when(init_EIR_pf == 0 & init_EIR_pv != 0 ~ "Pv",
                                 init_EIR_pf != 0 & init_EIR_pv == 0 ~ "Pf",
                                 init_EIR_pf != 0 & init_EIR_pv != 0 ~ "Both")) %>%
  mutate(prop_pf = cases_pf/cases,
         prop_pf_groups = case_when(prop_pf <0.5 ~ "Majority Pv cases",
                                    prop_pf >=0.5 ~ "Majority Pf cases"),
         prop_pf_groups = factor(prop_pf_groups, levels = 
                                   c("Majority Pv cases", "Majority Pf cases"))) %>%
  select(init_EIR_pf, init_EIR_pv, init_EIR, scenario, eir_groups_who, pf_pv_group, B, b, prop, final_cov)

# Add budget 0
allocation_df <- rbind(allocation_df, 
                       data.frame(distinct(allocation_df, init_EIR_pf, init_EIR_pv, init_EIR, scenario, 
                    eir_groups_who, pf_pv_group)) %>%
  mutate(B=0.000001, b=0, prop = 0, final_cov = 0))  # set B to very small value for division

# Proportion of budget by total WHO transmission intensity
p1 <- ggplot(allocation_df %>% group_by(B, eir_groups_who) %>% summarise(b=sum(b))) +
  geom_area(aes(y=b/B*100, x=B/1000000000, group=eir_groups_who, fill=eir_groups_who), size=1) +
  scale_fill_manual(values=transmission_intensity_group_colours) +  
  labs(y= "Proportion of total budget (%)", x= "Budget (billion $)", 
       fill="Transmission intensity") +
  theme_classic() + 
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 

# ITN usage by total WHO transmission intensity
p2 <- ggplot(allocation_df %>% group_by(B, eir_groups_who) %>% summarise(usage=mean(final_cov))) +
  geom_line(aes(y=usage*100, x=B/1000000000, group=eir_groups_who, colour=eir_groups_who), size=1) +
  scale_colour_manual(values=transmission_intensity_group_colours) + 
  labs(y= "Funded mean ITN usage (%)", x= "Budget (billion $)", 
       colour="Transmission intensity") +
  theme_classic() +
  ylim(0,80) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = "none") 

# Proportion of budget by settings endemic for Pf only, Pv only or both
p3 <- ggplot(allocation_df %>% group_by(B, pf_pv_group) %>% summarise(b=sum(b))) +
  geom_area(aes(y=b/B*100, x=B/1000000000, group=pf_pv_group, fill=pf_pv_group), size=1) +
  scale_fill_manual(values = pf_pv_group_colours[c(2,1,3)],
                    labels = c("Pf" = expression(paste(italic('P. falciparum'), " only")),
                               "Pv" = expression(paste(italic('P. vivax'), " only")),
                               "Both" = expression(paste(italic('P. falciparum'), " and ",
                                                         italic('P. vivax'))))) +
  labs(y= "Proportion of total budget (%)", x= "Budget (billion $)", 
       fill="Endemic species") +
  theme_classic() +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) 

# ITN usage by settings endemic for Pf only, Pv only or both
p4 <- ggplot(allocation_df %>% group_by(B, pf_pv_group) %>% summarise(usage=mean(final_cov))) +
  geom_line(aes(y=usage*100, x=B/1000000000, group=pf_pv_group, colour=pf_pv_group), size=1) +
  scale_colour_manual(values = pf_pv_group_colours[c(2,1,3)],
                      labels = c("Pf" = expression(paste(italic('P. falciparum'), " only")),
                                 "Pv" = expression(paste(italic('P. vivax'), " only")),
                                 "Both" = expression(paste(italic('P. falciparum'), " and ",
                                                           italic('P. vivax'))))) +
  labs(y= "Funded mean ITN usage (%)", x= "Budget (billion $)", 
       colour="Endemic species") +
  theme_classic() +
  ylim(0,80) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position="none") 

(p1 + p2) / (p3 + p4) + plot_layout(guides = "collect", nrow=2) + 
  plot_annotation(tag_levels = c("A", "B", "C", "D")) 

# ggsave(paste0(outdir,'allocation_pattern_optimal_pf_pv.pdf'), 
#        width = 25, height = 15, dpi=300, units = 'cm')

# Summary statistics on allocation pattern for Pf Pv optimiser ----
df <- allocation_df %>% 
  group_by(B, eir_groups_who, pf_pv_group) %>% 
  summarise(b=sum(b), 
            mean_usage=mean(final_cov)) %>%
  mutate(prop = b/B)

# Proportion of budget allocated to high burden P. falciparum only settings (i.e. Africa)
budget_prop <- seq(0,1,0.05)
prop_allocated_high_pf <- approxfun(x=filter(df, eir_groups_who == "High" & pf_pv_group == "Pf")$B,
                                    y=filter(df, eir_groups_who == "High" & pf_pv_group == "Pf")$prop, method = "linear")
prop_allocated_moderate <- approxfun(x=filter(df, eir_groups_who == "Moderate")$B,
                                     y=filter(df, eir_groups_who == "Moderate")$prop, method = "linear")
prop_allocated_low_pf_pv <- approxfun(x=filter(df, eir_groups_who == "Low" & pf_pv_group == "Both")$B,
                                      y=filter(df, eir_groups_who == "Low" & pf_pv_group == "Both")$prop, method = "linear")
prop_allocated_very_low_pf_pv <- approxfun(x=filter(df, eir_groups_who == "Very low" & pf_pv_group == "Both")$B,
                                           y=filter(df, eir_groups_who == "Very low" & pf_pv_group == "Both")$prop, method = "linear")
# All high transmission intensity is P. falciparum only so that proportion is the same
prop_allocated_summary <- data.frame(budget_prop=budget_prop,
                                     prop_for_high_transmission_pf = 
                                       prop_allocated_high_pf(budget_prop*max(df$B)),
                                     prop_for_moderate_transmission =
                                       prop_allocated_moderate(budget_prop*max(df$B)),
                                     prop_for_low_transmission_pf_pv =
                                       prop_allocated_low_pf_pv(budget_prop*max(df$B)),
                                     prop_for_very_low_transmission_pf_pv =
                                       prop_allocated_very_low_pf_pv(budget_prop*max(df$B))) %>%
  mutate(sum_high_moderate = prop_for_high_transmission_pf + prop_for_moderate_transmission,
         sum_low_very_low_pf_pv = prop_for_low_transmission_pf_pv + prop_for_very_low_transmission_pf_pv)


## Figure S2. Model simulations over time-----------------------------------------
# Pf
output_pf <- readRDS(paste0(datadir, "dmm_output_foi_bednets_yearly_pf.rds"))

# choose EIRS to plot
readRDS(paste0(datadir, "dmm_output_foi_bednets_wpop_pf.rds")) %>%
  ungroup() %>%
  filter(itn_cov==0) %>% 
  select(init_EIR, pop) %>%
  arrange(init_EIR) %>% 
  mutate(t = sum(pop), ct = cumsum(pop), p = ct / t * 100)
ITNset <- seq(0,1,0.01)
EIRset <- c(0.001, 0.05, 140)

A <- ggplot(data = output_pf %>%
              filter(itn_cov %in% ITNset & init_EIR %in% EIRset) %>% 
              rename(EIR = init_EIR)) +
  geom_line(aes(x=years, y=clinical_incidence*100, color=-itn_cov*100, group=-itn_cov*100)) + 
  labs(x='Year', y="Clinical incidence per 100 person-years", color="ITN usage (%)") + 
  guides(color = guide_legend(reverse=TRUE)) +
  facet_wrap(. ~ EIR, scales = 'free_y', labeller = label_both) + 
  scale_color_distiller(palette="YlGnBu", trans ="reverse",
                        labels = c(100,75,50,25,0)) +  
  theme_classic()

# vivax
# read in yearly bed net runs
output_pv <- readRDS(paste0(datadir, "dmm_output_foi_bednets_yearly_pv.rds"))

# choose EIRS to plot
readRDS(paste0(datadir, "dmm_output_foi_bednets_wpop_pv.rds")) %>%
  ungroup() %>%
  filter(itn_cov==0) %>% 
  select(init_EIR, pop) %>%
  arrange(init_EIR) %>% 
  mutate(t = sum(pop), ct = cumsum(pop), p = ct / t * 100)

# choose set to plot - just EIRs 0.1 and 80 right now
ITNset <- seq(0,1,0.01)
EIRset <- c(0.001, 0.1, 1.7)

B <- ggplot(data = output_pv %>%
              filter(itn_cov %in% ITNset & init_EIR %in% EIRset) %>% 
              rename(EIR = init_EIR)) +
  geom_line(aes(x=years, y=clinical_incidence*100, color=-itn_cov*100, group=-itn_cov*100)) + 
  labs(x='Year', y="Clinical incidence per 100 person-years", color="ITN usage (%)") + 
  guides(color = guide_legend(reverse=TRUE)) +
  facet_wrap(. ~ EIR, scales = 'free_y', labeller = label_both) + 
  scale_color_distiller(palette="YlGnBu", trans ="reverse",
                        labels = c(100,75,50,25,0)) +  
  theme_classic()

A / B + plot_layout(nrow=2, guides="collect") + plot_annotation(tag_levels = 'A')

#ggsave(paste0(outdir,'odin_equilibrium.png'), width = 25, height = 15, units = 'cm')

## Figure S3 is in nets_conversion script ---------------------------------------
## Figure S4. Surface plots ----------------------------------------------------
surfaceplot <- function(species, ITNuse, outcome, titlename){
  
  # pull in main simulation runs for either species
  if (species == 'pf') {
    file <- readRDS(paste0(datadir, 'dmm_output_foi_bednets_wpop_pf.rds'))
  }
  
  if (species == 'pv') {
    file <- readRDS(paste0(datadir, 'dmm_output_foi_bednets_wpop_pv.rds'))
  }
  
  # set max ITN coverage
  file <- file %>% filter(itn_cov <= ITNuse/100)
  
  if (ITNuse == 80) {
    netz <- readRDS(paste0(datadir, "netz_median_use_rate.rds"))
  }
  
  if (ITNuse == 90) {
    netz <- readRDS(paste0(datadir, "netz_max_use_rate.rds"))
  }
  
  file <- file %>%
    left_join(netz, by=c('itn_cov'='usage'))
  
  # create surface
  dimx <- 100
  dimy <- 100

  x <- file$itn_cov
  
  y <- file$init_EIR
  
  if (outcome == 1) {
    file <- file %>% filter(!is.na(clinical_incidence_bednet))
    z <- file$clinical_incidence_bednet
    
  }
  
  if (outcome == 2) {
    file <- file %>% filter(!is.na(prev))
    z <- file$prev
  }
  
  spline <- akima::interp(x,y,z, yo=seq(min(y), max(y),length.out=dimy), xo=seq(0,1,length.out=dimx), duplicate = 'mean', linear = TRUE)
  
  p <- plot_ly() %>%
    add_trace(x=spline$x, y=spline$y, z=t(spline$z), type="surface") %>%
    layout(
      showlegend=T,
      scene = list(
        xaxis = list(title = "ITN usage", range=c(0,ITNuse/100)),
        yaxis = list(title = "EIR"),
        zaxis = list(title = titlename)
      )) %>% hide_colorbar()
  
  p
  
}

# call function and plot
ppf <- surfaceplot("pf", 80,  1, "Clinical incidence (all ages, per year)"); ppf
ppv <- surfaceplot("pv", 80, 1, "Clinical incidence (all ages, per year)"); ppv

## Figure S5+S6. Prev to EIR --------------------------------------------------------
# country prevalence match to EIR

output_pf <- readRDS(file.path(datadir, "country_match_pf.rds")) %>%
  rename(prev = prev_pf)
output_pv <- readRDS(file.path(datadir, "country_match_pv.rds")) %>%
  rename(prev = prev_pv)


# creating pop + eir + PfPR2-10 matches
eirPR_match <- function(data, ylab, overlaps = 10, base_size = 10){
  
  match <- data %>%
    mutate(pop_admins0 = ifelse(pop_admins0==0,NA_real_,pop_admins0),
           # set population size groups
           pop_group = case_when(pop_admins0<1000000~1,
                                 pop_admins0>=1000000 &  pop_admins0 <= 10000000 ~ 2,
                                 pop_admins0>10000000 & pop_admins0 <= 100000000 ~ 3,
                                 pop_admins0>100000000 ~ 4)) %>%
    # change region for a few mismatches by hand
    mutate(region = case_when(Country=='Jordan'~'Eastern Mediterranean',
                              Country=='Uruguay'~'Americas',
                              Country=='Chile'~'Americas',
                              Country=='Mongolia'~'Western Pacific',
                              Country=='Russian Federation'~'Europe',
                              Country=='Taiwan, Province of China'~'Western Pacific',
                              Country=='United States'~'Americas',
                              TRUE~region)) %>%
    filter(!is.na(pop_admins0))
  
  # create plot matching countries to EIR and prevalence
  A <- ggplot(match) +
    geom_point(aes(x=init_EIR, y=prev,
                   color=region, size=pop_group), alpha=0.3) +
    labs(x="EIR", y=ylab, size='Population',color="WHO region") +
    geom_text_repel(aes(x=init_EIR,y=prev, label=ISO), size=2, box.padding = 0,
                    max.overlaps = overlaps) +
    scale_size_continuous(range = c(2, 5), labels = c("< 1M","1M to 10M", "10M to 100M", "> 100M")) +
    theme_classic(base_size=12) +
    theme(legend.position = c(.75,.4))
  
  # making table of pop in each EIR to go with plot
  poptable <- match %>% group_by(init_EIR) %>%
    summarize(n=n(),
              pop = round(sum(pop_admins0),0)) %>%
    mutate(per = round(pop / sum(pop)*100,2),
           pop = format(pop,big.mark=","))
  
  # plot and save
  A + tableGrob(poptable, rows = NULL, cols = c('EIR', 'n', 'population', '% of global'), 
                theme=ttheme_default(base_size=base_size))
  
}

# call function and save plots
eirPR_match(output_pf, "P. falciparum prevalence, 2000 \n(2-10 years)", base_size = 7)
#ggsave(paste0(outdir,'country_EIR_match_pf.pdf'), width = 22, height = 18, units = 'cm')

eirPR_match(output_pv, "P. vivax prevalence, 2000 \n(0-99 years)", overlaps = 50,
            base_size = 10)
#ggsave(paste0(outdir,'country_EIR_match_pv.pdf'), width = 22, height = 18, units = 'cm')

## Figure S7. Compare impact with loess vs linear -------------------------------
manual_loess <- readRDS(paste0(outdir, "bednets_data_equilibrium_pf_pv_80_loess.rds"))
manual_linear <- readRDS(paste0(outdir, "bednets_data_equilibrium_pf_pv_80_linear.rds"))
optimizecases_loess <- readRDS(paste0(outdir, 'SA_63_pf_pv_80_loess.rds'))
optimizecases_linear <- readRDS(paste0(outdir, 'SA_63_pf_pv_80_linear.rds'))

output_loess <- combine_strategies_impact(manual = manual_loess,
                                          optimizecases = optimizecases_loess)

output_linear <- combine_strategies_impact(manual = manual_linear,
                                           optimizecases = optimizecases_linear)

# Transform budgets into % of manual
output_loess <- mutate(output_loess, 
                       tdollar_prop = tdollar/max(manual_loess$tdollar),
                       assumption = "loess")
output_linear <- mutate(output_linear, 
                        tdollar_prop = tdollar/max(manual_linear$tdollar),
                        assumption = "linear")

output <- rbind(output_loess,output_linear) %>% 
  mutate(assumption = factor(assumption, levels = c("loess", "linear")))

colors <- c('#F8766D', '#00B9E3','#619CFF','#00BA38','#D39200')

CP_A <- ggplot(output) + 
  geom_line(aes(x=tdollar_prop*100,y=tcumulative/1000000,color=scenario, 
                linetype = assumption),
            size=1.1, alpha=0.5, show.legend = F) +   
  labs(x="Budget (% of maximum)",
       y="Total malaria cases\n(millions)",
       color='Allocation strategy', linetype = "ITN cost-usage relationship",
       tag = "A") +
  scale_color_manual(guide = guide_legend(reverse = F), values = colors) +
  scale_linetype_discrete(labels = c("Non-linear","Linear")) +
  guides(colour = guide_legend(order = 1), 
         linetype = guide_legend(order = 2)) +
  theme_classic() 

CP_B <- ggplot(output) + 
  geom_step(aes(x=tdollar_prop*100,y=pcumulative/1000000000,color=scenario,
                linetype = assumption),
            size=1.1, alpha=0.5) +
  labs(x="Budget (% of maximum)",
       y="Total population at risk of malaria\n(billions)",
       color='Allocation strategy', linetype = "ITN cost-usage relationship",
       tag = "B") +
  scale_color_manual(guide = guide_legend(reverse = F), values = colors) +
  scale_linetype_discrete(labels = c("Non-linear","Linear")) +
  guides(colour = guide_legend(order = 1), 
         linetype = guide_legend(order = 2)) +
  theme_classic() 

x <- (CP_A + CP_B) + plot_layout(nrow=1) 
x
# ggsave(paste0(outdir,'sensitivity_model_compare_pf_pv_80_loess_vs_linear.png'), 
#        x, width = 25, height = 10, dpi=300, 
#        units = 'cm')
