# Find annual nets to distribute for a given net usage 
library(here)
library(tidyverse)
library(ggplot2)

# working directory
datadir <- here("files", "/")

# Load and prepare data  -------------------------------------------------------

# Inputs

# 2019 data of access, use rate and npc
cube_nat_level_annual <- read.csv("https://raw.github.com/bertozzivill/map-itn-gts/master/data/coverage_metrics/aggregated_predictions_2019.csv")
cube_nat_level_annual  <- cube_nat_level_annual[cube_nat_level_annual$iso3 != "AFR" &
                                                  cube_nat_level_annual$year == 2019,]

# 2020 data of access, use rate and npc
cube_nat_level <- read.csv("https://raw.github.com/bertozzivill/map-itn-cube/publication-2021/paper_figures/figure_data/fig_4_access_npc.csv")

# 2020 estimates of net half lives (median retention times) 
half_life_data <- read.csv("https://raw.github.com/bertozzivill/map-itn-cube/publication-2021/paper_figures/figure_data/fig_5_llin_half_lives.csv")

# Data Prep

# Aggregate annual data
# From this, annual mean use rate is used to convert target usage to target access
cube_nat_level_annual <- cube_nat_level_annual[is.na(cube_nat_level_annual$month),]
cube_nat_level_annual <- reshape(cube_nat_level_annual[,c("iso3", "mean", "year", "variable")], 
                                 idvar = c("iso3", "year"), timevar = "variable", direction = "wide",
                                 varying=list(c(unique(cube_nat_level_annual$variable))))

# Data by month
# Monthly mean is used to generate loess curve of access against nets per capita
cube_nat_level <- cube_nat_level[,c("iso3", "year","month", "access_mean", "percapita_nets_mean")]
names(cube_nat_level)[names(cube_nat_level) %in% 
                        c("access_mean", "percapita_nets_mean")] <- c("access", "percapita_nets")

# Explore data

# Country-specific half life estimates
hist(half_life_data$half_life)
# Country specific use rate estimates
hist(cube_nat_level_annual$use_rate)

# Loess curve of access vs nets per capita -------------------------------------

# Fit Loess curve of access against NPC based on data from all countries and all months
curve_fit <- loess(access ~ percapita_nets, data=cube_nat_level,
                   control=loess.control(surface="direct")) 
loess_for_prediction <- cube_nat_level[,c("iso3", "month", "access", "percapita_nets")]
loess_for_prediction$loess_predicted_access <- predict(curve_fit)
loess_for_prediction <- loess_for_prediction[order(loess_for_prediction$loess_predicted_access),]

# Extrapolate higher access by extending Loess curve
extrapolate_access <- data.frame(iso3=NA, month=NA, access=NA, percapita_nets=
             seq(ceiling(max(cube_nat_level$percapita_nets)*100)/100,1,0.01))
extrapolate_access$loess_predicted_access <- predict(curve_fit, extrapolate_access$percapita_nets)

# Extrapolate lower access by assuming all access below minimum observed
# requires the same nets per capita (i.e. the same cost)
extrapolate_access2 <- data.frame(iso3=NA, month=NA, access=NA, percapita_nets=0,
                                  loess_predicted_access=seq(0,floor(min(loess_for_prediction$loess_predicted_access)*100)/100,0.01))
extrapolate_access2$percapita_nets[extrapolate_access2$loess_predicted_access != 0] <-
  loess_for_prediction$percapita_nets[loess_for_prediction$loess_predicted_access==min(loess_for_prediction$loess_predicted_access)]

# Combine into 1 dataset
loess_for_prediction <- rbind(extrapolate_access2,loess_for_prediction, extrapolate_access)

# Plot access vs nets per capita
ggplot(loess_for_prediction) +
  geom_point(aes(x=percapita_nets, y = access)) +
  geom_line(aes(x=percapita_nets, y = loess_predicted_access), col = "red", size=1.1) +
  xlim(0,1) + ylim(0,1) +
  theme_classic()

# Convert net usage to nets distributed per person-year ------------------------

# Inputs/assumptions

# Choose a target usage (all values in simulations)
target_usage <- seq(0,1,0.01)

# Choose a net half life in days
half_life <- median(half_life_data$half_life)*365

# Create function for finding the ITN usage
find_usage <- function(use_rate){
  # Choose a distribution frequency in days 
  distribution_freq <- 3*365  
  
  ### Convert target usage to nets per capita ###
  
  result <- expand.grid(usage = target_usage, use_rate=use_rate, distribution_frequency=
                         distribution_freq)
  # Convert target usage to target access: access = usage/use rate
  # since use rate = the proportion of people with access to a net who slept under it
  # Access can only be calculated where usage <= use_rate (otherwise NA)
  result$access <- ifelse(result$usage/result$use_rate<=1, result$usage/result$use_rate, NA)
  # Convert access to nets per capita using the Loess curve 
  # For all usages where access does not exceed values on the curve
  result$percapita_nets <- NA
  result$percapita_nets[!(is.na(result$access)) & result$access<=max(loess_for_prediction$loess_predicted_access)] <- 
    approx(x=loess_for_prediction$loess_predicted_access,
           y=loess_for_prediction$percapita_nets, 
           xout=result$access[!(is.na(result$access)) & result$access<=max(loess_for_prediction$loess_predicted_access)])$y

  ### Convert equilibrium nets per capita to annual nets distributed per capita ###
  
  # Define smooth compact net loss function from Bertozzi-Villa, Amelia, et al. Nature communications 12.1 (2021): 1-12.
  # k is a fixed parameter from the paper
  net_loss_map <- function(t, k=20, hl = half_life) {
    
    # Convert half life into the time at which no nets are retained (nets=0) in days:
    l <- hl / sqrt(1 - k / (k - log(0.5)))
    
    prop_retained <- exp(k - k / (1 - (t / l)^2))
    prop_retained[t >= l] <- 0
    return(prop_retained)
  }
  
  # Plot net loss over time within a distribution cycle
  plot(x=seq(0,distribution_freq,1)/365, y = net_loss_map(t=seq(0,distribution_freq,1)),
       xlab="Years", ylab="Proportion of nets retained", ylim=c(0,1))
  
  # Use net loss function to calculate annual nets distributed per capita to maintain the 
  # input equilibrium nets per capita 
  result$annual_percapita_nets_distributed <- result$percapita_nets / 
    ((distribution_freq/365) * mean(net_loss_map(seq(0,distribution_freq,1))))
  
  # Note this represents the equilibrium nets to distribute annually for the given 
  # distribution cycle, but how this is applied can be chosen.
  # E.g. assume for a distribution cycle of 3 years, the annual nets distributed =
  # 0.2 * a population at risk of 10,000 = 2,000
  # So either distribute 2,000 nets in years 1, 2 and 3 each,
  # or distribute 6,000 nets in year 1 and 0 in years 2 and 3
  
  return(result)

}

# Choose a use rate: MEDIAN ----
use_rate <- median(cube_nat_level_annual$use_rate)
result <- find_usage(use_rate)

# Plot usage vs nets distributed
ggplot(result) +
  geom_line(aes(x=annual_percapita_nets_distributed, y = usage, 
                group=as.factor(round(use_rate,2)), colour=as.factor(round(use_rate,2))),
            size=1.1) +
  labs(x="Annual nets distributed per capita", y = "Net usage",
       colour= "Net use rate") +
  theme_classic() 

# save with one record per ITN use rate
result <- result %>% mutate(usage=round(usage,2))
saveRDS(result, paste0(datadir, "netz_median_use_rate.rds"))

# Reverse: save a dataset with usage values for a given value of nets distributed
conversion_usage_pcnets <- result[, c("usage", "annual_percapita_nets_distributed")]
conversion_usage_pcnets$annual_pcnets <- round(conversion_usage_pcnets$annual_percapita_nets_distributed, digits = 3)
conversion_usage_pcnets <- conversion_usage_pcnets[!is.na(conversion_usage_pcnets$annual_pcnets),]

pc_nets_annual <- seq(min(conversion_usage_pcnets$annual_pcnets), max(conversion_usage_pcnets$annual_pcnets), by= 0.001) #set a sequence of possible values for annual pcnets

usage <- c(rep(0,25))
for(i in 26: length(pc_nets_annual)){
  x <- pc_nets_annual[i]
  usage[i] <- conversion_usage_pcnets[max(which(abs(conversion_usage_pcnets$annual_pcnets - x) == 
                                                  min(abs(conversion_usage_pcnets$annual_pcnets - x)))), "usage"]
}

conversion_usage_pcnets <- as.data.frame(cbind(usage, pc_nets_annual))
saveRDS(conversion_usage_pcnets, paste0(datadir, "conversion_usage_pcnets_median_use_rate.rds"))

# Choose a use rate: MAXIMUM ----
use_rate <- max(cube_nat_level_annual$use_rate)
result <- find_usage(use_rate)
result <- result %>% mutate(usage=round(usage,2))
saveRDS(result, paste0(datadir, "netz_max_use_rate.rds"))

# Reverse: save a dataset with usage values for a given value of nets distributed
conversion_usage_pcnets <- result[, c("usage", "annual_percapita_nets_distributed")]
conversion_usage_pcnets$annual_pcnets <- round(conversion_usage_pcnets$annual_percapita_nets_distributed, digits = 3)
conversion_usage_pcnets <- conversion_usage_pcnets[!is.na(conversion_usage_pcnets$annual_pcnets),]

pc_nets_annual <- seq(min(conversion_usage_pcnets$annual_pcnets), max(conversion_usage_pcnets$annual_pcnets), by= 0.001) #set a sequence of possible values for annual pcnets

usage <- c(rep(0,25))
for(i in 26: length(pc_nets_annual)){
  x <- pc_nets_annual[i]
  usage[i] <- conversion_usage_pcnets[max(which(abs(conversion_usage_pcnets$annual_pcnets - x) == 
                                                  min(abs(conversion_usage_pcnets$annual_pcnets - x)))), "usage"]
}

conversion_usage_pcnets <- as.data.frame(cbind(usage, pc_nets_annual))
saveRDS(conversion_usage_pcnets, paste0(datadir, "conversion_usage_pcnets_max_use_rate.rds"))

# Assumptions -----------------------------------------------------------------

# Section 2: Extrapolation of higher access through continuation of Loess curve, 
# assumes access cannot exceed ~95%. 
# Section 2: Assume nets per capita for lower than observed access 
# (~9.5%) to be constant (i.e. low levels of usage incur the same fixed cost)
# Section 3: Net half life = median half life in Africa
# Section 3: Net use rate = median use rate in Africa
# Section 3: Nets are distributed every 3 years
# Section 3: Nets are lost according to a smooth compact function (MAP)


