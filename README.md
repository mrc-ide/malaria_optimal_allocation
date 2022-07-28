# Malaria Eradication Project :mosquito: :x:

## Summary

This modelling study explores optimal strategies to allocate limited resources to achieve global malaria eradication. Essentially, we attempt to address the following question: should resources initially focus on high burden countries, elimination countries, or a balance between the two?

The mathematical model used here adapts the [deterministic, compartmental model](https://github.com/mrc-ide/deterministic-malaria-model) of *P. falciparum* transmission previously fit to age-stratified data across a variety of transmission settings in sub-Saharan Africa ([Griffin et al. 2014](https://www.nature.com/articles/ncomms4136)). We use data from simulations of increasing insecticide-treated net usage in 25 different locations representing the world's range of *PfPR*<sub>2-10</sub> values and corresponding population sizes.

## Objectives

:heavy_check_mark: Create an optimization algorithm to minimize total global cases and total global population at risk of infection, with increasing budget values.

:heavy_check_mark: Compare optimization results with three standard strategies for resource allocation: 1) prioritizing high burden countries ("high burden to high impact"), 2) elimination countries ("shrinking the map"), and 3) a balance between the two ("proportional allocation").

:heavy_check_mark: Explore how competing aims in the short-term hamper achieving long-term eradication goals.
