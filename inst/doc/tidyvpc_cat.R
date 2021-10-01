## ---- warning = FALSE, echo = FALSE, message = FALSE--------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7.5, fig.height = 5, fig.align = "center")
options(datatable.print.nrows = 8)
library(tidyvpc)
library(ggplot2)
library(magrittr)
library(data.table)
set.seed(1014)

## ----tidyvpc_setup------------------------------------------------------------
library(tidyvpc)
library(data.table)

obs_cat_data <- tidyvpc::obs_cat_data
sim_cat_data <- tidyvpc::sim_cat_data

obs_cat_data <- obs_cat_data[order(PID_code, agemonths)]
sim_cat_data <- sim_cat_data[order(Replicate, PID_code, IVAR)]


## ----echo=FALSE---------------------------------------------------------------
sim_cat_data <- sim_cat_data[Replicate <= 30]

## ----binning_1----------------------------------------------------------------

vpc <- observed(obs_cat_data, x = agemonths, yobs = zlencat) %>%
  simulated(sim_cat_data, ysim = DV) %>%
  binning(bin = round(agemonths, 0)) %>%
  vpcstats(vpc.type = "categorical")

plot(vpc, facet = TRUE, legend.position = "bottom", facet.scales = "fixed")

## ----binning_2----------------------------------------------------------------

vpc <- observed(obs_cat_data, x = agemonths, yobs = zlencat) %>%
  simulated(sim_cat_data, ysim = DV) %>%
  stratify(~ Country_ID_code) %>%
  binning(bin = "pam", nbins = 6) %>%
  vpcstats(vpc.type = "categorical", conf.level = .9, quantile.type = 6)

plot(vpc, facet = TRUE, legend.position = "bottom", facet.scales = "fixed")

## -----------------------------------------------------------------------------

vpc <- observed(obs_cat_data, x = agemonths, yobs = zlencat) %>%
  simulated(sim_cat_data, ysim = DV) %>%
  binless(optimize = TRUE) %>%
  vpcstats(vpc.type = "categorical", quantile.type = 6)

plot(vpc, facet = TRUE, legend.position = "bottom", facet.scales = "fixed")



## -----------------------------------------------------------------------------

vpc <- observed(obs_cat_data, x = agemonths, yobs = zlencat) %>%
  simulated(sim_cat_data, ysim = DV) %>%
  stratify(~ Country_ID_code) %>%
  binless(optimize = TRUE, optimization.interval = c(0,300)) %>%
  vpcstats(vpc.type = "categorical")

plot(vpc, facet = TRUE, legend.position = "bottom", facet.scales = "fixed")



## -----------------------------------------------------------------------------
sort(unique(obs_cat_data$zlencat))

## -----------------------------------------------------------------------------

sp_user <- list(p0 = 300,
                p1 = 50,
                p2 = 100)


## -----------------------------------------------------------------------------
vpc <- observed(obs_cat_data, x = agemonths, yobs = zlencat) %>%
  simulated(sim_cat_data, ysim = DV) %>%
  binless(optimize = FALSE, sp = sp_user) %>%
  vpcstats(vpc.type = "categorical", quantile.type = 6)

plot(vpc, facet = TRUE, legend.position = "bottom", facet.scales = "fixed")


## -----------------------------------------------------------------------------
sort(unique(obs_cat_data$Country_ID_code))
sort(unique(obs_cat_data$zlencat))

## -----------------------------------------------------------------------------
user_sp <- list(
           Country1_prob0 = 100,
           Country1_prob1 = 3,
           Country1_prob2 = 4,
           Country2_prob0 = 90,
           Country2_prob1 = 3,
           Country2_prob2 = 4,
           Country3_prob0 = 55,
           Country3_prob1 = 3,
           Country3_prob2 = 200)

## -----------------------------------------------------------------------------
vpc <- observed(obs_cat_data, x = agemonths, yobs = zlencat) %>%
  simulated(sim_cat_data, ysim = DV) %>%
  stratify(~ Country_ID_code) %>%
  binless(optimize = FALSE, sp = user_sp) %>%
  vpcstats(vpc.type = "categorical"
           , conf.level = 0.9
           , quantile.type = 6
  )
plot(vpc, facet = TRUE)


## -----------------------------------------------------------------------------
library(dplyr)

obs_cat_data <- obs_cat_data %>%
  mutate(gender = ifelse(PID_code %% 2 == 1, "male", "female"))

## -----------------------------------------------------------------------------
vpc <- observed(obs_cat_data, x = agemonths, yobs = zlencat) %>%
  simulated(sim_cat_data, ysim = DV) %>%
  stratify(~ gender + Country_ID_code)


## -----------------------------------------------------------------------------
sort(unique(obs_cat_data$gender))
sort(unique(obs_cat_data$Country_ID_code))
sort(unique(obs_cat_data$zlencat))


## -----------------------------------------------------------------------------
user_sp <- list(
  female.1.prob0 = 1,
  female.1.prob1 = 3,
  female.1.prob2 = 9,
  female.2.prob0 = 5,
  female.2.prob1 = 10,
  female.2.prob2 = 12,
  female.3.prob0 = 33,
  female.3.prob1 = 44,
  female.3.prob2 = 88,
  male.1.prob0 = 4,
  male.1.prob1 = 12,
  male.1.prob2 = 15,
  male.2.prob0 = 800,
  male.2.prob1 = 19,
  male.2.prob2 = 28,
  male.3.prob0 = 22,
  male.3.prob1 = 88,
  male.3.prob2 = 11
)

## -----------------------------------------------------------------------------
vpc <- observed(obs_cat_data, x = agemonths, yobs = zlencat) %>%
  simulated(sim_cat_data, ysim = DV) %>%
  stratify(~ gender + Country_ID_code) %>%
  binless(optimize = FALSE, sp = user_sp) %>%
  vpcstats(vpc.type = "categorical", conf.level = 0.9, quantile.type = 6)

plot(vpc, facet=TRUE)


