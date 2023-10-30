## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE--------------------------------------------------------------
data.table::setDTthreads(2)

## ----setup--------------------------------------------------------------------
suppressPackageStartupMessages(library(tidyvpc, quietly = TRUE))
library(nlmixr2, quietly = TRUE)
library(magrittr)

## ----model-estimate-----------------------------------------------------------
one_compartment <- function() {
  ini({
    tka <- log(1.57); label("Ka")
    tcl <- log(2.72); label("Cl")
    tv <- log(31.5); label("V")
    eta_ka ~ 0.6
    eta_cl ~ 0.3
    eta_v ~ 0.1
    add_sd <- 0.7
  })
  # and a model block with the error specification and model specification
  model({
    ka <- exp(tka + eta_ka)
    cl <- exp(tcl + eta_cl)
    v <- exp(tv + eta_v)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add_sd)
  })
}

data_model <- theo_sd
data_model$WTSTRATA <- ifelse(data_model$WT < median(data_model$WT), "Low weight", "High weight")

fit <- nlmixr2(one_compartment, data_model, est="saem", saemControl(print=0))

## ----vpcsim-------------------------------------------------------------------
fit_vpcsim <- vpcSim(object = fit, keep = "WTSTRATA")

## ----vpc-std-setup------------------------------------------------------------
obs_data <- data_model[data_model$EVID == 0,]

vpc <-
  observed(obs_data, x=TIME, y=DV) %>%
  simulated(fit_vpcsim, x=time, y=sim) %>%
  stratify(~ WTSTRATA) %>%
  binning(bin = "jenks") %>%
  vpcstats()

## ----vpc-std-plot-------------------------------------------------------------
plot(vpc)

## ----vpc-predcorr-setup-------------------------------------------------------
# Add PRED to observed data
data_pred <- data_model[data_model$EVID == 0, ]
data_pred$PRED <- fit$PRED

vpc_predcorr <-
  observed(data_pred, x=TIME, y=DV) %>%
  simulated(fit_vpcsim, x=time, y=sim) %>%
  stratify(~ WTSTRATA) %>%
  binning(bin = "jenks") %>%
  predcorrect(pred=PRED) %>%
  vpcstats()

## ----vpc-predcorr-plot--------------------------------------------------------
plot(vpc_predcorr)

