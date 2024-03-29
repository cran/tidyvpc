## ---- warning = FALSE, echo = FALSE, message = FALSE--------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
options(datatable.print.nrows = 8)
library(tidyvpc)
library(ggplot2)
library(magrittr)
set.seed(1014)

## ---- echo=FALSE--------------------------------------------------------------
data.table::setDTthreads(2)

## ----message=FALSE------------------------------------------------------------
obs_data <- data.table::as.data.table(tidyvpc::obs_data)
head(obs_data)


## ----message=FALSE------------------------------------------------------------
sim_data <- data.table::as.data.table(tidyvpc::sim_data)
head(sim_data)

## ----echo=FALSE---------------------------------------------------------------
sim_data <- sim_data[REP <= 30]

## -----------------------------------------------------------------------------
obs_data <- obs_data[MDV == 0]
sim_data <- sim_data[MDV == 0]

## -----------------------------------------------------------------------------
obs_data$PRED <- sim_data[REP == 1, PRED]

## ----message = FALSE----------------------------------------------------------
library(tidyvpc)
vpc <- observed(obs_data, x = TIME, y = DV)

## -----------------------------------------------------------------------------
vpc <- observed(obs_data, x = TIME, y = DV) %>%
  simulated(sim_data, y = DV)

## ----message=FALSE, fig.width = 9, fig.height = 6, out.width=640--------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    binning(bin = NTIME)

## -----------------------------------------------------------------------------
vpc <- observed(obs_data, x = TIME, y = DV) %>%
  simulated(sim_data, y = DV) %>%
  binning(bin = "ntile", nbins = 9)

## -----------------------------------------------------------------------------
vpc <- observed(obs_data, x = TIME, y = DV) %>%
  simulated(sim_data, y = DV) %>%
  binning(bin = "breaks", breaks = c(1,5,7,9,10))

## ----warning=FALSE------------------------------------------------------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    binless()

## ----warning = FALSE----------------------------------------------------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    binless(qpred = c(0.1, 0.5, 0.9), optimize = FALSE, lambda = c(1,3,2))

## ----fig.width = 9, fig.height = 6, out.width=640-----------------------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    binning(bin = NTIME) %>%
    vpcstats()

plot(vpc)

## ----fig.width = 9, fig.height = 6, out.width=640-----------------------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    stratify(~ GENDER) %>%
    binning(bin = "pam", nbins = 7) %>%
    vpcstats()

plot(vpc)

## ----warning = FALSE, fig.width = 9, fig.height = 6, out.width=640------------

vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    stratify(~ GENDER + STUDY) %>%
    binless() %>%
    vpcstats()

plot(vpc)

## ----fig.width = 9, fig.height = 6, out.width=640, warning = FALSE------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    censoring(blq=(DV < 25), lloq=25) %>%
    binning(bin = "jenks", nbins = 5) %>%
    vpcstats()

plot(vpc, censoring.type = "blq")

## ----fig.width = 9, fig.height = 6, out.width=640, warning = FALSE------------
obs_data$LLOQ <- obs_data[, ifelse(STUDY == "Study A", 50, 25)]

vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    censoring(blq=(DV < LLOQ), lloq=LLOQ) %>%
    stratify(~ STUDY) %>%
    binning(bin = "pam", nbins = 4) %>%
    vpcstats(qpred = c(0.1, 0.5, 0.9))

plot(vpc, censoring.type = "blq", facet.scales = "fixed")

## ----fig.width = 9, fig.height = 6, out.width=640, warning = FALSE, eval = FALSE----
#  obs_data$ULOQ <- obs_data[, ifelse(STUDY == "Study A", 125, 100)]
#  
#  vpc <- observed(obs_data, x = TIME, y = DV) |>
#    simulated(sim_data, y = DV) |>
#    censoring(alq = DV > ULOQ, uloq = ULOQ) |>
#    stratify(~ STUDY) |>
#    binning(bin = NTIME) |>
#    vpcstats(qpred = c(0.1, 0.5, 0.9))

## ----fig.width = 9, fig.height = 6, out.width=640, warning = FALSE------------
obs_data$LLOQ <- obs_data[, ifelse(STUDY == "Study A", 50, 25)]
obs_data$ULOQ <- obs_data[, ifelse(STUDY == "Study A", 125, 100)]

vpc <- observed(obs_data, x = TIME, y = DV) |>
  simulated(sim_data, y = DV) |>
  censoring(blq = DV < LLOQ, lloq = LLOQ,  alq = DV > ULOQ, uloq = ULOQ) |>
  stratify(~ STUDY) |>
  binning(bin = NTIME) |>
  vpcstats(qpred = c(0.1, 0.5, 0.9))

plot(vpc, censoring.type = "both", facet.scales = "fixed")

## ----fig.width = 9, fig.height = 6, out.width=640-----------------------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    stratify(~GENDER) %>%
    binning(bin = NTIME) %>%
    predcorrect(pred=PRED) %>%
    vpcstats()

plot(vpc)

## ----warning=FALSE, fig.width = 9, fig.height = 6, out.width=640--------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    stratify(~GENDER) %>%
    binless(qpred = c(0.1, 0.5, 0.9), optimize = TRUE) %>%
    predcorrect(pred=PRED) %>%
    vpcstats()

plot(vpc)

## ----fig.width = 9, fig.height = 6, out.width=640, warning = FALSE------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    predcorrect(pred=PRED) %>%
    binless(optimize = FALSE, lambda = c(.95,3,1.2), span = .6) %>%
    vpcstats()

plot(vpc)

## ----fig.width = 9, fig.height = 6, out.width=640-----------------------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    stratify(~ GENDER + STUDY) %>%
    binning(stratum = list(GENDER = "M", STUDY = "Study A"), bin = "jenks", nbins = 5, by.strata = T) %>%
    binning(stratum = list(GENDER = "F", STUDY = "Study A"), bin = "centers", centers = c(0.5,3,5,10,15), by.strata = T) %>%
    binning(stratum = list(GENDER = "M", STUDY = "Study B"), bin = "kmeans", by.strata = T) %>%
    binning(stratum = list(GENDER = "F", STUDY = "Study B"), bin = "pam", nbins = 5, by.strata = T) %>%
    predcorrect(pred=PRED) %>%
    vpcstats()

plot(vpc)

## ----warning = FALSE, fig.width = 9, fig.height = 6, out.width=640------------
user_lambda <- data.frame(GENDER_F = c(2,4,2), GENDER_M = c(1.9,3,2.25) )

vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    stratify(~ GENDER) %>%
    binless(optimize = FALSE, lambda = user_lambda, span = c(.6, .85)) %>%
    predcorrect(pred=PRED) %>%
    vpcstats(qpred = c(0.1, 0.5, 0.9))

plot(vpc)

## ----fig.width = 9, fig.height = 6, out.width=640, warning = FALSE------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    binning(bin = "jenks", nbins = 7)

plot(vpc)


## -----------------------------------------------------------------------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    binning(bin = "jenks", nbins = 4) %>%
    vpcstats()

bin_information <- bininfo(vpc)
head(bin_information)


## ----fig.width = 9, fig.height = 6, out.width=640, warning=FALSE--------------
library(ggplot2)
obs_data$LLOQ <- obs_data[, ifelse(STUDY == "Study A", 50, 25)]

vpc <- observed(obs_data, x = TIME, y = DV) %>%
  simulated(sim_data, y = DV) %>%
  censoring(blq = DV < LLOQ, lloq = LLOQ) %>%
  stratify(~STUDY) %>%
  binning(bin = NTIME) %>%
  vpcstats(qpred = c(0.1, 0.5, 0.9))

ggplot(vpc$stats, aes(x = xbin)) + 
  facet_grid(~STUDY, scales = "free", as.table = FALSE) + 
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = qname, col = qname, group = qname),alpha = 0.1, col = NA) + 
  geom_line(aes(y = md, col = qname, group = qname)) +
  geom_line(aes(y = y, linetype = qname), size = 1) + 
  geom_hline(data=unique(obs_data[, .(STUDY, LLOQ)]), aes(yintercept=LLOQ), linetype="dotted", size=1) +
  geom_text(data = unique(vpc$data[, .(LLOQ), by = "STUDY"]), 
            aes(x = 10, y = LLOQ, label = paste("LLOQ", LLOQ, sep = "="), ), vjust = 1, hjust = 1) +
  scale_colour_manual(name = "Simulated Percentiles\nMedian (lines) 95% CI (areas)",
                    breaks = c("q0.1", "q0.5", "q0.9"), 
                    values = c("red", "blue", "red"), 
                    labels = c("10%", "50%", "90%")) + 
  scale_fill_manual(name = "Simulated Percentiles\nMedian (lines) 95% CI (areas)", 
                    breaks = c("q0.1", "q0.5", "q0.9"), 
                    values = c("red", "blue", "red"), 
                    labels = c("10%", "50%", "90%")) + 
  scale_linetype_manual(name = "Observed Percentiles\nMedian (lines) 95% CI (areas)", 
                    breaks = c("q0.1", "q0.5", "q0.9"), 
                    values = c("dotted", "solid", "dashed"), 
                    labels = c("10%", "50%", "90%")) + 
  guides(fill = guide_legend(order = 2), colour = guide_legend(order = 2), linetype = guide_legend(order = 1)) + 
  theme(legend.position = "top", legend.key.width = grid::unit(1, "cm")) + 
  labs(x = "TIME", y = "Concentration") + 
  geom_point(data = vpc$obs, aes(x = x, y = y), size = 1, alpha = 0.1, show.legend = FALSE) + 
  geom_vline(data = bininfo(vpc)[, .(x = sort(unique(c(xleft, xright)))), by = names(vpc$strat)],aes(xintercept = x), size = rel(0.5), col = "gray80") + 
  theme(panel.grid = element_blank()) + 
  geom_rug(data = bininfo(vpc)[, .(x = sort(unique(c(xleft, xright)))), by = names(vpc$strat)],aes(x = x), sides = "t", size = 1)


## -----------------------------------------------------------------------------
vpc <- observed(obs_data, x=TIME, y=DV) %>%
    simulated(sim_data, y=DV) %>%
    binning(bin = "jenks", nbins = 4) %>%
    vpcstats()

#Get vpcstats df
stats <- vpc$stats
#Get bininfo df
bin_information <- bininfo(vpc)
#Left join bin_info to vpcstats on bin
bin_information <- stats[bin_information, on = "bin"]
#Generate ymin
bin_information <- bin_information[, ymin := min(y), by = "bin"]
#Generate ymax
bin_information <- bin_information[, ymax := max(y), by = "bin"]
head(bin_information)

## ----fig.width = 9, fig.height = 6, out.width=640-----------------------------
ggplot(bin_information, aes(x = xbin)) + 
  geom_line(aes(y = md, col = qname, group = qname)) +
  geom_line(aes(y = y, linetype = qname), size = 1) +
  geom_rect(aes(xmin= xleft,xmax= xright, ymin =  ymin, ymax =  ymax),alpha = .1, col = "black", fill = "green") +
  geom_point(data = vpc$obs, aes(x = x, y = y), size = 1, alpha = 0.1, show.legend = FALSE) + 
  scale_colour_manual(name = "Simulated Percentiles\nMedian (lines) 95% CI (areas)",
                    breaks = c("q0.05", "q0.5", "q0.95"), 
                    values = c("red", "blue", "red"), 
                    labels = c("5%", "50%", "95%")) + 
  scale_linetype_manual(name = "Observed Percentiles\nMedian (lines) 95% CI (areas)", 
                    breaks = c("q0.05", "q0.5", "q0.95"), 
                    values = c("dotted", "solid", "dashed"), 
                    labels = c("5%", "50%", "95%")) + 
  geom_vline(data = bin_information[, .(x = sort(unique(c(xleft, xright))))],aes(xintercept = x), size = rel(0.5), col = "gray80") +
  geom_rug(data = bin_information[, .(x = sort(unique(c(xleft, xright))))],aes(x = x), sides = "t", size = 1) +
  guides(fill = guide_legend(order = 2), colour = guide_legend(order = 2), linetype = guide_legend(order = 1)) + 
  theme(legend.position = "top", legend.key.width = grid::unit(1, "cm")) + 
  labs(x = "TIME", y = "Concentration")

      

## ----fig.width = 9, fig.height = 6, out.width=640-----------------------------
obs_data$LLOQ <- obs_data[, ifelse(STUDY == "Study A", 50, 25)]

vpc <- observed(obs_data, x = TIME, y = DV) %>%
  simulated(sim_data, y = DV) %>%
  censoring(blq = DV < LLOQ, lloq = LLOQ) %>%
  stratify(~STUDY) %>%
  binning(bin = NTIME) %>%
  vpcstats(qpred = c(0.1, 0.5, 0.9))


ggplot(vpc$stats[bininfo(vpc), on=c("STUDY", "bin")], aes(x = xbin)) + 
  facet_grid(~STUDY, scales = "free", as.table = FALSE) + 
  geom_rect(aes(xmin = xleft, xmax = xright, ymin = lo, ymax = hi, fill = qname, col = qname, group = qname),alpha = 0.1, col = NA) + 
  geom_segment(aes(x = xleft, xend = xright, y = md, yend = md, col = qname, group = qname)) +
  geom_segment(aes(x = xleft, xend = xright, y = y, yend = y, linetype = qname), size = 1) +
  geom_line(aes(y = md, col = qname, group = qname)) +
  geom_line(aes(y = y, linetype = qname), size = 1) + 
  geom_hline(data=unique(obs_data[, .(STUDY, LLOQ)]), aes(yintercept=LLOQ), linetype="dotted", size=1) +
  geom_text(data = unique(vpc$data[, .(LLOQ), by = "STUDY"]), 
            aes(x = 10, y = LLOQ, label = paste("LLOQ", LLOQ, sep = "="), ), vjust = 1, hjust = 1) +
  scale_colour_manual(name = "Simulated Percentiles\nMedian (lines) 95% CI (areas)",
                      breaks = c("q0.1", "q0.5", "q0.9"), 
                      values = c("red", "blue", "red"), 
                      labels = c("10%", "50%", "90%")) + 
  scale_fill_manual(name = "Simulated Percentiles\nMedian (lines) 95% CI (areas)", 
                    breaks = c("q0.1", "q0.5", "q0.9"), 
                    values = c("red", "blue", "red"), 
                    labels = c("10%", "50%", "90%")) + 
  scale_linetype_manual(name = "Observed Percentiles\nMedian (lines) 95% CI (areas)", 
                        breaks = c("q0.1", "q0.5", "q0.9"), 
                        values = c("dotted", "solid", "dashed"), 
                        labels = c("10%", "50%", "90%")) + 
  guides(fill = guide_legend(order = 2), colour = guide_legend(order = 2), linetype = guide_legend(order = 1)) + 
  theme(legend.position = "top", legend.key.width = grid::unit(1, "cm")) + 
  labs(x = "TIME", y = "Concentration") + 
  geom_point(data = vpc$obs, aes(x = x, y = y), size = 1, alpha = 0.1, show.legend = FALSE) + 
  geom_vline(data = bininfo(vpc)[, .(x = sort(unique(c(xleft, xright)))), by = names(vpc$strat)],aes(xintercept = x), size = rel(0.5), col = "gray80") + 
  theme(panel.grid = element_blank()) + 
  geom_rug(data = bininfo(vpc)[, .(x = sort(unique(c(xleft, xright)))), by = names(vpc$strat)],aes(x = x), sides = "t", size = 1)


## ----fig.width = 9, fig.height = 6, out.width=640, warning = FALSE------------
obs_data$LLOQ <- obs_data[, ifelse(STUDY == "Study A", 50, 25)]

vpc <- observed(obs_data, x = TIME, y = DV) %>%
  simulated(sim_data, y = DV) %>%
  censoring(blq = DV < LLOQ, lloq = LLOQ) %>%
  stratify(~STUDY) %>%
  binning(bin = NTIME) %>%
  vpcstats(qpred = c(0.1, 0.5, 0.9))

ggplot(vpc$pctblq) + 
  facet_grid(~STUDY) +
  geom_ribbon(aes(x = xbin, ymin= lo, ymax = hi), fill = "red", alpha = .2) + 
  geom_line(aes(x = xbin, y = y)) + 
  geom_line(aes(x = xbin, y = md), color = "red") + 
  labs(x= "TIME", y= "% BLQ")

