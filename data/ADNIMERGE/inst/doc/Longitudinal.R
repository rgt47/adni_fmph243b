## ----knitropts, results = 'hide', echo = FALSE, message = FALSE---------------
knitr::opts_chunk$set(fig.path='figure/enrollment-', echo = FALSE, message = FALSE, warning = FALSE, fig.width=10, fig.height=8/1.6,  fig.align = 'center', tidy = TRUE, comment = NA, cache = FALSE, cache.path='cache/longitudinal', results = 'markup')

## ----loadpackages, results = 'hide', echo = TRUE------------------------------
options(digits = 3)
library(knitr)
library(ADNIMERGE)
library(ggplot2)
library(Hmisc)
library(gridExtra)
library(RColorBrewer)
library(patchwork)

library(tidyverse)
library(kableExtra)
source("https://adni.bitbucket.io/myfunctions.R")
theme_set(theme_bw())
dx3pal <- c("#2b8cbe", "#fc8d59", "#e34a33")
scale_colour_discrete <- function(...) scale_colour_manual(..., values = dx3pal) 
scale_fill_discrete <- function(...) scale_fill_manual(..., values = dx3pal)

## ----data---------------------------------------------------------------------
dd <- adnimerge
dd <- dd %>% mutate(across(.cols = c("ABETA","TAU", "PTAU"), ~as.numeric(gsub(pattern = ">", replacement = "", .))))

dd$AGEC <- cut(dd$AGE, breaks = seq(39, 110, by = 10))

label(dd$PTEDUCAT) <- "Education"
label(dd$AGE) <- "Age"
label(dd$AGEC) <- "Age"
label(dd$PTGENDER) <- "Sex"
label(dd$PTETHCAT) <- "Ethnicity"
label(dd$PTRACCAT) <- "Race"
label(dd$PTMARRY) <- "Marital"
label(dd$CDRSB) <- "CDR-SB"
label(dd$ADAS11) <- "ADAS 11"
label(dd$ADAS13) <- "ADAS 13"
label(dd$MMSE) <- "MMSE"

# indicate the bl withdrawals:
BLTDS <- subset(treatdis, VISCODE == "bl" & WDRAWTYPE == "Full")$RID
dd$VISCODE <- with(dd, ifelse(VISCODE == "bl" & RID %in% BLTDS, "blwd", VISCODE))

## ----longsums, results = 'asis', echo = FALSE---------------------------------
longsums(dd, "ADAS13" )
longsums(dd, "MMSE")
longsums(dd, "CDRSB")
longsums(dd, "EcogPtMem")
longsums(dd, "EcogSPMem")
longsums(dd, "EcogPtTotal")
longsums(dd, "EcogSPTotal")
cat("\n\n")

## -----------------------------------------------------------------------------
# How to handle values such as "> 1700"?

# longsums(dd, "PTAU")
# longsums(dd, "ABETA")
# longsums(dd, "TAU")

## ----results = 'asis', echo = FALSE-------------------------------------------
longsums(dd, "Hippocampus")
longsums(dd, "WholeBrain")
longsums(dd, "Entorhinal")

