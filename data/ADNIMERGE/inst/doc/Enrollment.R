## ----knitropts, results = 'hide', echo = FALSE, message = FALSE---------------
knitr::opts_chunk$set(fig.path='figure/enrollment-', echo = FALSE, message = FALSE, warning = FALSE, fig.width=6, fig.height=8/1.6, fig.align = 'center', tidy = TRUE, comment = NA, cache = FALSE, cache.path='cache/enrollment-', results = 'markup')

## ----loadpackages, results = 'hide', echo = TRUE------------------------------
options(digits = 3)
library(knitr)
library(ADNIMERGE)
library(ggplot2)
library(Hmisc)
library(gridExtra)
library(RColorBrewer)
library(tidyverse)
library(kableExtra)
source("https://adni.bitbucket.io/myfunctions.R")
theme_set(theme_bw())
dxpal <- c("#2b8cbe", "#74a9cf", "#fdcc8a", "#fc8d59","#e34a33")
scale_colour_discrete <- function(...) scale_colour_manual(..., values = dxpal) 
scale_fill_discrete <- function(...) scale_fill_manual(..., values = dxpal)

## ----data---------------------------------------------------------------------
dd <- adnimerge
dd$AGEC <- cut(dd$AGE, breaks = seq(39, 110, by = 10))

dd <- mutate(dd,
             AMYLOID = case_when(
  PIB.bl > 1.43 | AV45.bl > 1.11 | FBB.bl > 1.08 | ABETA.bl < 900 ~ "Positive", 
!is.na(PIB.bl) | !is.na(AV45.bl) | !is.na(FBB.bl) | !is.na(ABETA.bl)  ~ "Negative",
        is.na(PIB.bl) & is.na(AV45.bl) & is.na(FBB.bl) & is.na(ABETA.bl) ~ NA_character_
    ))

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
label(dd$AMYLOID) <- "Amyloid Status"
# indicate the bl withdrawals:
BLTDS <- subset(treatdis, VISCODE == "bl" & WDRAWTYPE == "Full")$RID
dd$VISCODE <- with(dd, ifelse(VISCODE == "bl" & RID %in% BLTDS, "blwd", VISCODE))

## ----baselinesummary1, results = 'asis'---------------------------------------
# Subsetting on RID < 2000 excludes ADNIGO/2 subjects
tab <- arsenal::tableby(DX.bl ~ AGE + AGEC + PTGENDER + PTEDUCAT + PTMARRY +  PTETHCAT + PTRACCAT + AMYLOID + CDRSB + ADAS13 + MMSE, 
  data = dd,subset = VISCODE %in% c('bl', 'blwd')) %>%
  summary( pfootnote = TRUE)

## ----baselinesummarygo2, results = 'asis'-------------------------------------
# Subsetting on RID > 2000 excludes ADNI1 subjects
arsenal::tableby(DX.bl ~ AGE + AGEC + PTGENDER + PTEDUCAT + PTMARRY + 
  PTETHCAT + PTRACCAT + AMYLOID + CDRSB + ADAS13 + MMSE, 
  data = dd, subset = VISCODE == 'bl') %>%
  summary( pfootnote = TRUE)

## ----adasmmsebp, fig=TRUE-----------------------------------------------------
dd1 <- subset(dd, VISCODE == "bl" & !is.na(DX.bl), 
  c("RID", "DX.bl", "ADAS13", "CDRSB", "MMSE", 
  "EcogPtMem", "EcogSPMem", "EcogPtTotal", "EcogSPTotal"))
ggplot(dd1, aes(DX.bl, ADAS13)) + 
  geom_boxplot(aes(fill = DX.bl), outlier.shape = NA) + 
  geom_jitter(size = 1, position = position_jitter(h=0)) + 
  theme(legend.position = "none")
ggplot(dd1, aes(DX.bl, MMSE)) + 
  geom_boxplot(aes(fill = DX.bl), outlier.shape = NA) + 
  geom_jitter(size = 1, position = position_jitter(h=0)) + 
  theme(legend.position = "none")

## ----cdrbp, fig=TRUE----------------------------------------------------------
ggplot(dd1, aes(DX.bl, CDRSB)) + 
  geom_boxplot(aes(fill = DX.bl), outlier.shape = NA) + 
  geom_jitter(size = 1, position = position_jitter(h=0)) + 
  theme(legend.position = "none")

## ----ecogmembp, fig=TRUE------------------------------------------------------
ggplot(dd1, aes(DX.bl, EcogPtMem)) + 
  geom_boxplot(aes(fill = DX.bl), outlier.shape = NA) + 
  geom_jitter(size = 1, position = position_jitter(h=0)) + 
  theme(legend.position = "none") +
  ylab("ECog Self Memory")
ggplot(dd1, aes(DX.bl, EcogSPMem)) + 
  geom_boxplot(aes(fill = DX.bl), outlier.shape = NA) + 
  geom_jitter(size = 1, position = position_jitter(h=0)) + 
  theme(legend.position = "none") +
  ylab("ECog Partner Memory")

## ----ecogtotbp, fig=TRUE------------------------------------------------------
ggplot(dd1, aes(DX.bl, EcogPtTotal)) + 
  geom_boxplot(aes(fill = DX.bl), outlier.shape = NA) + 
  geom_jitter(size = 1, position = position_jitter(h=0)) + 
  theme(legend.position = "none") +
  ylab("ECog Self Total")
ggplot(dd1, aes(DX.bl, EcogSPTotal)) + 
  geom_boxplot(aes(fill = DX.bl), outlier.shape = NA) + 
  geom_jitter(size = 1, position = position_jitter(h=0)) + 
  theme(legend.position = "none") +
  ylab("ECog Partner Total")

## ----results = 'asis'---------------------------------------------------------
cols <- c('ADAS13', 'MMSE', 'CDRSB', 'EcogPtMem', 'EcogPtTotal', 'EcogSPTotal')
for(cc in cols){
  dd1$x <- dd1[, cc]
  temp <- dd1 %>%
    as_tibble() %>%
    filter(!is.na(get(cc)))
  
  N <- temp %>%
    group_by(DX.bl) %>%
    dplyr::summarise(N = sum(!is.na(x))) %>%
    .$N
  
  tab <- with(temp, Hmisc::summarize(x, 
    by = DX.bl, FUN = summary.default, stat.name = 'Min.'))
  tab <- cbind(tab, N)
  cat("\n\n")
  cat("### ", label(dd1[, cc]), "\n")
  cat("\n\n")
  rownames(tab) <- tab[, 1]
  colnames <- c("N", "Q")
  print(kable(tab[, -1], format = 'simple'))
}

