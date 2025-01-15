## ----knitropts, results = 'hide', echo = FALSE, message = FALSE---------------
knitr::opts_chunk$set(fig.path='figure/progression-', echo = FALSE, message = FALSE, warning = FALSE, fig.width=6, fig.height=8/1.6, out.width='600px', fig.align = 'center', tidy = TRUE, comment = NA, cache = FALSE, cache.path='cache/enrollment-', results = 'markup')

## ----loadpackages, results = 'hide', echo = TRUE------------------------------
options(digits = 3)
library(knitr)
library(ADNIMERGE)
library(Hmisc)
library(msSurv)
library(zoo)
library(ggplot2)
theme_set(theme_bw())
source("http://adni.bitbucket.io/myfunctions.R")

cbbPalette <- c('#0072B2', '#D55E00', '#CC79A7', '#E69F00', '#56B4E9', '#009E73', '#F0E442', '#000000')
scale_colour_discrete <- function(...) scale_colour_manual(..., values = cbbPalette) 
scale_fill_discrete <- function(...) scale_fill_manual(..., values = cbbPalette)

## ----data, echo= FALSE--------------------------------------------------------
mmse_bl <- mmse
mmse_bl$VISCODE <- with(mmse_bl, ifelse(VISCODE == "sc", "bl", VISCODE))
neurobat_bl <- subset(neurobat, VISCODE != "bl")
neurobat_bl$VISCODE <- with(neurobat_bl, ifelse(VISCODE == "sc", "bl", VISCODE))
cdr_bl <- cdr
cdr_bl$VISCODE <- with(cdr_bl, ifelse(VISCODE == "sc", "bl", VISCODE))
dd <- merge(dxsum[, c("RID", "COLPROT", "VISCODE", "DXCHANGE", 
  paste("DXMPTR", 1:6, sep = ""))], 
  mmse_bl[, c("RID", "COLPROT", "VISCODE", "MMSCORE")], 
  by = c("RID", "COLPROT", "VISCODE"), all = TRUE)
dd <- merge(dd, neurobat_bl[, c("RID", "COLPROT", "VISCODE", "LDELTOTAL")], 
  by = c("RID", "COLPROT", "VISCODE"), all = TRUE)
dd <- merge(dd, cdr_bl[, c("RID", "COLPROT", "VISCODE", "CDGLOBAL", "CDMEMORY")], 
  by = c("RID", "COLPROT", "VISCODE"), all = TRUE)
dd <- merge(dd, ecogpt[, c("RID", "COLPROT", "VISCODE", "CONCERN")], 
  by = c("RID", "COLPROT", "VISCODE"), all = TRUE)
dd <- merge(dd, registry[, c("RID", "COLPROT", "VISCODE", "EXAMDATE")], 
  by = c("RID", "COLPROT", "VISCODE"), all = TRUE)

dd <- merge(dd, subset(registry, VISCODE == "bl" & COLPROT == ORIGPROT, 
  c("RID", "COLPROT", "VISCODE", "EXAMDATE")), 
  by = "RID", all = TRUE, suffixes = c("", ".bl"))
dd <- merge(dd, subset(arm, COLPROT == ORIGPROT, c("RID", "DX")), 
  by = "RID", all = TRUE)
dd <- merge(dd, subset(ptdemog, VISCODE == "sc" & COLPROT == ORIGPROT, c("RID", "PTEDUCAT")), 
 by = "RID", all = TRUE)

COHORT <- unique(subset(registry, 
  (RGCONDCT == 'Yes' | VISTYPE == "Standard") & VISCODE == 'bl')$RID)
dd <- subset(dd, RID %in% COHORT & !VISCODE %in% c("sc", "scmri") & !is.na(DXCHANGE))

dd$Year <- as.numeric(with(dd, EXAMDATE - EXAMDATE.bl)/365.25)

dd$Dx <- as.character(dd$DXCHANGE)
dd$Dx <- gsub("Stable: ", "", dd$Dx)
dd$Dx <- gsub("Conversion: ", "", dd$Dx)
dd$Dx <- gsub("Reversion: ", "", dd$Dx)
dd$Dx <- gsub("MCI to ", "", dd$Dx)
dd$Dx <- gsub("Dementia to ", "", dd$Dx)
dd$Dx <- gsub("NL to ", "", dd$Dx)

dd$Dx.sc <- gsub("AD", "Dementia", dd$DX)

dd$obj_memory_loss <- apply(dd[, c('PTEDUCAT', 'LDELTOTAL')], 1, obj_memory)
dd$subj_memory_comp <- with(dd, 
  DXMPTR1 == "Yes" | DXMPTR2 == "Yes" | CONCERN == "Yes" | Dx.sc == "EMCI")

dd <- subset(dd, !is.na(Year))
dd <- dd[order(dd$RID, dd$Year), ]
# LOCF for ojective memory complaint since Log Memory II not done at m06 and m18
for(rid in unique(dd$RID[!is.na(dd$obj_memory_loss)])){
	dd[dd$RID == rid, "obj_memory_loss"] <- 
    na.locf(dd[dd$RID == rid, "obj_memory_loss"], na.rm = TRUE)
}

# MMSCORE %in% 24:30 & subj_memory_comp & 
# CDGLOBAL == 0.5 & CDMEMORY >= 0.5

dd$Dx2 <- with(dd, 
  ifelse(Dx == "MCI" & 
  obj_memory_loss %in% c("EMCI", "NL"), "EMCI", Dx))
dd$Dx2 <- ifelse(is.na(dd$Dx2), dd$Dx, dd$Dx2)
dd$stage4 <- as.numeric(as.character(factor(dd$Dx2, 
  levels = c("NL", "EMCI", "MCI", "Dementia"),
  labels = 1:4)))
# Dx3, stage3: no EMCI allowed at followup
dd$Dx3 <- with(dd, ifelse(Year == 0, Dx2, Dx))
dd$stage3 <- as.numeric(as.character(factor(dd$Dx3, 
  levels = c("NL", "EMCI", "MCI", "Dementia"),
  labels = 1:4)))

# remove duplicated RID EXAMDATEs
dd <- subset(dd, !duplicated(paste(RID, EXAMDATE)))

# first to final transition (including transitions to EMCI)
dd_first_to_final_4 <- with(dd, make_first_to_final_trans(RID, stage4, Year))
dd_all_int_final_4 <- with(dd, make_all_intermediate_trans(RID, stage4, Year))

# first to final transition (excluding transitions to EMCI)
dd_first_to_final_3 <- with(dd, make_first_to_final_trans(RID, stage3, Year))
dd_all_int_final_3 <- with(dd, make_all_intermediate_trans(RID, stage3, Year))

# Nodes <- c("1", "2", "3", "4")
# Edges <- list("1" = list(edges = c("2", "3", "4")),
#               "2" = list(edges = c("1", "3", "4")),
#               "3" = list(edges = c("1", "2", "4")),
#               "4" = list(edges = c("3")))
# treeobj <- new("graphNEL", nodes = Nodes, edgeL = Edges,
#             edgemode = "directed")
# fit4 <- msSurv(dd_all_int_final_4, treeobj, bs = TRUE, LT = TRUE)

Nodes <- c("1", "2", "3", "4")
Edges <- list("1" = list(edges = c("3", "4")),
              "2" = list(edges = c("1", "3", "4")),
              "3" = list(edges = c("1", "4")),
              "4" = list(edges = c("3")))
treeobj <- new("graphNEL", nodes = Nodes, edgeL = Edges,
            edgemode = "directed")
fit3 <- msSurv(dd_all_int_final_3, treeobj, bs = TRUE, LT = TRUE)

## ----echo= FALSE--------------------------------------------------------------
dd3 <- dd_all_int_final_3
dd3$end.stage <- with(dd3, ifelse(end.stage == 0, start.stage, end.stage))
dd3$start.stage <- factor(dd3$start.stage, levels = c(1:4), labels = c("NL", "EMCI", "LMCI", "AD"))
dd3$end.stage <- factor(dd3$end.stage, levels = c(1:4), labels = c("NL", "EMCI", "LMCI", "AD"))
addmargins(with(subset(dd3, end.stage != 0), table(start.stage, end.stage)))

## -----------------------------------------------------------------------------
ddf <- dd_first_to_final_3
ddf$start.stage <- factor(ddf$start.stage, levels = c(1:4), labels = c("NL", "EMCI", "LMCI", "AD"))
ddf$end.stage <- factor(ddf$end.stage, levels = c(1:4), labels = c("NL", "EMCI", "LMCI", "AD"))
addmargins(with(ddf, table(start.stage, end.stage)))

## ----results = 'hide'---------------------------------------------------------
dx <- c("NL", "EMCI", "LMCI", "AD")
tab1 <- Pst(fit3, s = 0, t = 1.1)
tab2 <- Pst(fit3, s = 0, t = 2.1)
tab3 <- Pst(fit3, s = 0, t = 3.1)
tab4 <- Pst(fit3, s = 0, t = 4.1)

## -----------------------------------------------------------------------------
dimnames(tab1$Pst) <- list(dx, dx)
addmargins(tab1$Pst * 100)

## -----------------------------------------------------------------------------
dimnames(tab2$Pst) <- list(dx, dx)
addmargins(tab2$Pst * 100)

## -----------------------------------------------------------------------------
dimnames(tab3$Pst) <- list(dx, dx)
addmargins(tab3$Pst * 100)

## -----------------------------------------------------------------------------
dimnames(tab4$Pst) <- list(dx, dx)
addmargins(tab4$Pst * 100)

## -----------------------------------------------------------------------------
my_msSurv_plot(fit3, plot.type="transprob", trans=c("1 2", "1 3", "1 4")  , 
  state.names = c("NL", "EMCI", "LMCI", "AD"),
  timelab = "Year")
my_msSurv_plot(fit3, plot.type="transprob", trans=c("2 1", "2 3", "2 4")  , 
  state.names = c("NL", "EMCI", "LMCI", "AD"),
  timelab = "Year")
my_msSurv_plot(fit3, plot.type="transprob", trans=c("3 1", "3 2", "3 4")  , 
  state.names = c("NL", "EMCI", "LMCI", "AD"),
  timelab = "Year")

## -----------------------------------------------------------------------------
my_msSurv_plot(fit3, 
  state.names = c("NL", "EMCI", "LMCI", "AD"),
  timelab = "Year")

