#' @title Generate ADNI modified versions of the Preclinical Alzheimer's Cognitive Composite (PACC)
#'
#' @description The Preclinical Alzheimer's Cognitive Composite (\href{https://doi.org/10.1001/jamaneurol.2014.803}{PACC}) 
#' is a baseline standardized composite of
#'
#' @details 
#' \itemize{
#'   \item Free and Cued Selective Reminding Test (FCSRT)
#'   \item Logical Memory IIa Delayed Recall (LM)
#'   \item Digit Symbol Substitituion Test (DSST)
#'   \item Mini-Mental State Examination (MMSE)
#' }
#'
#' This function generates two modified versions of PACC scores based on ADNI data. FCSRT is not used in ADNI, so we use the
#' Delayed Recall portion of the Alzheimer's Disease Assessment Scale (ADAS) as a proxy. Score \code{mPACCdigit} uses the
#' DSST when available (ADNI1). \code{mPACCtrailsB} uses (log transformed) Trails B as a proxy for DSST. Raw component scores 
#' standardized according to the mean and standard deviation of baseline scores of ADNI subjects with normal cognition 
#' to create \code{Z} scores for each component \code{(Z=(raw - mean(raw.bl))/sd(raw.bl))}. 
#' The Z scores are reoriented if necessary so that greater scores reflect better performance. 
#' The composite is the sum of these Z scores. 
#'
#' \strong{Missing components:} At least two components must be present to produce a score. If more than two components are missing, 
#' the PACC will be \code{NA}.
#'
#' @param dd data frame to score with variables ADASQ4, LDELTOTAL, DIGITSCOR, MMSE, TRABSCOR, DX.bl, and VISCODE
#' @param keepComponents should components and z scores be returned? Defaults to FALSE so that just the composite scores are returned.
#'
#' @return The \code{data.frame} with appended columns for \code{mPACCdigit} and \code{mPACCtrailsB}
#' @references 
#' \itemize{
#'   \item Donohue MC, et al. The Preclinical Alzheimer Cognitive Composite: Measuring Amyloid-Related Decline. \emph{JAMA Neurol}. 2014;71(8):961–970. doi:10.1001/jamaneurol.2014.803 \url{http://dx.doi.org/10.1001/jamaneurol.2014.803}
#'   \item Donohue MC, Sperling RA, Petersen R, Sun C, Weiner MW, Aisen PS, for the Alzheimer’s Disease Neuroimaging Initiative. Association Between Elevated Brain Amyloid and Subsequent Cognitive Decline Among Cognitively Normal Persons. \emph{JAMA}. 2017;317(22):2305–2316. \url{http://dx.doi.org/10.1001/jama.2017.6669}
#' }
#' @author Michael Donohue \email{mdonohue@@usc.edu}
#' @seealso \code{\link{adnimerger}}, \code{\link{adnimerge}}
#' @examples
#' library(nlme)
#' library(dplyr)
#' library(multcomp)
#' library(Hmisc)
#' 
#' csf2numeric <- function(x){
#'   as.numeric(gsub("<", "", gsub(">", "", x)))
#' }
#' 
#' dd <- subset(adnimerge, DX.bl %in% c('CN', 'SMC') & !is.na(mPACCtrailsB)) %>%
#'   mutate(
#'     ABETA = csf2numeric(ABETA)
#'   )
#'
#' # identify those with elevated PIB PET at ANY visit OR 
#' # elevated BASELINE AV45 PET OR elevated BASELINE FBB PET OR
#' # low BASELINE CSF Abeta
#' # AV45 ~ PIB regression from Landau et al 2012:
#' elevatedAmyloid <- unique(c(
#'   subset(dd, 0.67*PIB + 0.15 > 1.11)$RID,
#'   subset(dd, VISCODE == 'bl' & AV45 > 1.11)$RID,
#'   subset(dd, VISCODE == 'bl' & FBB > 1.08)$RID,
#'   subset(dd, VISCODE == 'bl' & ABETA < 900)$RID
#' ))
#' anyAmyloid <- unique(subset(dd, !is.na(AV45.bl) | !is.na(PIB) | 
#'  | !is.na(FBB.bl) | !is.na(ABETA.bl))$RID)
#' dd <- mutate(dd,
#'   ElevatedAmyloid =
#'     ifelse(RID %in% elevatedAmyloid, 1,
#'     ifelse(RID %in% anyAmyloid, 0, NA)),
#'   m = Month.bl, m2 = Month.bl^2,
#'   APOEe4 = APOE4>0)
#'
#' summary(ElevatedAmyloid ~ APOEe4 + AGE + PTEDUCAT, data = dd,
#'   subset=VISCODE=='bl', method='reverse', overall=TRUE)
#'
#' # Quadratic time model:
#' fit <- lme(mPACCtrailsB ~ mPACCtrailsB.bl + APOEe4 + AGE + PTEDUCAT + m + m2 + (m + m2):ElevatedAmyloid,
#'   random = ~m|RID, data = dd, na.action=na.omit)
#'
#' Months <- seq(12,96,12)
#' elevated.design <- model.matrix(
#'   mPACCtrailsB ~ mPACCtrailsB.bl + APOEe4 + AGE + PTEDUCAT + m + m2 + (m + m2):ElevatedAmyloid,
#'   data = data.frame(mPACCtrailsB = 0, mPACCtrailsB.bl=0, APOEe4=TRUE, AGE=75, PTEDUCAT=12, ElevatedAmyloid=1, m=Months, m2=Months^2)
#' )
#'
#' normal.design <- model.matrix(
#'   mPACCtrailsB ~ mPACCtrailsB.bl + APOEe4 + AGE + PTEDUCAT + m + m2 + (m + m2):ElevatedAmyloid,
#'   data = data.frame(mPACCtrailsB = 0, mPACCtrailsB.bl=0, APOEe4=TRUE, AGE=75, PTEDUCAT=12, ElevatedAmyloid=0, m=Months, m2=Months^2)
#' )
#'
#' contrast.data <- elevated.design - normal.design
#' summary(multcomp::glht(fit, linfct = contrast.data))
#' @export
pacc <- function(dd, keepComponents = FALSE)
{
  require(dplyr)
  require(tidyr)
  holdNames <- colnames(dd)
	dd <- mutate(dd, log.TRABSCOR = log(TRABSCOR+1))
	
	bl.summary <- filter(dd, VISCODE == 'bl') %>% select(DX.bl, ADASQ4, LDELTOTAL, DIGITSCOR, log.TRABSCOR, MMSE) %>%
	  gather(VARIABLE, SCORE, -DX.bl) %>%
		filter(!is.na(SCORE)) %>%
		group_by(DX.bl, VARIABLE) %>%
		dplyr::summarize(
			N=n(),
			mean = mean(SCORE),
			sd = sd(SCORE))

  zscore <- function(x, var){ 
    (x - filter(bl.summary, DX.bl=='CN' & VARIABLE==var)$mean)/
		     filter(bl.summary, DX.bl=='CN' & VARIABLE==var)$sd
	}
	
	dd <- mutate(dd, 
			ADASQ4.z       = -zscore(ADASQ4, 'ADASQ4'),
			LDELTOTAL.z    =  zscore(LDELTOTAL, 'LDELTOTAL'),
			DIGITSCOR.z    =  zscore(DIGITSCOR, 'DIGITSCOR'),
			log.TRABSCOR.z = -zscore(log.TRABSCOR, 'log.TRABSCOR'),
			MMSE.z         =  zscore(MMSE, 'MMSE'))

  # check that all measure are positively correlated:
  corTest <- cor(select(dd, ADASQ4.z, LDELTOTAL.z, DIGITSCOR.z, log.TRABSCOR.z, MMSE.z), use='pairwise.complete.obs')
  if(any(corTest<0)) stop('Some PACC z scores are negatively correlated!')
  # GGally::ggpairs(select(dd, ADASQ4.z, LDELTOTAL.z, DIGITSCOR.z, log.TRABSCOR.z, MMSE.z))

  compscore <- function(x, n.components=4, n.missing=2){
		ifelse(sum(is.na(x))>n.missing, NA, mean(x, na.rm=TRUE))*n.components
	}
	
	dd$mPACCdigit <- apply(dd[, c('ADASQ4.z', 'LDELTOTAL.z', 'DIGITSCOR.z', 'MMSE.z')], 1, compscore)
	dd$mPACCtrailsB <- apply(dd[, c('ADASQ4.z', 'LDELTOTAL.z', 'log.TRABSCOR.z', 'MMSE.z')], 1, compscore)
  if(!keepComponents) dd <- dd[, c(holdNames, 'mPACCdigit', 'mPACCtrailsB')]  
  as.data.frame(dd)
}

                             
#' @title Score Montreal Cognitive Assessment (MOCA) without Education Adjustment.
#'
#' @description This function scores the MOCA.
#'
#' @param qlist Item list
#' @param trueCoding String for correct responses (default is `Correct`)
#' @param delayedRecallTrueCoding  defaults to `Correct with No Cue``
#' @param qScore defaults to `TRUE`
#' @param dd data frame to score.
#'
#' @return a data.frame with MOCA scores.
#' @references \url{http://www.parkinsons.va.gov/consortium/MOCA.asp}
#' @author Chung-Kai Sun \email{chungkas@@usc.edu}
#' @export
scoreMOCA <- function(dd, qlist = list(visuo = c(i1 = 'TRAILS', i2 = 'CUBE', 
                                                  i3Contour = 'CLOCKCON', 
                                                  i3Numbers = 'CLOCKNO', 
                                                  iHands = 'CLOCKHAN'), 
                                        naming = c(i4Lion = 'LION', 
                                                   i4Rhino = 'RHINO', 
                                                   i4Camel = 'CAMEL'), 
                                        memory = paste('IMMT', rep(1:2, c(5, 5)),
                                                       paste('W', 1:5, sep=''), sep=''), 
                                        attention = c(i6ForDigitSpan = 'DIGFOR', 
                                                      i6BackDigitSpan = 'DIGBACK', 
                                                      i6Vigilance = 'LETTERS', 
                                                      i6Serial1 = 'SERIAL1',
                                                      i6Serial2 = 'SERIAL2',
                                                      i6Serial3 = 'SERIAL3', 
                                                      i6Serial4 = 'SERIAL4', 
                                                      i6Serial5 = 'SERIAL5'), 
                                        language = c(i7SentenceReptition = paste('REPEAT', 1:2, 
                                                                                 sep=''), 
                                                     i8VerbalFluency = 'FFLUENCY'), 
                                        abstraction = c(i9Abstraction1 = 'ABSTRAN', 
                                                        i9Abstraction2 = 'ABSMEAS'), 
                                        delayedRecall = paste('DELW', 1:5, sep=''), 
                                        orientation = c(i11Date = 'DATE', i11Month = 'MONTH', 
                                                        i11Year = 'YEAR', i11Day = 'DAY',
                                                        i11Place = 'PLACE', i11City = 'CITY')), 
                           trueCoding = 'Correct', 
                           delayedRecallTrueCoding = 'Correct with No Cue', qScore = TRUE) { 

   if(!is.data.frame(dd)) {
     stop('dd is not a data.frame')  
   } else if (any(!unlist(qlist) %in% colnames(dd))) {
     stop(paste('Check qlist. dd does not contain', setdiff(unlist(qlist), colnames(dd)), sep=''))    
   } else {
     # non Binary variables
     # c('LETTERS', paste('SERIAL', 1:5, sep=''), 'FFLUENCY', paste('DELW', 1:5, sep=''))
     .i6Vigilance <- qlist$attention['i6Vigilance'] 
     .i6Serial <- qlist$attention[paste('i6Serial', 1:5, sep='')]
     .i8VerbalFluency <- qlist$language['i8VerbalFluency']
     .delayedRecall <- qlist$delayedRecall

     .notBinaryVariables <- c(.i6Vigilance, .i6Serial, .i8VerbalFluency, .delayedRecall)
     .binaryVariables <- setdiff(unlist(qlist), .notBinaryVariables)

     dd[, c(.binaryVariables, .i6Serial)] <-  sapply(dd[, c(.binaryVariables, .i6Serial)], 
            function(x) {
            # for binary variables
            # score is one if coding is trueCoding; missing score is NA; other coding is zero 
               x <- ifelse(x %in% trueCoding, 1, 
                           ifelse(is.na(x), NA, 0))
               return(x)                              
            })

     # i6Vigilance: 1, 0~1 error 
     dd[, .i6Vigilance] <- ifelse(dd[, .i6Vigilance] %in% c(0, 1), 1, 
                                  ifelse(is.na(dd[, .i6Vigilance]), NA, 0))

     # i6Serial7:  0, no correct subtractions; 1, 1 correct substraction; 
     #             2, 2-3 correct substraction; 3, 4-5 correct substraction.
     .i6SerialScore <- rowSums(dd[, .i6Serial])
     .i6SerialDic <- list(list(ans=0, score=0), list(ans=1, score=1), list(ans=c(2, 3), score=2), 
                           list(ans=c(4, 5), score=3))
     for (i in 1:length(.i6SerialDic)) {
        .ans <- .i6SerialDic[[i]]$ans
        .score <- .i6SerialDic[[i]]$score
        .i6SerialScore <- ifelse(.i6SerialScore %in% .ans, .score, .i6SerialScore)
     }
     dd[, 'i6Serial']  <- .i6SerialScore

     # i8VerbalFluency:  1, generates 11 words or more; otherwise 0
     dd[, .i8VerbalFluency] <- ifelse(is.numeric(dd[, .i8VerbalFluency]) 
                                      & dd[, .i8VerbalFluency] >= 11, 1, 
                                      ifelse(is.na(dd[, .i8VerbalFluency]), NA, 0))

     # delayedRecall: 1: Correct with No Cue; 
     #                0: No points are allocated for words recalled with a cue.
     dd[, .delayedRecall] <- sapply(dd[, .delayedRecall], 
                                    function(x) {
                                       x <- ifelse(x %in% delayedRecallTrueCoding, 1, 
                                                   ifelse(is.na(x), NA, 0))
                                       return(x)                              
                                    })

     if (qScore) {
        dd[, 'visuo'] <- rowSums(dd[, qlist$visuo]) 
        dd[, 'naming'] <- rowSums(dd[, qlist$naming])
        # dd[, memory] <- rowSums(dd[, memory])
        dd[, 'attention'] <- rowSums(dd[, c(qlist$attention[1:3], 'i6Serial')])   
        dd[, 'language'] <- rowSums(dd[, qlist$language])
        dd[, 'abstraction'] <- rowSums(dd[, qlist$abstraction])
        dd[, 'delayedRecall'] <- rowSums(dd[, .delayedRecall])
        dd[, 'orientation'] <- rowSums(dd[, qlist$orientation])
        # dd$MOCA1 <- rowSums(dd[, c('visuo', 'naming', 'attention', 'language', 'abstraction',
        #                            'delayedRecall', 'orientation')])
     }

    dd$MOCA <- rowSums(dd[, c(setdiff(.binaryVariables, qlist$memory), 
                               .i6Vigilance, 'i6Serial', .i8VerbalFluency, .delayedRecall)]) 
   return(dd)
   }
}

#' @title Score Everyday Cognition E-Cog
#'
#' @description This function scores the ECOG
#' following: http://www.sciencedirect.com/science/article/pii/S1552526011000896
#'
#' @param x row of data to score
#'
#' @return scored value
#' @references \url{http://www.sciencedirect.com/science/article/pii/S1552526011000896}
#' @author Michael Donohue \email{mdonohue@@usc.edu}
#' @export
ecog.score <- function(x){ 

  # coded data:
  # 1=1- Better or no change; 
  # 2=2- Questionable/occasionally worse; 
  # 3=3- Consistently a little worse; 
  # 4=4- Consistently much worse; 
  # 9=9- I don't know
  
  # coerce to numeric
  x <- as.numeric(unlist(lapply(strsplit(as.character(x), '- ', fixed = TRUE),
                                FUN = function(x) x[1])))

  # convert 9s to NA
  x[x == 9] <- NA

  missing <- sum(is.na(x))
  if(missing/length(x) < 0.5){
    return(mean(x, na.rm = TRUE))
  }else{
    return(NA)
  }
}

#' @title Function to merge commonly used ADNI data tables
#' 
#' @description This function merges commonly used ADNI data tables which come packaged in ADNIMERGE.
#'
#' @param csf use Roche Elecsys ('Roche', default) or INNO-BIA AlzBio3 immunoassay ('INNO')
#'
#' @usage adnimerger()
#'
#' @return Returns an analysis-ready dataset.
#' @export
adnimerger <- function(csf=c('Roche', 'INNO')){
  require(dplyr)
  require(tidyr)
  
  csf <- match.arg(csf)
  #--------------------------------------DEFINE INTERNAL FUNCTIONS-------------------------------------------

  # Merge recursively
  # Recursively merge data frames
  #
  # @arguments list of data frames to merge
  # @seealso \code{\link{merge_all}}
  # @keyword internal
  merge_recurse <- function(dfs, ...) {
    getCommonCols <- function(x, y) {
       intersect(colnames(x), colnames(y)) 
    }

    if (length(dfs) == 2) {
      merge(dfs[[1]], dfs[[2]], by = getCommonCols(dfs[[1]], dfs[[2]]), all=TRUE, sort=FALSE, ...)
    } else {
      merge(dfs[[1]], Recall(dfs[-1]), all=TRUE, sort=FALSE, ...)
    }
  }

  ordered <- function(x) x[with(x, order(RID, EXAMDATE)), ]

  # Check if data.frame has duplications
  is.duplicated <- function(x, bycols=c('RID', 'EXAMDATE')) {
   .return <- any(duplicated(Reduce(function(x, y) paste(x, y, sep='_'), x[, bycols])))
   if (.return) {
      .dd <- x[duplicated(Reduce(function(x, y) paste(x, y, sep='_'), x[, bycols])), bycols]
      .x <- x[Reduce(function(x, y) paste(x, y, sep='_'), x[, bycols]) %in% 
              Reduce(function(x, y) paste(x, y, sep='_'), .dd[, bycols]), ] 
      .x <- .x[eval(parse(text=paste('order(', paste('.x$', bycols, collapse=', '), ')', sep=''))), ]
      return(.x)
   } else { 
      return(.return)
   }
  }

  # sort ADNIMERGE datatable and remove duplications
  removeDups.adnimerge <- function(dd) {
     .dd <- dd[order(dd$RID, match(dd$COLPROT, c('ADNI1', 'ADNIGO', 'ADNI2', 'ADNI3'))), ]
     .dd <- .dd[!duplicated(paste(.dd$RID, .dd$VISCODE, sep='')), ]
     return(.dd) 
  }

  #--------------------------------------END OF FUNCTIONS-------------------------------------------

  CDRs <- c('CDMEMORY', 'CDORIENT', 'CDJUDGE', 'CDCOMMUN', 'CDHOME', 'CDCARE')
  cdr$CDRSB <- apply(cdr[, CDRs], 1, sum)
  adas$ADAS11 <- adas[, 'TOTSCORE']
  adas$ADAS13 <- adas[, 'TOTAL13']
  adas$ADASQ4 <- adas[, 'Q4SCORE']
  mmse$MMSE <- mmse[, 'MMSCORE']
  faq$FAQ <- faq[, 'FAQTOTAL']
  # MOCA Score without Education adjustment
  moca <- scoreMOCA(moca)
  moca <- moca[!is.na(moca$MOCA), ]

  roster$SITE <- unlist(lapply(strsplit(
                               as.character(roster$PTID), "_S_", fixed = TRUE), 
                               function(x){ x[1] }))
  apoeres$APOE4 <- rowSums(apoeres[, c("APGEN1","APGEN2")] == 4)
  apoego2$APOE4 <- rowSums(apoego2[, c("APGEN1","APGEN2")] == 4)
  apoe3$APOE4 <- rowSums(apoe3[, c("APGEN1","APGEN2")] == 4)
  apoe <- rbind(subset(apoeres, VISCODE == 'sc', c("RID", "APOE4")), apoego2[, c("RID", "APOE4")])
  apoe <- rbind(apoe, apoe3[, c("RID", "APOE4")])
  apoe <- apoe[!duplicated(apoe$RID, fromLast = TRUE),]

  pibpetsuvr$PIB   <- rowMeans(pibpetsuvr[, c('FRC', 'ACG', 'PRC', 'PAR')])
  ucberkeleyav45$AV45 <- ucberkeleyav45$SUMMARYSUVR_WHOLECEREBNORM
  ucberkeleyfbb$FBB <- ucberkeleyfbb$SUMMARYSUVR_WHOLECEREBNORM
  
  ucberkeleyfdg$ROINAME <- as.character(ucberkeleyfdg$ROINAME)
  ucberkeleyfdg1 <- ucberkeleyfdg %>%
    select(RID:UID, ROINAME, MEAN) %>%
    pivot_wider(names_from = 'ROINAME', values_from = 'MEAN') %>%
    mutate(FDG = metaroi/`pons-vermis`)
	if(any(with(ucberkeleyfdg1, duplicated(paste(RID, VISCODE))))){ warning("ucberkeleyfdg: duplicated RID/VISCODE") }
	ucberkeleyfdg1 <- subset(ucberkeleyfdg1, !duplicated(paste(RID, VISCODE)))
  fdg1 <- ucberkeleyfdg1

  # Estevez-Gonzalez, A., Kulisevsky, J., Boltes, A., Otermin, P., & Garcia-Sanchez, C. (2003). 
  # Rey verbal learning test is a useful tool for differential diagnosis in the preclinical phase
  # of Alzheimer's disease: comparison with mild cognitive impairment and normal aging. 
  # International Journal of Geriatric Psychiatry. 18 (11), 1021.

  neurobat$RAVLT.immediate <- rowSums(neurobat[, paste('AVTOT', 1:5, sep='')], na.rm=FALSE)
  neurobat$RAVLT.learning <- neurobat$AVTOT5-neurobat$AVTOT1
  neurobat$RAVLT.forgetting <- neurobat$AVTOT5-neurobat$AVDEL30MIN
  neurobat$RAVLT.perc.forgetting <- 100*neurobat$RAVLT.forgetting / neurobat$AVTOT5
  neurobat$RAVLT.perc.forgetting <- ifelse(neurobat$AVTOT5 == 0, NA, neurobat$RAVLT.perc.forgetting)

	# E-Cog
	# Change VISSPAT5 to LANG5  Aug.15 2012
	mem <- c('MEMORY1', 'MEMORY2', 'MEMORY3', 'MEMORY4', 'MEMORY5', 'MEMORY6', 'MEMORY7', 'MEMORY8')
	lang <-  c('LANG1', 'LANG2', 'LANG3', 'LANG4', 'LANG5', 'LANG6', 'LANG7', 'LANG8', 'LANG9')
	visspat <- c('VISSPAT1', 'VISSPAT2', 'VISSPAT3', 'VISSPAT4', 'VISSPAT6', 'VISSPAT7', 'VISSPAT8')
	plan <- c('PLAN1', 'PLAN2', 'PLAN3', 'PLAN4', 'PLAN5')
	organ <- c('ORGAN1', 'ORGAN2', 'ORGAN3', 'ORGAN4', 'ORGAN5', 'ORGAN6')
	divatt <- c('DIVATT1', 'DIVATT2', 'DIVATT3', 'DIVATT4')
 
  # ecogpt   Everyday Cognition - Participant Self-Report
  ecogpt$EcogPtMem <- apply(ecogpt[, mem], 1, ecog.score)
  ecogpt$EcogPtLang <- apply(ecogpt[, lang], 1, ecog.score)
  ecogpt$EcogPtVisspat <- apply(ecogpt[, visspat], 1, ecog.score)
  ecogpt$EcogPtPlan <- apply(ecogpt[, plan], 1, ecog.score)
  ecogpt$EcogPtOrgan <- apply(ecogpt[, organ], 1, ecog.score)
  ecogpt$EcogPtDivatt <- apply(ecogpt[, divatt], 1, ecog.score)
  ecogpt$EcogPtTotal <- apply(ecogpt[, c(mem, lang, visspat, plan, organ, divatt)], 1, ecog.score)
  
  # ecogsp   Everyday Cognition - Study Partner Report
  ecogsp$EcogSPMem <- apply(ecogsp[, mem], 1, ecog.score)
  ecogsp$EcogSPLang <- apply(ecogsp[, lang], 1, ecog.score)
  ecogsp$EcogSPVisspat <- apply(ecogsp[, visspat], 1, ecog.score)
  ecogsp$EcogSPPlan <- apply(ecogsp[, plan], 1, ecog.score)
  ecogsp$EcogSPOrgan <- apply(ecogsp[, organ], 1, ecog.score)
  ecogsp$EcogSPDivatt <- apply(ecogsp[, divatt], 1, ecog.score)
  ecogsp$EcogSPTotal <- apply(ecogsp[, c(mem, lang, visspat, plan, organ, divatt)], 1, ecog.score)
 
  #simplify DXCHANGE labels 
  dxsum$DX <- dxsum$DIAGNOSIS

  COLIDs <- c("COLPROT", "ORIGPROT", "RID", "VISCODE")
  KEYIDs <- c('RID', 'VISCODE')
  # RGCONDCT: Was this visit conducted? Yes
  registry <- filter(registry, VISTYPE != 'Not done' | RGCONDCT == 'Yes' | (COLPROT == 'ADNI1' & VISCODE == 'sc'))
  registry <- registry[order(registry$RID, match(registry$COLPROT, c('ADNI1', 'ADNIGO', 'ADNI2', 'ADNI3'))), ]
  # remove duplicates 
  registry <- removeDups.adnimerge(registry)
  
  # ADNI3 does not have arm table.
  arm$DX.bl <- arm$DX
  adni3arm <- subset(dxsum, ORIGPROT %in% "ADNI3" & COLPROT %in% "ADNI3" & VISCODE %in% 'bl', c("RID", "DIAGNOSIS"))
  colnames(adni3arm)[colnames(adni3arm) %in% "DIAGNOSIS"] <- "DX.bl"
  adni3arm$DX.bl <- as.character(adni3arm$DX.bl)
  # adni3 MCI is including both LMCI and EMCI, keep it as MCI for now, parse out LMCI and EMCI below.
  adni3arm$DX.bl <- sub("CN", "NL", sub("Dementia", "AD", adni3arm$DX.bl))

  dd <- merge(registry[, c(COLIDs, "EXAMDATE", "RGCONDCT", "RGSTATUS")],
              rbind(subset(arm, COLPROT == ORIGPROT, select = c("RID", "DX.bl")), adni3arm), 
              all.x = TRUE)

  dd$Baselined <- ifelse(dd$RID %in% 
    subset(registry, VISCODE %in% c("bl", "init") & (RGCONDCT == "Yes" | RGSTATUS == "Yes" | VISTYPE != 'Not done'))$RID, 1, 0)
  
  dd <- dd[dd$Baselined %in% 1, ]
  dd <- merge(dd,
              subset(ptdemog, COLPROT == ORIGPROT & VISCODE == "sc", 
               select = c("RID", "AGE", "PTGENDER", "PTEDUCAT", "PTETHCAT", "PTRACCAT", "PTMARRY")),
              all.x = TRUE)
  dd <- merge(dd,
              subset(roster, COLPROT == ORIGPROT, 
              select = c("RID", "PTID", "SITE")),
              all.x = TRUE)
  dd <- merge(dd, apoe, by = "RID", all.x = TRUE)

  pibpetsuvr <- ordered(pibpetsuvr)
  ucberkeleyav45 <- ordered(ucberkeleyav45)
  ucberkeleyfbb <- ordered(ucberkeleyfbb)

  dd <- merge(dd, subset(fdg1, , c("RID", "VISCODE", "FDG")), by = c("RID", "VISCODE"), all.x = TRUE)
  dd <- merge(dd, subset(pibpetsuvr, , c("RID", "VISCODE", "PIB")), by = c("RID", "VISCODE"), all.x = TRUE)
  dd <- merge(dd, subset(ucberkeleyav45, , c("RID", "VISCODE", "AV45")), by = c("RID", "VISCODE"), all.x = TRUE)
  dd <- merge(dd, subset(ucberkeleyfbb, , c("RID", "VISCODE", "FBB")), by = c("RID", "VISCODE"), all.x = TRUE)
  
  if(csf=='INNO'){
	  upb <- upennbiomk_master
	  upb <- mutate(upb, BATCH = gsub('UPENNBIOMK', '', BATCH)) %>%
	    mutate(BATCH = ifelse(BATCH=='', 1, BATCH)) %>%
	    select(RID, VISCODE, BATCH, ABETA, TAU, PTAU, ABETA_RAW, TAU_RAW, PTAU_RAW) %>%
	    gather(ASSAY, VALUE, -RID, -VISCODE, -BATCH) %>%
	    mutate(ASSAY = gsub('_', '.', ASSAY)) %>%
	    unite(ASSAY, ASSAY, BATCH, sep='.') %>%
	    spread(ASSAY, VALUE)  
	  dd <- merge(dd, upb, by = c("RID", "VISCODE"), all.x = TRUE)
  }else if(csf=='Roche'){
	  upb <- upennbiomk9
	  upb <- select(upb, RID, VISCODE, ABETA, TAU, PTAU)
	  dd <- merge(dd, upb, by = c("RID", "VISCODE"), all.x = TRUE)  	
  }else{
	stop("Unrecognized option for csf argument")
  }

  process.fs <- function(ucsfvol) {
    ucsfvol$ICV <- ucsfvol$ST10CV
    ucsfvol$Hippocampus <- (ucsfvol$ST29SV + ucsfvol$ST88SV)
    ucsfvol$Entorhinal <- (ucsfvol$ST24CV + ucsfvol$ST83CV)
    ucsfvol$Fusiform <- (ucsfvol$ST26CV + ucsfvol$ST85CV)
    ucsfvol$MidTemp <- (ucsfvol$ST40CV + ucsfvol$ST99CV)
    ucsfvol$Ventricles <- (ucsfvol$ST30SV +
                           ucsfvol$ST37SV + 
                           ucsfvol$ST89SV + 
                           ucsfvol$ST96SV)
    ucsfvol$WholeBrain <- apply(ucsfvol[, c('ST128SV', 'ST17SV', 'ST18SV', 'ST61SV', 
                                            'ST16SV', 'ST53SV', 'ST42SV', 'ST29SV', 
                                            'ST12SV', 'ST11SV', 'ST65SV', 'ST76SV', 
                                            'ST77SV', 'ST120SV', 'ST75SV', 'ST112SV', 
                                            'ST101SV', 'ST88SV', 'ST71SV', 'ST70SV',
                                            'ST124SV', 'ST147SV', 'ST148SV', 'ST150SV', 
                                            'ST151SV')], 1, sum)

    if('LHIPQC' %in% colnames(ucsfvol)){
      # make bilateral hippocampal qc variable
      ucsfvol$HIPQC <- with(ucsfvol, ifelse(LHIPQC=="Fail" | RHIPQC=="Fail", "Fail", "Pass"))
    }else{
      # use TEMPQC
      ucsfvol$HIPQC <- ucsfvol$TEMPQC
    }
    ucsfvol$Hippocampus <- ifelse(ucsfvol$HIPQC=="Fail", NA, ucsfvol$Hippocampus)
    ucsfvol$Entorhinal <- ifelse(ucsfvol$TEMPQC=="Fail", NA, ucsfvol$Entorhinal)
    ucsfvol$Fusiform <- ifelse(ucsfvol$TEMPQC=="Fail", NA, ucsfvol$Fusiform)
    ucsfvol$MidTemp <- ifelse(ucsfvol$TEMPQC=="Fail", NA, ucsfvol$MidTemp)
    ucsfvol$Ventricles <- ifelse(ucsfvol$VENTQC=="Fail", NA, ucsfvol$Ventricles)
    ucsfvol$WholeBrain <- ifelse(ucsfvol$OVERALLQC=="Fail", NA, ucsfvol$WholeBrain)
    ucsfvol
  }
	
  ## ucsffsx ----
  MRIs <- c('Hippocampus', 'Entorhinal', 'Fusiform', 'MidTemp', 'Ventricles', 'WholeBrain', 'ICV')
  colnames(ucsffsx)[colnames(ucsffsx) %in% 'ST19SV'] <- 'ST147SV'
  colnames(ucsffsx)[colnames(ucsffsx) %in% 'ST20SV'] <- 'ST150SV'
  colnames(ucsffsx)[colnames(ucsffsx) %in% 'ST78SV'] <- 'ST148SV'
  colnames(ucsffsx)[colnames(ucsffsx) %in% 'ST79SV'] <- 'ST151SV'
  ucsfvol1 <- process.fs(ucsffsx)
  ucsfvol1$FSVERSION <- subset(datadic, TBLNAME == 'UCSFFSX')[1,'CRFNAME']
  colnames(ucsffsx)[colnames(ucsffsx) %in% 'ST147SV'] <- 'ST19SV' 
  colnames(ucsffsx)[colnames(ucsffsx) %in% 'ST150SV'] <- 'ST20SV'
  colnames(ucsffsx)[colnames(ucsffsx) %in% 'ST148SV'] <- 'ST78SV'
  colnames(ucsffsx)[colnames(ucsffsx) %in% 'ST151SV'] <- 'ST79SV'
  ucsfvol1 <- ucsfvol1[!apply(ucsfvol1[, c('FSVERSION', MRIs)], 1, function(x) all(is.na(x))), ]
  ucsfvol1 <- ucsfvol1[!ucsfvol1$VISCODE %in% 'f', ]
  ucsfvol1 <- ucsfvol1[order(ucsfvol1$RID, ucsfvol1$EXAMDATE, 
		match(ucsfvol1$OVERALLQC, c('Pass', 'Partial')),
		-apply(ucsfvol1[, c('TEMPQC', 'FRONTQC', 'PARQC', 'INSULAQC', 'OCCQC', 'BGQC', 'CWMQC', 'VENTQC')] == 'Pass', 
		  1, sum, na.rm = TRUE),
		-as.numeric(ucsfvol1$RUNDATE)), ]
  # Remove duplicates in ucsfvol1 with the same RID, VISCODE
  ucsfvol1 <- ucsfvol1[!duplicated(with(ucsfvol1, paste(RID, VISCODE, sep=''))), ]

  ## ucsffsx51 ----
  ucsfvol2 <- process.fs(ucsffsx51)
  ucsfvol2$FSVERSION <- subset(datadic, TBLNAME == 'UCSFFSX51')[1,'CRFNAME']
  # Removed MRIs all NAs in ucsfvol2 
  ucsfvol2 <- ucsfvol2[!apply(ucsfvol2[, c('FSVERSION', MRIs)], 1, function(x) all(is.na(x))), ]
  ucsfvol2 <- ucsfvol2[order(ucsfvol2$RID, ucsfvol2$EXAMDATE, 
		match(ucsfvol2$OVERALLQC, c('Pass', 'Partial')),
		-apply(ucsfvol2[, c('TEMPQC', 'FRONTQC', 'PARQC', 'INSULAQC', 'OCCQC', 'BGQC', 'CWMQC', 'VENTQC', 'LHIPQC', 'RHIPQC')] == 'Pass', 
		  1, sum, na.rm = TRUE),
		match(ucsfvol2$IMAGETYPE, c('Accelerated T1', 'Non-Accelerated T1')),
		-as.numeric(ucsfvol2$RUNDATE)), ]
  # Remove duplicates in ucsfvol2 with the same RID, VISCODE 
  ucsfvol2[ucsfvol2$VISCODE %in% 'scmri', 'VISCODE'] <- 'sc'
  ucsfvol2 <- ucsfvol2[!duplicated(with(ucsfvol2, paste(RID, VISCODE, sep=''))), ]

  ## ucsffsx6 ----
  ucsfvol3 <- process.fs(ucsffsx6)
  ucsfvol3$FSVERSION <- subset(datadic, TBLNAME == 'UCSFFSX6')[1,'CRFNAME']
  # Removed MRIs all NAs in ucsfvol3 
  ucsfvol3 <- ucsfvol3[!apply(ucsfvol3[, c('FSVERSION', MRIs)], 1, function(x) all(is.na(x))), ]
  ucsfvol3 <- ucsfvol3[order(ucsfvol3$RID, ucsfvol3$EXAMDATE, 
    match(ucsfvol3$OVERALLQC, c('Pass', 'Partial')),
    -apply(ucsfvol3[, c('TEMPQC', 'FRONTQC', 'PARQC', 'INSULAQC', 'OCCQC', 'BGQC', 
      'CWMQC', 'VENTQC', 'HIPPOQC')] == 'Pass', 
      1, sum, na.rm = TRUE),
    -as.numeric(as.Date(ucsfvol3$RUNDATE, format='%m/%d/%Y'))), ]
  # Remove duplicates in ucsfvol2 with the same RID, VISCODE 
  ucsfvol3 <- ucsfvol3[!duplicated(with(ucsfvol3, paste(RID, VISCODE, sep=''))), ]
  if(!'LONISID' %in% colnames(ucsfvol3)){
    ucsfvol3$LONISID <- NA
  } 
  COLs <- c('RID', 'VISCODE', 'EXAMDATE', 'FSVERSION', 'LONISID', 'LONIUID', 'IMAGEUID', MRIs)
	# prefer records from ucsfvol2 over ucsfvol1
  ucsfvol <- rbind(
    ucsfvol1[!with(ucsfvol1, paste(RID, VISCODE)) %in% with(ucsfvol2, paste(RID, VISCODE)), COLs],
    ucsfvol2[, COLs], ucsfvol3[, COLs])
  mrinv <- subset(inventory, MODALITY_LONI == 'MRI' & QCSTAT_LONI == 'Available')
  mrinv$FLDSTRENG <- mrinv$EVENT
  mrinv <- subset(mrinv, !duplicated(paste(RID, IMAGEUID_LONI, FLDSTRENG)))
  dups <- mrinv$STUDYID_LONI[which(duplicated(mrinv$STUDYID_LONI))]
  mrinv <- subset(mrinv, !STUDYID_LONI %in% dups)
  nrow0 <- nrow(ucsfvol)
  ucsfvol <- merge(ucsfvol, mrinv[, c('RID', 'STUDYID_LONI', 'FLDSTRENG')], 
    by.x = c('RID', 'LONISID'), by.y = c('RID', 'STUDYID_LONI'),
    all.x = TRUE)
  if(nrow0 != nrow(ucsfvol)) stop('Problem in adnimerger with merging FLDSTRENG from inventory to UCSF FREESURFER data')
  ucsfvol <- ucsfvol[, c('RID', 'VISCODE', 'EXAMDATE', 'FLDSTRENG', 'FSVERSION', 'LONISID', 'LONIUID', 'IMAGEUID', MRIs)]

  # Refine levels & labels
  # change dd$DX.bl to reflect baseline diagnosis instead of screening diagnosis
  dd$DX.bl <- as.character(dd$DX.bl)
  dd$DX.bl <- with(dd,
    ifelse(RID %in% subset(dxsum, VISCODE == "bl" & DXCHANGE == "Conversion: MCI to Dementia")$RID, 
    "AD", DX.bl))
  dd$DX.bl <- with(dd, 
    ifelse(RID %in% subset(dxsum, VISCODE == "bl" & DXCHANGE == "Conversion: NL to Dementia")$RID,
    "AD", DX.bl))
  dd$DX.bl <- with(dd, 
    ifelse(RID %in% subset(dxsum, VISCODE == "bl" & DXCHANGE == "Conversion: NL to MCI")$RID & DX.bl != "EMCI",
    "LMCI", DX.bl))
  dd$DX.bl <-with(dd, 
    ifelse(RID %in% subset(dxsum, VISCODE == "bl" & DXCHANGE == "Reversion: MCI to NL")$RID & DX.bl != "EMCI",
    "NL", DX.bl))
  dd$DX.bl <- with(dd,
    ifelse(RID %in% subset(dxsum, VISCODE == "bl" & DXCHANGE == "Reversion: Dementia to MCI")$RID, 
    "LMCI", DX.bl))
    
  dd$SITE <- as.factor(dd$SITE)

  dd$PTETHCAT <- factor(dd$PTETHCAT, 
                        levels = c('Unknown', 'Not Hispanic or Latino', 'Hispanic or Latino'),
                        labels = c('Unknown', 'Not Hisp/Latino', 'Hisp/Latino'))
  dd$PTRACCAT <- factor(dd$PTRACCAT,
                        levels = c('American Indian or Alaskan Native', 'Asian', 
                                   'Native Hawaiian or Other Pacific Islander', 'Black or African American', 
                                   'White', 'More than one race', 'Unknown'),
                        labels =   c('Am Indian/Alaskan', 'Asian', 
                                     'Hawaiian/Other PI', 'Black', 
                                     'White', 'More than one', 'Unknown'))

  CDRs <- 'CDRSB'; CDR.LBs <- 'CDR-SB' 
  ADASs <- c('ADAS11', 'ADAS13', 'ADASQ4'); ADAS.LBs <- c('ADAS 11', 'ADAS 13 (including Delayed Word Recall and Number Cancellation)', 'ADAS Delayed Word Recall')
  MMSEs <- 'MMSE'; MMSE.LBs <- 'MMSE' 

  NEUROBATs <- c('RAVLT.immediate', 'RAVLT.learning', 'RAVLT.forgetting', 'RAVLT.perc.forgetting', 
    'LDELTOTAL', 'DIGITSCOR', 'TRABSCOR')  
  NEUROBAT.LBs <- c('RAVLT Immediate (sum of 5 trials)',
    'RAVLT Learning (trial 5 - trial 1)',
    'RAVLT Forgetting (trial 5 - delayed)',
    'RAVLT Percent Forgetting',
    'Logical Memory - Delayed Recall',
    'Digit Symbol Substitution',
    'Trails B')
  FAQs <- "FAQ"; FAQ.LBs <- 'FAQ'
  PACCs <- c("mPACCdigit", "mPACCtrailsB")
  PACC.LBs <- c("ADNI modified Preclinical Alzheimer's Cognitive Composite (PACC) with Digit Symbol Substitution", 
                "ADNI modified Preclinical Alzheimer's Cognitive Composite (PACC) with Trails B")
  UCSFVOLs <- c('FLDSTRENG', 'FSVERSION', 'IMAGEUID',
    "Ventricles", "Hippocampus", "WholeBrain", "Entorhinal", "Fusiform", "MidTemp", "ICV")
  UCSF.LBs <- c('MRI Field Strength', 'FreeSurfer Software Version', 'LONI Image ID',
    'UCSF Ventricles', 'UCSF Hippocampus', 'UCSF WholeBrain', 'UCSF Entorhinal', 'UCSF Fusiform', 'UCSF Med Temp', 'UCSF ICV')
  MOCAs <- c('MOCA'); MOCA.LBs <- c('MOCA')
  ECOGPTs <- paste('EcogPt', c('Mem', 'Lang', 'Visspat', 'Plan', 'Organ', 'Divatt', 'Total'), sep = '')
  ECOGPT.LBs <- paste('Pt ECog', c('Mem', 'Lang', 'Vis/Spat', 'Plan', 'Organ', 'Div atten', 'Total'), sep = ' - ')
  ECOGSPs <- c('EcogPtMem', 'EcogPtVisspat', 'EcogPtPlan', 'EcogPtOrgan', 'EcogPtDivatt', 'EcogPtLang', 'EcogPtTotal')
  ECOGSPs <- paste('EcogSP', c('Mem', 'Lang', 'Visspat', 'Plan', 'Organ', 'Divatt', 'Total'), sep = '')
  ECOGSP.LBs <- paste('SP ECog', c('Mem', 'Lang', 'Vis/Spat', 'Plan', 'Organ', 'Div atten', 'Total'), sep = ' - ')
  UPENNBIOMKs <- setdiff(colnames(upb), c('RID', 'VISCODE'))
  UPENNBIOMK.LBs <- gsub('.RAW', ' raw', UPENNBIOMKs)
  UPENNBIOMK.LBs <- gsub('.', ' batch ', UPENNBIOMK.LBs, fixed = TRUE)
  UPENNBIOMK.LBs <- paste('CSF', gsub(' batch MEDIAN', ' median of batches', UPENNBIOMK.LBs))
 
  covariates <- c(CDRs, ADASs, MMSEs, NEUROBATs, FAQs, PACCs, UCSFVOLs, MOCAs, ECOGPTs, ECOGSPs, UPENNBIOMKs)
  vlabels <- c(CDR.LBs, ADAS.LBs, MMSE.LBs, NEUROBAT.LBs, FAQ.LBs, PACC.LBs, UCSF.LBs, MOCA.LBs, ECOGPT.LBs, ECOGSP.LBs, UPENNBIOMK.LBs) 
   
  cdr <- cdr[!cdr$VISCODE %in% 'bl', ]
  cdr$VISCODE <- ifelse(cdr$VISCODE %in% 'sc', 'bl', cdr$VISCODE) 
  mmse <- mmse[!mmse$VISCODE %in% 'bl', ]
  mmse$VISCODE <- ifelse(mmse$VISCODE %in% 'sc', 'bl', mmse$VISCODE) 
  ucsfvol <- ucsfvol[!ucsfvol$VISCODE %in% 'bl', ]
  ucsfvol$VISCODE <- ifelse(ucsfvol$VISCODE %in% 'sc', 'bl', ucsfvol$VISCODE) 
  
  adas <- adas[!adas$VISCODE %in% 'sc', ]
  # carry screening LDELTOTAL forward to baseline
	neurobat <- full_join(neurobat, filter(neurobat, VISCODE == 'sc') %>% 
    select(RID, LDELTOTAL), by='RID', suffix=c('', '.sc')) %>%
  	mutate(LDELTOTAL = ifelse(VISCODE == 'bl', LDELTOTAL.sc, LDELTOTAL)) %>%
    filter(!VISCODE %in% 'sc')
  faq <- faq[!faq$VISCODE %in% 'sc', ]
  ecogpt <- ecogpt[!ecogpt$VISCODE %in% 'sc', ]
  ecogsp <- ecogsp[!ecogsp$VISCODE %in% 'sc', ] 

  ldf <- list(cdr=cdr[, c(COLIDs, CDRs)],
              adas=adas[, c(COLIDs, ADASs)],
              mmse=mmse[, c(COLIDs, MMSEs)], 
              neurobat=neurobat[, c(COLIDs, NEUROBATs)], 
              faq=faq[, c(COLIDs, FAQs)], 
              moca=moca[, c(COLIDs, MOCAs)],
              ecogpt=ecogpt[, c(COLIDs, ECOGPTs)],
              ecogsp=ecogsp[, c(COLIDs, ECOGSPs)])
              # ucsfvol=ucsfvol[, c(COLIDs, UCSFVOLs)])
  ldf <- lapply(ldf, removeDups.adnimerge)             
  ldf <- lapply(ldf, function(x) {
                         x <- x[, -grep('PROT', colnames(x))]
                         return(x) 
                      })

  mdf <- merge_recurse(ldf) 
  dd <- merge(dd, mdf, by=KEYIDs, all.x=TRUE)
  dup <- is.duplicated(dd, bycols=KEYIDs)
  # make Education adjustment for MOCA; 
  # Add one point for an individual who has 12 years or fewer of formal education, for a possible maximum of 30 points
  dd$MOCA <- ifelse(dd$PTEDUCAT <= 12, dd$MOCA + 1, dd$MOCA)
  dd$MOCA <- ifelse(dd$MOCA > 30, 30, dd$MOCA)
  dd <- merge(dd, ucsfvol[c('RID', 'VISCODE', UCSFVOLs)], by=c('RID', 'VISCODE'), all.x=TRUE)
  dd <- dd[!dd$VISCODE %in% c('sc', 'f', 'scmri', 'uns1'), ] 
  # dd <- merge(dd, subset(dxsum, VISCODE %in% c('sc', 'bl'), c('RID', 'VISCODE', 'DX')), by=c('RID', 'VISCODE'), all.x=TRUE)
  dd <- merge(dd, subset(dxsum, !duplicated(paste(RID, VISCODE)), c('RID', 'VISCODE', 'DX')), by=c('RID', 'VISCODE'), all.x=TRUE)

  # Parse out DX.bl LCMI and MCI for ADNI3
  # EMCI Subjects: delayed recall of one paragraph from Wechsler Memory Scale Logical Memory II:
  #   \geq 16 years: 9-11;
  #   8-15 years: 5-9
  #   0-7 years: 3-6
  # LMCI Subjects: delayed recall of one paragraph from Wechsler Memory Scale Logical Memory II:
  #   \geq 16 years: \leq 8;
  #   8-15 years: \leq 4;
  #   0-7 years: \leq 2 
  adni3LMCIids <- subset(dd, (VISCODE == 'bl' & DX.bl=='MCI' & ORIGPROT=='ADNI3') & ( 
    (PTEDUCAT>=16 &                LDELTOTAL<=8) |
    (PTEDUCAT>=8  & PTEDUCAT<=15 & LDELTOTAL<=4) |
    (PTEDUCAT<=7  &                LDELTOTAL<=2)))$RID
  dd$DX.bl <- with(dd, ifelse(DX.bl=='MCI' & RID %in% adni3LMCIids, 'LMCI', DX.bl))
  dd$DX.bl <- with(dd, ifelse(DX.bl=='MCI' & !RID %in% adni3LMCIids, 'EMCI', DX.bl))
  
  # Parse out DX.bl SMC from CN for ADNI3
  # SMC: CCI12TOT >= 16
  highCCI <- subset(cci, CCI12TOT >= 16 & COLPROT=='ADNI3' & VISCODE == 'sc')$RID
  dd$DX.bl <- with(dd, ifelse(DX.bl=='NL' & ORIGPROT=='ADNI3' & RID %in% highCCI, 'SMC', DX.bl))

  dd$DX.bl <- factor(dd$DX.bl, levels = c('NL', 'SMC', 'EMCI', 'LMCI', 'AD'),
                  labels = c('CN', 'SMC', 'EMCI', 'LMCI', 'AD'))

  dd <- mutate(dd, DIGITSCOR    = as.numeric(DIGITSCOR))
  dd <- pacc(dd)

  label(dd$RID) <- 'Participant roster ID'
  label(dd$VISCODE) <- 'Visit code'
  label(dd$COLPROT) <- 'Study protocol of data collection'
  label(dd$ORIGPROT) <- 'Original study protocol'
  label(dd$PTID) <- 'Original study protocol'
  label(dd$SITE) <- 'Site'
  label(dd$FDG) <- 'FDG-PET metaROI'
  label(dd$PIB) <- 'Average PIB SUVR of frontal cortex, anterior cingulate, precuneus cortex, and parietal cortex'
  label(dd$AV45) <- 'AV45 ratio (cortical grey matter/whole cerebellum) Summary florbetapir cortical SUVR normalized by whole cerebellum. See Jagust lab PDF on LONI for details and cutoff info'
  label(dd$FBB) <- 'FBB ratio (cortical grey matter/whole cerebellum) Summary florbetaben cortical SUVR normalized by whole cerebellum. See Jagust lab PDF on LONI for details and cutoff info'
  
  label(dd$DX.bl) <- "Baseline Dx"
  label(dd$PTEDUCAT) <- "Education"
  label(dd$AGE) <- "Age"
  label(dd$PTGENDER) <- "Sex"
  label(dd$PTETHCAT) <- "Ethnicity"
  label(dd$PTRACCAT) <- "Race"
  label(dd$PTMARRY) <- "Marital"
  label(dd$APOE4) <- 'Number of APOEe4 alleles'

  for(i in 1:length(covariates)) label(dd[, covariates[i]]) <- vlabels[i]
  for(i in 1:length(MRIs)) units(dd[, MRIs[i]]) <- "mm3"
 
  dd.bl <- subset(dd, COLPROT==ORIGPROT & VISCODE %in% c('bl'), 
    c('RID', 'EXAMDATE', covariates, 'FDG', 'PIB', 'AV45', 'FBB')) 
  dd <- merge(dd, dd.bl, by='RID', all.x=TRUE, suffixes=c('', '.bl'))
  
  dd$Years.bl <- as.numeric(dd$EXAMDATE - dd$EXAMDATE.bl)/365.25
  dd$Month.bl <- as.numeric(dd$EXAMDATE - dd$EXAMDATE.bl)/30.5
  dd$Month <- cut(dd$Month.bl, breaks=c(-1, 0, 4.5, seq(9, 597, 6)), labels=c(0, 3, seq(6, 594, 6)))
  # Months <- as.numeric(levels(dd$Month))[dd$Month]
  # dd$CompVis <- factor(dd$Month, labels= c('bl', 'm03', 'm06', paste('m', sort(unique(Months[Months > 10])), sep='')))

  dd$M <- as.numeric(gsub("m", "", gsub("bl", 0, dd$VISCODE)))
  dd$m <- as.factor(dd$M)
	
  label(dd$Month) <- 'Months since baseline'
  label(dd$M) <- 'Month since baseline'
  
	if(any(with(dd, duplicated(paste(RID, VISCODE))))){ warning("adnimerge: duplicated RID/VISCODE") }
	dd <- subset(dd, !duplicated(paste(RID, VISCODE)))
	
  dd[with(dd, order(RID, EXAMDATE)), !colnames(dd) %in% c('Baselined', 'm', 'RGCONDCT', 'RGSTATUS')]
}

