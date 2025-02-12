% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataDoc.R
\docType{data}
\name{reccmeds}
\alias{reccmeds}
\title{Concurrent Medications Log}
\format{
A data frame with 64362 rows and 29 variables
}
\usage{
data(reccmeds)
}
\description{
Concurrent Medications Log
}
\details{
A dataset with variables as follows:

\itemize{
  \item CMEVNUM T  If Adverse Event, enter event number:
  \item CMMEDID N NA medlist ID
  \item CMUNITSID T NA 
  \item CMFREQID T NA 
  \item CMROUTEID T NA 
  \item RID N  Participant roster ID
  \item SITEID N  Site ID
  \item VISCODE T  Visit code
  \item USERDATE S  Date record created
  \item RECNO N  
  \item EXAMDATE D  Examination Date
  \item CMHFS N  HFS Code
  \item CMBGN D  Date Began
  \item CMCONT N  Is the medication continuing?
  \item CMEND D  Date Ended
  \item USERDATE2 S  Date record last updated
  \item CMMED TRUE NA Free text meds entered into Concurrent Medications Log
  \item DONEPEZIL TRUE NA DONEPEZIL
  \item GALANTAMINE TRUE NA GALANTAMINE
  \item RIVASTIGMINE TRUE NA RIVASTIGMINE
  \item MEMANTINE TRUE NA MEMANTINE
  \item TACRINE TRUE NA TACRINE
  \item THYROXINE TRUE NA THYROXINE
  \item METHIMAZOLE TRUE NA METHIMAZOLE
  \item NSAID TRUE NA Nonsteroidal anti-inflammatory drugs
  \item BP TRUE NA Blood pressure drugs
  \item PPAR TRUE NA Peroxisome proliferator-activated receptors
}
}
\keyword{datasets}
