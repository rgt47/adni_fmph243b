% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataDoc.R
\docType{data}
\name{baipetnmrcfdg}
\alias{baipetnmrcfdg}
\title{Banner Alzheimer's Institute PET NMRC FDG Summaries}
\format{
A data frame with 3684 rows and 12 variables
}
\usage{
data(baipetnmrcfdg)
}
\description{
Banner Alzheimer's Institute PET NMRC FDG Summaries
}
\details{
A dataset with variables as follows:

\itemize{
  \item ORIGPROT T -4 Original study protocol
  \item RID N -4 Participant roster ID
  \item VISCODE T -4 Visit code
  \item VISCODE2 T -4 Translated visit code
  \item EXAMDATE D -4 Examination Date
  \item SCANDATE -4 -4 -4
  \item LONIUID N -4 LONI Unique Image Identifier
  \item RUNDATE D -4 DATE analysis performed
  \item STATUS T -4 Complete/Incomplete
  \item MODALITY T -4 The type of pet scan used for analysis
  \item HCI N unitless Hypometabolic Convergence Index
  \item SROI.AD -4 -4 -4
  \item SROI.MCI -4 -4 -4
}
}
\examples{
\dontrun{
describe(baipetnmrcfdg)
}
}
\seealso{
\url{https://adni.bitbucket.io/reference/docs/BAIPETNMRCFDG/BAI_HCI.pdf}, \url{https://adni.bitbucket.io/reference/docs/BAIPETNMRCFDG/BAI-PET-NMRC-FDG-Processing-Methods.pdf}, \url{https://adni.bitbucket.io/reference/docs/BAIPETNMRCFDG/Banner_HCIsROI\%20trajectories\%20prior\%20to\%20disease\%20onset\%202016\%2008.pdf}
}
\keyword{datasets}
