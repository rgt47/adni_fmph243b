% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataDoc.R
\docType{data}
\name{ucd_adni3_wmh}
\alias{ucd_adni3_wmh}
\title{Univerity of California Davis White Matter Hyperintensity Volumes (ADNI3)}
\format{
A data frame with 720 rows and 21 variables
}
\usage{
data(ucd_adni3_wmh)
}
\description{
Univerity of California Davis White Matter Hyperintensity Volumes (ADNI3)
}
\details{
A dataset with variables as follows:

\itemize{
  \item RID N  Participant roster ID
  \item VISCODE -4 -4 Visit code
  \item VISCODE2 -4 -4 Translated visit code
  \item EXAMDATE D  Examination Date
  \item TCV N cc Supratentorial Cranial Volume
  \item TCB N cc Supratentorial Brain Volume
  \item TCC N cc Supratentorial CSF Volume
  \item CEREBRUM_GRAY N cc Cerebral Gray Matter Volume
  \item CEREBRUM_WHITE N cc Cerebral White Matter Volume
  \item HIPPOCAMPUS_L N cc Left Hippocampal Volume
  \item HIPPOCAMPUS_R N cc Right Hippocampal Volume
  \item HIPPOCAMPUS_TOT N cc Total Hippocampal Volume
  \item CSF N cc Intracranial CSF volume
  \item GRAY N cc Intracranial gray matter volume
  \item WHITE N cc Intracranial white matter volume
  \item WHITMATHYP N cc WMH volume
  \item ICV N cc Intracranial volume including posterior fossa
  \item MANUFACTURER T  MRI Machine Manufacturer
  \item MANUFACTURERS_MODEL T  MRI Machine Model
  \item SEGPROC T  Images used for segmentation
  \item SEGPROCIMG T  Analysis method used for segmentation
}
}
\examples{
\dontrun{
describe(ucd_adni3_wmh)
}
}
\keyword{datasets}
