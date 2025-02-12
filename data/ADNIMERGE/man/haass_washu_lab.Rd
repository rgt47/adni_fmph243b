% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataDoc.R
\docType{data}
\name{haass_washu_lab}
\alias{haass_washu_lab}
\title{Haass and Wash U lab summaries}
\format{
A data frame with 1878 rows and 15 variables
}
\usage{
data(haass_washu_lab)
}
\description{
Haass and Wash U lab summaries
}
\details{
A dataset with variables as follows:

\itemize{
  \item RID -4 -4 Participant roster ID
  \item VISCODE -4  Visit code
  \item VISCODE2 -4  Translated visit code
  \item USERDATE -4  Date record created
  \item EXAMDATE -4  Examination Date
  \item WU_STREM2 -4 pg/mL WashU CSF sTREM2 (pg/mL)
  \item WU_STREM2_CV -4 % WashU CSF sTREM2 Intraplate Coefficient of Variation (%)
  \item WU_STREM2CORRECTED -4 pg/mL WashU CSF sTREM2 (pg/mL) corrected values
  \item MSD_STREM2 -4 pg/mL Haass CSF sTREM2 (pg/mL)
  \item MSD_STREM2_CV -4 % Haass CSF sTREM2 Intraplate Coefficient of Variation (%)
  \item MSD_STREM2CORRECTED -4 pg/mL Haass CSF sTREM2 (pg/mL) corrected values
  \item TREM2OUTLIER -4  Correlation between the WashU and Haass CSF sTREM2 measurements
  \item MSD_PGRN -4 pg/mL Haass CSF PGRN (pg/mL)
  \item MSD_PGRN_CV -4 % Haass CSF PGRN Intraplate Coefficient of Variation (%)
  \item MSD_PGRNCORRECTED -4 pg/mL Haass CSF PGRN (pg/mL) corrected values
}
}
\examples{
\dontrun{
describe(haass_washu_lab)
}
}
\seealso{
\url{https://adni.bitbucket.io/reference/docs/HAASS_WASHU_LAB/ADNI_CSF_sTREM2_PGRN_methods_doc.pdf}
}
\keyword{datasets}
