% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataDoc.R
\docType{data}
\name{spap_avid_florbetapir}
\alias{spap_avid_florbetapir}
\title{sPAP Avid ADNI Florbetapir summaries}
\format{
A data frame with 604 rows and 19 variables
}
\usage{
data(spap_avid_florbetapir)
}
\description{
sPAP Avid ADNI Florbetapir summaries
}
\details{
A dataset with variables as follows:

\itemize{
  \item RID T  Participant roster ID
  \item VISCODE -4 -4 Visit code
  \item VISCODE2 -4 -4 Translated visit code
  \item EXAMDATE -4 -4 Examination Date
  \item DIAGRP2 T  clinical diagnosis close to the AV45 scan available during the data download in August 2012
  \item SPAP_GLOBAL_SUVR N Florbetapir (AV45) ratio (cortical grey matter/whole cerebellum), using sPAP software Average six region SUVr normalized by entire cerebellum for flroetbapir (AV-45) images calcualted using sPAP software. Refer to sPAP_Avid_Florbetapir_Analysis_Methods PDF on LONI for details
  \item SPAP_FRONTAL_SUVR N av45 uptake in cortical grey matter medial orbital frontal normalized to entire cerebellum in sPAP. Refer to sPAP_Avid_Florbetapir_Analysis_Methods PDF on LONI for details
  \item SPAP_TEMPORAL_SUVR N av45 uptake in cortical grey matter lateral temporal normalized to entire cerebellum in sPAP. Refer to sPAP_Avid_Florbetapir_Analysis_Methods PDF on LONI for details
  \item SPAP_ANTERIOR_CINGULATE_SUVR N av45 uptake in cortical grey matter anterior cingulate normalized to entire cerebellum in sPAP. Refer to sPAP_Avid_Florbetapir_Analysis_Methods PDF on LONI for details
  \item SPAP_POSTERIOR_CINGULATE_SUVR N av45 uptake in cortical grey matter posterior cingulate normalized to entire cerebellum is sPAP. Refer to sPAP_Avid_Florbetapir_Analysis_Methods PDF on LONI for details
  \item SPAP_PARIETAL_SUVR N av45 uptake in cortical grey matter parietal region normalized to entire cerebellum is sPAP. Refer to sPAP_Avid_Florbetapir_Analysis_Methods PDF on LONI for details
  \item SPAP_PRECUNEUS_SUVR N av45 uptake in cortical grey matter precuneus region normalized to entire cerebellum is sPAP. Refer to sPAP_Avid_Florbetapir_Analysis_Methods PDF on LONI for details
  \item AVID_STAGE_4_GLOBAL_SUVR N Florbetapir (AV45) ratio (cortical grey matter/whole cerebellum), using Avid semi-automated quantification method Average six region SUVr normalized by entire cerebellmu for florbetapir (AV-45) images calcualted using AVID semi-automated quantification. Refer to sPAP_Avid_Florbetapir_Analysis_Methods PDF on LONI for details
  \item AVID_STAGE_4_FRONTAL_MEDIAL_ORBITAL_SUVR N av45 uptake in cortical grey matter medial orbital frontal normalized to entire cerebellum using Avid semi-automated method. Refer to sPAP_Avid_Florbetapir_Analysis_Methods PDF on LONI for details
  \item AVID_STAGE_4_TEMPORAL_SUVR N av45 uptake in cortical grey matter lateral temporal normalized to entire cerebellum using Avid semi-automated method. Refer to sPAP_Avid_Florbetapir_Analysis_Methods PDF on LONI for details
  \item AVID_STAGE_4_PARIETAL_SUVR N av45 uptake in cortical grey matter parietal normalized to entire cerebellum using Avid semi-automated method. Refer to sPAP_Avid_Florbetapir_Analysis_Methods PDF on LONI for details
  \item AVID_STAGE_4_PRECUNEUS_SUVR N av45 uptake in cortical grey matter precuneus normalized to entire cerebellum using Avid semi-automated method. Refer to sPAP_Avid_Florbetapir_Analysis_Methods PDF on LONI for details
  \item AVID_STAGE_4_ANTERIOR_CINGULATE_SUVR N av45 uptake in cortical grey matter anterior cingulate normalized to entire cerebellum using Avid semi-automated method. Refer to sPAP_Avid_Florbetapir_Analysis_Methods PDF on LONI for details
  \item AVID_STAGE_4_POSTERIOR_CINGULATE_SUVR N av45 uptake in cortical grey matter s
}
}
\examples{
\dontrun{
describe(spap_avid_florbetapir)
}
}
\seealso{
\url{https://adni.bitbucket.io/reference/docs/SPAP_AVID_FLORBETAPIR/sPAP_Avid_Florbetapir_Analysis_Methods.pdf}
}
\keyword{datasets}
