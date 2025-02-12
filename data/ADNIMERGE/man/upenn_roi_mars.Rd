% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataDoc.R
\docType{data}
\name{upenn_roi_mars}
\alias{upenn_roi_mars}
\title{UPENN - Hierarchical Parcellation of MRI Using Multi-atlas Labeling Methods}
\format{
A data frame with 840 rows and 269 variables
}
\usage{
data(upenn_roi_mars)
}
\description{
UPENN - Hierarchical Parcellation of MRI Using Multi-atlas Labeling Methods
}
\details{
A dataset with variables as follows:

\itemize{
  \item RID -4 -4 Participant roster ID
  \item VISCODE -4 -4 Visit code
  \item EXAMDATE -4 -4 Examination Date
  \item VERSION -4 -4 -4
  \item IMAGE_UID -4 -4 -4
  \item RUNDATE -4 -4 -4
  \item STATUS -4 -4 -4
  \item CBICA_ID -4 -4 -4
  \item DATE -4 -4 -4
  \item R702 -4 -4 ICV
  \item R701 -4 -4 TOTALBRAIN
  \item R601 -4 -4 GM
  \item R604 -4 -4 WM
  \item R606 -4 -4 GM_L
  \item R607 -4 -4 WM_L
  \item R613 -4 -4 GM_R
  \item R614 -4 -4 WM_R
  \item R501 -4 -4 CORPUS_CALLOSUM
  \item R502 -4 -4 CEREBELLUM
  \item R503 -4 -4 DEEP_WM_GM
  \item R504 -4 -4 FRONTAL
  \item R505 -4 -4 LIMBIC
  \item R506 -4 -4 OCCIPITAL
  \item R507 -4 -4 PARIETAL
  \item R508 -4 -4 TEMPORAL
  \item R509 -4 -4 VENTRICLE
  \item R510 -4 -4 CEREBELLUM_L
  \item R511 -4 -4 DEEP_WM_GM_L
  \item R512 -4 -4 FRONTAL_L
  \item R513 -4 -4 LIMBIC_L
  \item R514 -4 -4 OCCIPITAL_L
  \item R515 -4 -4 PARIETAL_L
  \item R516 -4 -4 TEMPORAL_L
  \item R517 -4 -4 VENTRICLE_L
  \item R518 -4 -4 CEREBELLUM_R
  \item R519 -4 -4 DEEP_WM_GM_R
  \item R520 -4 -4 FRONTAL_R
  \item R521 -4 -4 LIMBIC_R
  \item R522 -4 -4 OCCIPITAL_R
  \item R523 -4 -4 PARIETAL_R
  \item R524 -4 -4 TEMPORAL_R
  \item R525 -4 -4 VENTRICLE_R
  \item R401 -4 -4 BASAL_GANGLIA
  \item R402 -4 -4 DEEP_GM
  \item R403 -4 -4 DEEP_WM
  \item R404 -4 -4 FRONTAL_GM
  \item R405 -4 -4 FRONTAL_WM
  \item R406 -4 -4 LIMBIC_GM
  \item R407 -4 -4 OCCIPITAL_GM
  \item R408 -4 -4 OCCIPITAL_WM
  \item R409 -4 -4 PARIETAL_GM
  \item R410 -4 -4 PARIETAL_WM
  \item R411 -4 -4 TEMPORAL_GM
  \item R412 -4 -4 TEMPORAL_WM
  \item R413 -4 -4 BASAL_GANGLIA_L
  \item R414 -4 -4 DEEP_GM_L
  \item R415 -4 -4 DEEP_WM_L
  \item R416 -4 -4 FRONTAL_GM_L
  \item R417 -4 -4 FRONTAL_WM_L
  \item R418 -4 -4 LIMBIC_GM_L
  \item R419 -4 -4 OCCIPITAL_GM_L
  \item R420 -4 -4 OCCIPITAL_WM_L
  \item R421 -4 -4 PARIETAL_GM_L
  \item R422 -4 -4 PARIETAL_WM_L
  \item R423 -4 -4 TEMPORAL_GM_L
  \item R424 -4 -4 TEMPORAL_WM_L
  \item R425 -4 -4 BASAL_GANGLIA_R
  \item R426 -4 -4 DEEP_GM_R
  \item R427 -4 -4 DEEP_WM_R
  \item R428 -4 -4 FRONTAL_GM_R
  \item R429 -4 -4 FRONTAL_WM_R
  \item R430 -4 -4 LIMBIC_GM_R
  \item R431 -4 -4 OCCIPITAL_GM_R
  \item R432 -4 -4 OCCIPITAL_WM_R
  \item R433 -4 -4 PARIETAL_GM_R
  \item R434 -4 -4 PARIETAL_WM_R
  \item R435 -4 -4 TEMPORAL_GM_R
  \item R436 -4 -4 TEMPORAL_WM_R
  \item R301 -4 -4 FRONTAL_INFERIOR_GM
  \item R302 -4 -4 FRONTAL_INSULAR_GM
  \item R303 -4 -4 FRONTAL_LATERAL_GM
  \item R304 -4 -4 FRONTAL_MEDIAL_GM
  \item R305 -4 -4 FRONTAL_OPERCULAR_GM
  \item R306 -4 -4 LIMBIC_CINGULATE_GM
  \item R307 -4 -4 LIMBIC_MEDIALTEMPORAL_GM
  \item R308 -4 -4 OCCIPITAL_INFERIOR_GM
  \item R309 -4 -4 OCCIPITAL_LATERAL_GM
  \item R310 -4 -4 OCCIPITAL_MEDIAL_GM
  \item R311 -4 -4 PARIETAL_LATERAL_GM
  \item R312 -4 -4 PARIETAL_MEDIAL_GM
  \item R313 -4 -4 TEMPORAL_INFERIOR_GM
  \item R314 -4 -4 TEMPORAL_LATERAL_GM
  \item R315 -4 -4 TEMPORAL_SUPRATEMPORAL_GM
  \item R316 -4 -4 FRONTAL_INFERIOR_GM_L
  \item R317 -4 -4 FRONTAL_INSULAR_GM_L
  \item R318 -4 -4 FRONTAL_LATERAL_GM_L
  \item R319 -4 -4 FRONTAL_MEDIAL_GM_L
  \item R320 -4 -4 FRONTAL_OPERCULAR_GM_L
  \item R321 -4 -4 LIMBIC_CINGULATE_GM_L
  \item R322 -4 -4 LIMBIC_MEDIALTEMPORAL_GM_L
  \item R323 -4 -4 OCCIPITAL_INFERIOR_GM_L
  \item R324 -4 -4 OCCIPITAL_LATERAL_GM_L
  \item R325 -4 -4 OCCIPITAL_MEDIAL_GM_L
  \item R326 -4 -4 PARIETAL_LATERAL_GM_L
  \item R327 -4 -4 PARIETAL_MEDIAL_GM_L
  \item R328 -4 -4 TEMPORAL_INFERIOR_GM_L
  \item R329 -4 -4 TEMPORAL_LATERAL_GM_L
  \item R330 -4 -4 TEMPORAL_SUPRATEMPORAL_GM_L
  \item R331 -4 -4 FRONTAL_INFERIOR_GM_R
  \item R332 -4 -4 FRONTAL_INSULAR_GM_R
  \item R333 -4 -4 FRONTAL_LATERAL_GM_R
  \item R334 -4 -4 FRONTAL_MEDIAL_GM_R
  \item R335 -4 -4 FRONTAL_OPERCULAR_GM_R
  \item R336 -4 -4 LIMBIC_CINGULATE_GM_R
  \item R337 -4 -4 LIMBIC_MEDIALTEMPORAL_GM_R
  \item R338 -4 -4 OCCIPITAL_INFERIOR_GM_R
  \item R339 -4 -4 OCCIPITAL_LATERAL_GM_R
  \item R340 -4 -4 OCCIPITAL_MEDIAL_GM_R
  \item R341 -4 -4 PARIETAL_LATERAL_GM_R
  \item R342 -4 -4 PARIETAL_MEDIAL_GM_R
  \item R343 -4 -4 TEMPORAL_INFERIOR_GM_R
  \item R344 -4 -4 TEMPORAL_LATERAL_GM_R
  \item R345 -4 -4 TEMPORAL_SUPRATEMPORAL_GM_R
  \item R4 -4 -4 3rd Ventricle
  \item R11 -4 -4 4th Ventricle
  \item R23 -4 -4 Right Accumbens Area
  \item R30 -4 -4 Left Accumbens Area
  \item R31 -4 -4 Right Amygdala
  \item R32 -4 -4 Left Amygdala
  \item R35 -4 -4 Brain Stem
  \item R36 -4 -4 Right Caudate
  \item R37 -4 -4 Left Caudate
  \item R38 -4 -4 Right Cerebellum Exterior
  \item R39 -4 -4 Left Cerebellum Exterior
  \item R40 -4 -4 Right Cerebellum White Matter
  \item R41 -4 -4 Left Cerebellum White Matter
  \item R47 -4 -4 Right Hippocampus
  \item R48 -4 -4 Left Hippocampus
  \item R49 -4 -4 Right Inf Lat Vent
  \item R50 -4 -4 Left Inf Lat Vent
  \item R51 -4 -4 Right Lateral Ventricle
  \item R52 -4 -4 Left Lateral Ventricle
  \item R55 -4 -4 Right Pallidum
  \item R56 -4 -4 Left Pallidum
  \item R57 -4 -4 Right Putamen
  \item R58 -4 -4 Left Putamen
  \item R59 -4 -4 Right Thalamus Proper
  \item R60 -4 -4 Left Thalamus Proper
  \item R61 -4 -4 Right Ventral DC
  \item R62 -4 -4 Left Ventral DC
  \item R71 -4 -4 Cerebellar Vermal Lobules I-V
  \item R72 -4 -4 Cerebellar Vermal Lobules VI-VII
  \item R73 -4 -4 Cerebellar Vermal Lobules VIII-X
  \item R75 -4 -4 Left Basal Forebrain
  \item R76 -4 -4 Right Basal Forebrain
  \item R81 -4 -4 frontal lobe WM right
  \item R82 -4 -4 frontal lobe WM left
  \item R83 -4 -4 occipital lobe WM right
  \item R84 -4 -4 occipital lobe WM left
  \item R85 -4 -4 parietal lobe WM right
  \item R86 -4 -4 parietal lobe WM left
  \item R87 -4 -4 temporal lobe WM right
  \item R88 -4 -4 temporal lobe WM left
  \item R89 -4 -4 fornix right
  \item R90 -4 -4 fornix left
  \item R91 -4 -4 anterior limb of internal capsule right
  \item R92 -4 -4 anterior limb of internal capsule left
  \item R93 -4 -4 posterior limb of internal capsule inc. cerebral peduncle right
  \item R94 -4 -4 posterior limb of internal capsule inc. cerebral peduncle left
  \item R95 -4 -4 corpus callosum
  \item R100 -4 -4 Right ACgG  anterior cingulate gyrus
  \item R101 -4 -4 Left ACgG  anterior cingulate gyrus
  \item R102 -4 -4 Right AIns  anterior insula
  \item R103 -4 -4 Left AIns  anterior insula
  \item R104 -4 -4 Right AOrG  anterior orbital gyrus
  \item R105 -4 -4 Left AOrG  anterior orbital gyrus
  \item R106 -4 -4 Right AnG   angular gyrus
  \item R107 -4 -4 Left AnG   angular gyrus
  \item R108 -4 -4 Right Calc  calcarine cortex
  \item R109 -4 -4 Left Calc  calcarine cortex
  \item R112 -4 -4 Right CO    central operculum
  \item R113 -4 -4 Left CO    central operculum
  \item R114 -4 -4 Right Cun   cuneus
  \item R115 -4 -4 Left Cun   cuneus
  \item R116 -4 -4 Right Ent   entorhinal area
  \item R117 -4 -4 Left Ent   entorhinal area
  \item R118 -4 -4 Right FO    frontal operculum
  \item R119 -4 -4 Left FO    frontal operculum
  \item R120 -4 -4 Right FRP   frontal pole
  \item R121 -4 -4 Left FRP   frontal pole
  \item R122 -4 -4 Right FuG   fusiform gyrus
  \item R123 -4 -4 Left FuG   fusiform gyrus
  \item R124 -4 -4 Right GRe   gyrus rectus
  \item R125 -4 -4 Left GRe   gyrus rectus
  \item R128 -4 -4 Right IOG   inferior occipital gyrus
  \item R129 -4 -4 Left IOG   inferior occipital gyrus
  \item R132 -4 -4 Right ITG   inferior temporal gyrus
  \item R133 -4 -4 Left ITG   inferior temporal gyrus
  \item R134 -4 -4 Right LiG   lingual gyrus
  \item R135 -4 -4 Left LiG   lingual gyrus
  \item R136 -4 -4 Right LOrG  lateral orbital gyrus
  \item R137 -4 -4 Left LOrG  lateral orbital gyrus
  \item R138 -4 -4 Right MCgG  middle cingulate gyrus
  \item R139 -4 -4 Left MCgG  middle cingulate gyrus
  \item R140 -4 -4 Right MFC   medial frontal cortex
  \item R141 -4 -4 Left MFC   medial frontal cortex
  \item R142 -4 -4 Right MFG   middle frontal gyrus
  \item R143 -4 -4 Left MFG   middle frontal gyrus
  \item R144 -4 -4 Right MOG   middle occipital gyrus
  \item R145 -4 -4 Left MOG   middle occipital gyrus
  \item R146 -4 -4 Right MOrG  medial orbital gyrus
  \item R147 -4 -4 Left MOrG  medial orbital gyrus
  \item R148 -4 -4 Right MPoG  postcentral gyrus medial segment
  \item R149 -4 -4 Left MPoG  postcentral gyrus medial segment
  \item R150 -4 -4 Right MPrG  precentral gyrus medial segment
  \item R151 -4 -4 Left MPrG  precentral gyrus medial segment
  \item R152 -4 -4 Right MSFG  superior frontal gyrus medial segment
  \item R153 -4 -4 Left MSFG  superior frontal gyrus medial segment
  \item R154 -4 -4 Right MTG   middle temporal gyrus
  \item R155 -4 -4 Left MTG   middle temporal gyrus
  \item R156 -4 -4 Right OCP   occipital pole
  \item R157 -4 -4 Left OCP   occipital pole
  \item R160 -4 -4 Right OFuG  occipital fusiform gyrus
  \item R161 -4 -4 Left OFuG  occipital fusiform gyrus
  \item R162 -4 -4 Right OpIFG opercular part of the inferior frontal gyrus
  \item R163 -4 -4 Left OpIFG opercular part of the inferior frontal gyrus
  \item R164 -4 -4 Right OrIFG orbital part of the inferior frontal gyrus
  \item R165 -4 -4 Left OrIFG orbital part of the inferior frontal gyrus
  \item R166 -4 -4 Right PCgG  posterior cingulate gyrus
  \item R167 -4 -4 Left PCgG  posterior cingulate gyrus
  \item R168 -4 -4 Right PCu   precuneus
  \item R169 -4 -4 Left PCu   precuneus
  \item R170 -4 -4 Right PHG   parahippocampal gyrus
  \item R171 -4 -4 Left PHG   parahippocampal gyrus
  \item R172 -4 -4 Right PIns  posterior insula
  \item R173 -4 -4 Left PIns  posterior insula
  \item R174 -4 -4 Right PO    parietal operculum
  \item R175 -4 -4 Left PO    parietal operculum
  \item R176 -4 -4 Right PoG   postcentral gyrus
  \item R177 -4 -4 Left PoG   postcentral gyrus
  \item R178 -4 -4 Right POrG  posterior orbital gyrus
  \item R179 -4 -4 Left POrG  posterior orbital gyrus
  \item R180 -4 -4 Right PP    planum polare
  \item R181 -4 -4 Left PP    planum polare
  \item R182 -4 -4 Right PrG   precentral gyrus
  \item R183 -4 -4 Left PrG   precentral gyrus
  \item R184 -4 -4 Right PT    planum temporale
  \item R185 -4 -4 Left PT    planum temporale
  \item R186 -4 -4 Right SCA   subcallosal area
  \item R187 -4 -4 Left SCA   subcallosal area
  \item R190 -4 -4 Right SFG   superior frontal gyrus
  \item R191 -4 -4 Left SFG   superior frontal gyrus
  \item R192 -4 -4 Right SMC   supplementary motor cortex
  \item R193 -4 -4 Left SMC   supplementary motor cortex
  \item R194 -4 -4 Right SMG   supramarginal gyrus
  \item R195 -4 -4 Left SMG   supramarginal gyrus
  \item R196 -4 -4 Right SOG   superior occipital gyrus
  \item R197 -4 -4 Left SOG   superior occipital gyrus
  \item R198 -4 -4 Right SPL   superior parietal lobule
  \item R199 -4 -4 Left SPL   superior parietal lobule
  \item R200 -4 -4 Right STG   superior temporal gyrus
  \item R201 -4 -4 Left STG   superior temporal gyrus
  \item R202 -4 -4 Right TMP   temporal pole
  \item R203 -4 -4 Left TMP   temporal pole
  \item R204 -4 -4 Right TrIFG triangular part of the inferior frontal gyrus
  \item R205 -4 -4 Left TrIFG triangular part of the inferior frontal gyrus
  \item R206 -4 -4 Right TTG   transverse temporal gyrus
  \item R207 -4 -4 Left TTG   transverse temporal gyrus
}
}
\examples{
\dontrun{
describe(upenn_roi_mars)
}
}
\seealso{
\url{https://adni.bitbucket.io/reference/docs/UPENN_ROI_MARS/Multi-atlas_ROI_ADNI_Methods_mod_April2016.pdf}
}
\keyword{datasets}
