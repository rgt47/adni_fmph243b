% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataDoc.R
\docType{data}
\name{ucberkeleyav1451}
\alias{ucberkeleyav1451}
\title{UC Berkeley - AV1451 analysis}
\format{
A data frame with 1450 rows and 242 variables
}
\usage{
data(ucberkeleyav1451)
}
\description{
UC Berkeley - AV1451 analysis
}
\details{
A dataset with variables as follows:

\itemize{
  \item RID T  Participant roster ID
  \item VISCODE T  Visit code
  \item VISCODE2 T  Translated visit code
  \item EXAMDATE D  Examination Date
  \item INFERIORCEREBELLUM_SUVR N av1451 uptake - reference region Reference region - av1451 mean of inferior cerebellar grey matter
  \item INFERIORCEREBELLUM_VOLUME N mm^3 Inferior Cerebellar Grey Matter ROI size in mm^3
  \item ERODED_SUBCORTICALWM_SUVR N av1451 uptake - reference region Reference region - av1451 mean of subcortical white matter eroded away from cortex
  \item ERODED_SUBCORTICALWM_VOLUME N mm^3 Subcortical white matter eroded away from cortex size in mm^3
  \item BRAAK1_SUVR N av1451 uptake - composite region Weighted av1451 mean of Braak region 1
  \item BRAAK1_VOLUME N mm^3 Braak region 1 size in mm^3
  \item BRAAK34_SUVR N av1451 uptake - composite region Weighted av1451 mean of Braak region 3 and 4 
  \item BRAAK34_VOLUME N mm^3 Braak region 3+4 size in mm^3
  \item META_TEMPORAL_SUVR N av1451 uptake - composite region Weighted av1451 mean of meta-temporal ROI
  \item META_TEMPORAL_VOLUME N mm^3 Meta-temporal composite region size in mm^3
  \item BRAAK56_SUVR N av1451 uptake - composite region Weighted av1451 mean of Braak region 5 and 6
  \item BRAAK56_VOLUME N mm^3 Braak region 5+6 size in mm^3
  \item BRAINSTEM_SUVR N av1451 uptake - subregion brain-stem
  \item BRAINSTEM_VOLUME N mm^3 brain-stem ROI size in mm^3
  \item VENTRICLE_3RD_SUVR N av1451 uptake - subregion 3rd-ventricle
  \item VENTRICLE_3RD_VOLUME N mm^3 3rd-ventricle ROI size in mm^3
  \item VENTRICLE_4TH_SUVR N av1451 uptake - subregion 4th-ventricle
  \item VENTRICLE_4TH_VOLUME N mm^3 4th-ventricle ROI size in mm^3
  \item VENTRICLE_5TH_SUVR N av1451 uptake - subregion 5th-ventricle
  \item VENTRICLE_5TH_VOLUME N mm^3 5th-ventricle ROI size in mm^3
  \item CC_ANTERIOR_SUVR N av1451 uptake - subregion cc-anterior
  \item CC_ANTERIOR_VOLUME N mm^3 cc-anterior ROI size in mm^3
  \item CC_CENTRAL_SUVR N av1451 uptake - subregion cc-central
  \item CC_CENTRAL_VOLUME N mm^3 cc-central ROI size in mm^3
  \item CC_MID_ANTERIOR_SUVR N av1451 uptake - subregion cc-mid-anterior
  \item CC_MID_ANTERIOR_VOLUME N mm^3 cc-mid-anterior ROI size in mm^3
  \item CC_MID_POSTERIOR_SUVR N av1451 uptake - subregion cc-mid-posterior
  \item CC_MID_POSTERIOR_VOLUME N mm^3 cc-mid-posterior ROI size in mm^3
  \item CC_POSTERIOR_SUVR N av1451 uptake - subregion cc-posterior
  \item CC_POSTERIOR_VOLUME N mm^3 cc-posterior ROI size in mm^3
  \item CSF_SUVR N av1451 uptake - subregion csf
  \item CSF_VOLUME N mm^3 csf ROI size in mm^3
  \item CTX_LH_BANKSSTS_SUVR N av1451 uptake - subregion ctx-lh-bankssts
  \item CTX_LH_BANKSSTS_VOLUME N mm^3 ctx-lh-bankssts ROI size in mm^3
  \item CTX_LH_CAUDALANTERIORCINGULATE_SUVR N av1451 uptake - subregion ctx-lh-caudalanteriorcingulate
  \item CTX_LH_CAUDALANTERIORCINGULATE_VOLUME N mm^3 ctx-lh-caudalanteriorcingulate ROI size in mm^3
  \item CTX_LH_CAUDALMIDDLEFRONTAL_SUVR N av1451 uptake - subregion ctx-lh-caudalmiddlefrontal
  \item CTX_LH_CAUDALMIDDLEFRONTAL_VOLUME N mm^3 ctx-lh-caudalmiddlefrontal ROI size in mm^3
  \item CTX_LH_CUNEUS_SUVR N av1451 uptake - subregion ctx-lh-cuneus
  \item CTX_LH_CUNEUS_VOLUME N mm^3 ctx-lh-cuneus ROI size in mm^3
  \item CTX_LH_ENTORHINAL_SUVR N av1451 uptake - subregion ctx-lh-entorhinal
  \item CTX_LH_ENTORHINAL_VOLUME N mm^3 ctx-lh-entorhinal ROI size in mm^3
  \item CTX_LH_FRONTALPOLE_SUVR N av1451 uptake - subregion ctx-lh-frontalpole
  \item CTX_LH_FRONTALPOLE_VOLUME N mm^3 ctx-lh-frontalpole ROI size in mm^3
  \item CTX_LH_FUSIFORM_SUVR N av1451 uptake - subregion ctx-lh-fusiform
  \item CTX_LH_FUSIFORM_VOLUME N mm^3 ctx-lh-fusiform ROI size in mm^3
  \item CTX_LH_INFERIORPARIETAL_SUVR N av1451 uptake - subregion ctx-lh-inferiorparietal
  \item CTX_LH_INFERIORPARIETAL_VOLUME N mm^3 ctx-lh-inferiorparietal ROI size in mm^3
  \item CTX_LH_INFERIORTEMPORAL_SUVR N av1451 uptake - subregion ctx-lh-inferiortemporal
  \item CTX_LH_INFERIORTEMPORAL_VOLUME N mm^3 ctx-lh-inferiortemporal ROI size in mm^3
  \item CTX_LH_INSULA_SUVR N av1451 uptake - subregion ctx-lh-insula
  \item CTX_LH_INSULA_VOLUME N mm^3 ctx-lh-insula ROI size in mm^3
  \item CTX_LH_ISTHMUSCINGULATE_SUVR N av1451 uptake - subregion ctx-lh-isthmuscingulate
  \item CTX_LH_ISTHMUSCINGULATE_VOLUME N mm^3 ctx-lh-isthmuscingulate ROI size in mm^3
  \item CTX_LH_LATERALOCCIPITAL_SUVR N av1451 uptake - subregion ctx-lh-lateraloccipital
  \item CTX_LH_LATERALOCCIPITAL_VOLUME N mm^3 ctx-lh-lateraloccipital ROI size in mm^3
  \item CTX_LH_LATERALORBITOFRONTAL_SUVR N av1451 uptake - subregion ctx-lh-lateralorbitofrontal
  \item CTX_LH_LATERALORBITOFRONTAL_VOLUME N mm^3 ctx-lh-lateralorbitofrontal ROI size in mm^3
  \item CTX_LH_LINGUAL_SUVR N av1451 uptake - subregion ctx-lh-lingual
  \item CTX_LH_LINGUAL_VOLUME N mm^3 ctx-lh-lingual ROI size in mm^3
  \item CTX_LH_MEDIALORBITOFRONTAL_SUVR N av1451 uptake - subregion ctx-lh-medialorbitofrontal
  \item CTX_LH_MEDIALORBITOFRONTAL_VOLUME N mm^3 ctx-lh-medialorbitofrontal ROI size in mm^3
  \item CTX_LH_MIDDLETEMPORAL_SUVR N av1451 uptake - subregion ctx-lh-middletemporal
  \item CTX_LH_MIDDLETEMPORAL_VOLUME N mm^3 ctx-lh-middletemporal ROI size in mm^3
  \item CTX_LH_PARACENTRAL_SUVR N av1451 uptake - subregion ctx-lh-paracentral
  \item CTX_LH_PARACENTRAL_VOLUME N mm^3 ctx-lh-paracentral ROI size in mm^3
  \item CTX_LH_PARAHIPPOCAMPAL_SUVR N av1451 uptake - subregion ctx-lh-parahippocampal
  \item CTX_LH_PARAHIPPOCAMPAL_VOLUME N mm^3 ctx-lh-parahippocampal ROI size in mm^3
  \item CTX_LH_PARSOPERCULARIS_SUVR N av1451 uptake - subregion ctx-lh-parsopercularis
  \item CTX_LH_PARSOPERCULARIS_VOLUME N mm^3 ctx-lh-parsopercularis ROI size in mm^3
  \item CTX_LH_PARSORBITALIS_SUVR N av1451 uptake - subregion ctx-lh-parsorbitalis
  \item CTX_LH_PARSORBITALIS_VOLUME N mm^3 ctx-lh-parsorbitalis ROI size in mm^3
  \item CTX_LH_PARSTRIANGULARIS_SUVR N av1451 uptake - subregion ctx-lh-parstriangularis
  \item CTX_LH_PARSTRIANGULARIS_VOLUME N mm^3 ctx-lh-parstriangularis ROI size in mm^3
  \item CTX_LH_PERICALCARINE_SUVR N av1451 uptake - subregion ctx-lh-pericalcarine
  \item CTX_LH_PERICALCARINE_VOLUME N mm^3 ctx-lh-pericalcarine ROI size in mm^3
  \item CTX_LH_POSTCENTRAL_SUVR N av1451 uptake - subregion ctx-lh-postcentral
  \item CTX_LH_POSTCENTRAL_VOLUME N mm^3 ctx-lh-postcentral ROI size in mm^3
  \item CTX_LH_POSTERIORCINGULATE_SUVR N av1451 uptake - subregion ctx-lh-posteriorcingulate
  \item CTX_LH_POSTERIORCINGULATE_VOLUME N mm^3 ctx-lh-posteriorcingulate ROI size in mm^3
  \item CTX_LH_PRECENTRAL_SUVR N av1451 uptake - subregion ctx-lh-precentral
  \item CTX_LH_PRECENTRAL_VOLUME N mm^3 ctx-lh-precentral ROI size in mm^3
  \item CTX_LH_PRECUNEUS_SUVR N av1451 uptake - subregion ctx-lh-precuneus
  \item CTX_LH_PRECUNEUS_VOLUME N mm^3 ctx-lh-precuneus ROI size in mm^3
  \item CTX_LH_ROSTRALANTERIORCINGULATE_SUVR N av1451 uptake - subregion ctx-lh-rostralanteriorcingulate
  \item CTX_LH_ROSTRALANTERIORCINGULATE_VOLUME N mm^3 ctx-lh-rostralanteriorcingulate ROI size in mm^3
  \item CTX_LH_ROSTRALMIDDLEFRONTAL_SUVR N av1451 uptake - subregion ctx-lh-rostralmiddlefrontal
  \item CTX_LH_ROSTRALMIDDLEFRONTAL_VOLUME N mm^3 ctx-lh-rostralmiddlefrontal ROI size in mm^3
  \item CTX_LH_SUPERIORFRONTAL_SUVR N av1451 uptake - subregion ctx-lh-superiorfrontal
  \item CTX_LH_SUPERIORFRONTAL_VOLUME N mm^3 ctx-lh-superiorfrontal ROI size in mm^3
  \item CTX_LH_SUPERIORPARIETAL_SUVR N av1451 uptake - subregion ctx-lh-superiorparietal
  \item CTX_LH_SUPERIORPARIETAL_VOLUME N mm^3 ctx-lh-superiorparietal ROI size in mm^3
  \item CTX_LH_SUPERIORTEMPORAL_SUVR N av1451 uptake - subregion ctx-lh-superiortemporal
  \item CTX_LH_SUPERIORTEMPORAL_VOLUME N mm^3 ctx-lh-superiortemporal ROI size in mm^3
  \item CTX_LH_SUPRAMARGINAL_SUVR N av1451 uptake - subregion ctx-lh-supramarginal
  \item CTX_LH_SUPRAMARGINAL_VOLUME N mm^3 ctx-lh-supramarginal ROI size in mm^3
  \item CTX_LH_TEMPORALPOLE_SUVR N av1451 uptake - subregion ctx-lh-temporalpole
  \item CTX_LH_TEMPORALPOLE_VOLUME N mm^3 ctx-lh-temporalpole ROI size in mm^3
  \item CTX_LH_TRANSVERSETEMPORAL_SUVR N av1451 uptake - subregion ctx-lh-transversetemporal
  \item CTX_LH_TRANSVERSETEMPORAL_VOLUME N mm^3 ctx-lh-transversetemporal ROI size in mm^3
  \item CTX_LH_UNKNOWN_SUVR N av1451 uptake - subregion ctx-lh-unknown
  \item CTX_LH_UNKNOWN_VOLUME N mm^3 ctx-lh-unknown ROI size in mm^3
  \item CTX_RH_BANKSSTS_SUVR N av1451 uptake - subregion ctx-rh-bankssts
  \item CTX_RH_BANKSSTS_VOLUME N mm^3 ctx-rh-bankssts ROI size in mm^3
  \item CTX_RH_CAUDALANTERIORCINGULATE_SUVR N av1451 uptake - subregion ctx-rh-caudalanteriorcingulate
  \item CTX_RH_CAUDALANTERIORCINGULATE_VOLUME N mm^3 ctx-rh-caudalanteriorcingulate ROI size in mm^3
  \item CTX_RH_CAUDALMIDDLEFRONTAL_SUVR N av1451 uptake - subregion ctx-rh-caudalmiddlefrontal
  \item CTX_RH_CAUDALMIDDLEFRONTAL_VOLUME N mm^3 ctx-rh-caudalmiddlefrontal ROI size in mm^3
  \item CTX_RH_CUNEUS_SUVR N av1451 uptake - subregion ctx-rh-cuneus
  \item CTX_RH_CUNEUS_VOLUME N mm^3 ctx-rh-cuneus ROI size in mm^3
  \item CTX_RH_ENTORHINAL_SUVR N av1451 uptake - subregion ctx-rh-entorhinal
  \item CTX_RH_ENTORHINAL_VOLUME N mm^3 ctx-rh-entorhinal ROI size in mm^3
  \item CTX_RH_FRONTALPOLE_SUVR N av1451 uptake - subregion ctx-rh-frontalpole
  \item CTX_RH_FRONTALPOLE_VOLUME N mm^3 ctx-rh-frontalpole ROI size in mm^3
  \item CTX_RH_FUSIFORM_SUVR N av1451 uptake - subregion ctx-rh-fusiform
  \item CTX_RH_FUSIFORM_VOLUME N mm^3 ctx-rh-fusiform ROI size in mm^3
  \item CTX_RH_INFERIORPARIETAL_SUVR N av1451 uptake - subregion ctx-rh-inferiorparietal
  \item CTX_RH_INFERIORPARIETAL_VOLUME N mm^3 ctx-rh-inferiorparietal ROI size in mm^3
  \item CTX_RH_INFERIORTEMPORAL_SUVR N av1451 uptake - subregion ctx-rh-inferiortemporal
  \item CTX_RH_INFERIORTEMPORAL_VOLUME N mm^3 ctx-rh-inferiortemporal ROI size in mm^3
  \item CTX_RH_INSULA_SUVR N av1451 uptake - subregion ctx-rh-insula
  \item CTX_RH_INSULA_VOLUME N mm^3 ctx-rh-insula ROI size in mm^3
  \item CTX_RH_ISTHMUSCINGULATE_SUVR N av1451 uptake - subregion ctx-rh-isthmuscingulate
  \item CTX_RH_ISTHMUSCINGULATE_VOLUME N mm^3 ctx-rh-isthmuscingulate ROI size in mm^3
  \item CTX_RH_LATERALOCCIPITAL_SUVR N av1451 uptake - subregion ctx-rh-lateraloccipital
  \item CTX_RH_LATERALOCCIPITAL_VOLUME N mm^3 ctx-rh-lateraloccipital ROI size in mm^3
  \item CTX_RH_LATERALORBITOFRONTAL_SUVR N av1451 uptake - subregion ctx-rh-lateralorbitofrontal
  \item CTX_RH_LATERALORBITOFRONTAL_VOLUME N mm^3 ctx-rh-lateralorbitofrontal ROI size in mm^3
  \item CTX_RH_LINGUAL_SUVR N av1451 uptake - subregion ctx-rh-lingual
  \item CTX_RH_LINGUAL_VOLUME N mm^3 ctx-rh-lingual ROI size in mm^3
  \item CTX_RH_MEDIALORBITOFRONTAL_SUVR N av1451 uptake - subregion ctx-rh-medialorbitofrontal
  \item CTX_RH_MEDIALORBITOFRONTAL_VOLUME N mm^3 ctx-rh-medialorbitofrontal ROI size in mm^3
  \item CTX_RH_MIDDLETEMPORAL_SUVR N av1451 uptake - subregion ctx-rh-middletemporal
  \item CTX_RH_MIDDLETEMPORAL_VOLUME N mm^3 ctx-rh-middletemporal ROI size in mm^3
  \item CTX_RH_PARACENTRAL_SUVR N av1451 uptake - subregion ctx-rh-paracentral
  \item CTX_RH_PARACENTRAL_VOLUME N mm^3 ctx-rh-paracentral ROI size in mm^3
  \item CTX_RH_PARAHIPPOCAMPAL_SUVR N av1451 uptake - subregion ctx-rh-parahippocampal
  \item CTX_RH_PARAHIPPOCAMPAL_VOLUME N mm^3 ctx-rh-parahippocampal ROI size in mm^3
  \item CTX_RH_PARSOPERCULARIS_SUVR N av1451 uptake - subregion ctx-rh-parsopercularis
  \item CTX_RH_PARSOPERCULARIS_VOLUME N mm^3 ctx-rh-parsopercularis ROI size in mm^3
  \item CTX_RH_PARSORBITALIS_SUVR N av1451 uptake - subregion ctx-rh-parsorbitalis
  \item CTX_RH_PARSORBITALIS_VOLUME N mm^3 ctx-rh-parsorbitalis ROI size in mm^3
  \item CTX_RH_PARSTRIANGULARIS_SUVR N av1451 uptake - subregion ctx-rh-parstriangularis
  \item CTX_RH_PARSTRIANGULARIS_VOLUME N mm^3 ctx-rh-parstriangularis ROI size in mm^3
  \item CTX_RH_PERICALCARINE_SUVR N av1451 uptake - subregion ctx-rh-pericalcarine
  \item CTX_RH_PERICALCARINE_VOLUME N mm^3 ctx-rh-pericalcarine ROI size in mm^3
  \item CTX_RH_POSTCENTRAL_SUVR N av1451 uptake - subregion ctx-rh-postcentral
  \item CTX_RH_POSTCENTRAL_VOLUME N mm^3 ctx-rh-postcentral ROI size in mm^3
  \item CTX_RH_POSTERIORCINGULATE_SUVR N av1451 uptake - subregion ctx-rh-posteriorcingulate
  \item CTX_RH_POSTERIORCINGULATE_VOLUME N mm^3 ctx-rh-posteriorcingulate ROI size in mm^3
  \item CTX_RH_PRECENTRAL_SUVR N av1451 uptake - subregion ctx-rh-precentral
  \item CTX_RH_PRECENTRAL_VOLUME N mm^3 ctx-rh-precentral ROI size in mm^3
  \item CTX_RH_PRECUNEUS_SUVR N av1451 uptake - subregion ctx-rh-precuneus
  \item CTX_RH_PRECUNEUS_VOLUME N mm^3 ctx-rh-precuneus ROI size in mm^3
  \item CTX_RH_ROSTRALANTERIORCINGULATE_SUVR N av1451 uptake - subregion ctx-rh-rostralanteriorcingulate
  \item CTX_RH_ROSTRALANTERIORCINGULATE_VOLUME N mm^3 ctx-rh-rostralanteriorcingulate ROI size in mm^3
  \item CTX_RH_ROSTRALMIDDLEFRONTAL_SUVR N av1451 uptake - subregion ctx-rh-rostralmiddlefrontal
  \item CTX_RH_ROSTRALMIDDLEFRONTAL_VOLUME N mm^3 ctx-rh-rostralmiddlefrontal ROI size in mm^3
  \item CTX_RH_SUPERIORFRONTAL_SUVR N av1451 uptake - subregion ctx-rh-superiorfrontal
  \item CTX_RH_SUPERIORFRONTAL_VOLUME N mm^3 ctx-rh-superiorfrontal ROI size in mm^3
  \item CTX_RH_SUPERIORPARIETAL_SUVR N av1451 uptake - subregion ctx-rh-superiorparietal
  \item CTX_RH_SUPERIORPARIETAL_VOLUME N mm^3 ctx-rh-superiorparietal ROI size in mm^3
  \item CTX_RH_SUPERIORTEMPORAL_SUVR N av1451 uptake - subregion ctx-rh-superiortemporal
  \item CTX_RH_SUPERIORTEMPORAL_VOLUME N mm^3 ctx-rh-superiortemporal ROI size in mm^3
  \item CTX_RH_SUPRAMARGINAL_SUVR N av1451 uptake - subregion ctx-rh-supramarginal
  \item CTX_RH_SUPRAMARGINAL_VOLUME N mm^3 ctx-rh-supramarginal ROI size in mm^3
  \item CTX_RH_TEMPORALPOLE_SUVR N av1451 uptake - subregion ctx-rh-temporalpole
  \item CTX_RH_TEMPORALPOLE_VOLUME N mm^3 ctx-rh-temporalpole ROI size in mm^3
  \item CTX_RH_TRANSVERSETEMPORAL_SUVR N av1451 uptake - subregion ctx-rh-transversetemporal
  \item CTX_RH_TRANSVERSETEMPORAL_VOLUME N mm^3 ctx-rh-transversetemporal ROI size in mm^3
  \item CTX_RH_UNKNOWN_SUVR N av1451 uptake - subregion ctx-rh-unknown
  \item CTX_RH_UNKNOWN_VOLUME N mm^3 ctx-rh-unknown ROI size in mm^3
  \item LEFT_ACCUMBENS_AREA_SUVR N av1451 uptake - subregion left-accumbens-area
  \item LEFT_ACCUMBENS_AREA_VOLUME N mm^3 left-accumbens-area ROI size in mm^3
  \item LEFT_AMYGDALA_SUVR N av1451 uptake - subregion left-amygdala
  \item LEFT_AMYGDALA_VOLUME N mm^3 left-amygdala ROI size in mm^3
  \item LEFT_CAUDATE_SUVR N av1451 uptake - subregion left-caudate
  \item LEFT_CAUDATE_VOLUME N mm^3 left-caudate ROI size in mm^3
  \item LEFT_CEREBELLUM_CORTEX_SUVR N av1451 uptake - subregion left-cerebellum-cortex
  \item LEFT_CEREBELLUM_CORTEX_VOLUME N mm^3 left-cerebellum-cortex ROI size in mm^3
  \item LEFT_CEREBELLUM_WHITE_MATTER_SUVR N av1451 uptake - subregion left-cerebellum-white-matter
  \item LEFT_CEREBELLUM_WHITE_MATTER_VOLUME N mm^3 left-cerebellum-white-matter ROI size in mm^3
  \item LEFT_CEREBRAL_WHITE_MATTER_SUVR N av1451 uptake - subregion left-cerebral-white-matter
  \item LEFT_CEREBRAL_WHITE_MATTER_VOLUME N mm^3 left-cerebral-white-matter ROI size in mm^3
  \item LEFT_CHOROID_PLEXUS_SUVR N av1451 uptake - subregion left-choroid-plexus
  \item LEFT_CHOROID_PLEXUS_VOLUME N mm^3 left-choroid-plexus ROI size in mm^3
  \item LEFT_HIPPOCAMPUS_SUVR N av1451 uptake - subregion left-hippocampus
  \item LEFT_HIPPOCAMPUS_VOLUME N mm^3 left-hippocampus ROI size in mm^3
  \item LEFT_INF_LAT_VENT_SUVR N av1451 uptake - subregion left-inf-lat-vent
  \item LEFT_INF_LAT_VENT_VOLUME N mm^3 left-inf-lat-vent ROI size in mm^3
  \item LEFT_LATERAL_VENTRICLE_SUVR N av1451 uptake - subregion left-lateral-ventricle
  \item LEFT_LATERAL_VENTRICLE_VOLUME N mm^3 left-lateral-ventricle ROI size in mm^3
  \item LEFT_PALLIDUM_SUVR N av1451 uptake - subregion left-pallidum
  \item LEFT_PALLIDUM_VOLUME N mm^3 left-pallidum ROI size in mm^3
  \item LEFT_PUTAMEN_SUVR N av1451 uptake - subregion left-putamen
  \item LEFT_PUTAMEN_VOLUME N mm^3 left-putamen ROI size in mm^3
  \item LEFT_THALAMUS_PROPER_SUVR N av1451 uptake - subregion left-thalamus-proper
  \item LEFT_THALAMUS_PROPER_VOLUME N mm^3 left-thalamus-proper ROI size in mm^3
  \item LEFT_VENTRALDC_SUVR N av1451 uptake - subregion left-ventraldc
  \item LEFT_VENTRALDC_VOLUME N mm^3 left-ventraldc ROI size in mm^3
  \item LEFT_VESSEL_SUVR N av1451 uptake - subregion left-vessel
  \item LEFT_VESSEL_VOLUME N mm^3 left-vessel ROI size in mm^3
  \item NON_WM_HYPOINTENSITIES_SUVR N av1451 uptake - subregion non-wm-hypointensities
  \item NON_WM_HYPOINTENSITIES_VOLUME N mm^3 non-wm-hypointensities ROI size in mm^3
  \item OPTIC_CHIASM_SUVR N av1451 uptake - subregion optic-chiasm
  \item OPTIC_CHIASM_VOLUME N mm^3 optic-chiasm ROI size in mm^3
  \item RIGHT_ACCUMBENS_AREA_SUVR N av1451 uptake - subregion right-accumbens-area
  \item RIGHT_ACCUMBENS_AREA_VOLUME N mm^3 right-accumbens-area ROI size in mm^3
  \item RIGHT_AMYGDALA_SUVR N av1451 uptake - subregion right-amygdala
  \item RIGHT_AMYGDALA_VOLUME N mm^3 right-amygdala ROI size in mm^3
  \item RIGHT_CAUDATE_SUVR N av1451 uptake - subregion right-caudate
  \item RIGHT_CAUDATE_VOLUME N mm^3 right-caudate ROI size in mm^3
  \item RIGHT_CEREBELLUM_CORTEX_SUVR N av1451 uptake - subregion right-cerebellum-cortex
  \item RIGHT_CEREBELLUM_CORTEX_VOLUME N mm^3 right-cerebellum-cortex ROI size in mm^3
  \item RIGHT_CEREBELLUM_WHITE_MATTER_SUVR N av1451 uptake - subregion right-cerebellum-white-matter
  \item RIGHT_CEREBELLUM_WHITE_MATTER_VOLUME N mm^3 right-cerebellum-white-matter ROI size in mm^3
  \item RIGHT_CEREBRAL_WHITE_MATTER_SUVR N av1451 uptake - subregion right-cerebral-white-matter
  \item RIGHT_CEREBRAL_WHITE_MATTER_VOLUME N mm^3 right-cerebral-white-matter ROI size in mm^3
  \item RIGHT_CHOROID_PLEXUS_SUVR N av1451 uptake - subregion right-choroid-plexus
  \item RIGHT_CHOROID_PLEXUS_VOLUME N mm^3 right-choroid-plexus ROI size in mm^3
  \item RIGHT_HIPPOCAMPUS_SUVR N av1451 uptake - subregion right-hippocampus
  \item RIGHT_HIPPOCAMPUS_VOLUME N mm^3 right-hippocampus ROI size in mm^3
  \item RIGHT_INF_LAT_VENT_SUVR N av1451 uptake - subregion right-inf-lat-vent
  \item RIGHT_INF_LAT_VENT_VOLUME N mm^3 right-inf-lat-vent ROI size in mm^3
  \item RIGHT_LATERAL_VENTRICLE_SUVR N av1451 uptake - subregion right-lateral-ventricle
  \item RIGHT_LATERAL_VENTRICLE_VOLUME N mm^3 right-lateral-ventricle ROI size in mm^3
  \item RIGHT_PALLIDUM_SUVR N av1451 uptake - subregion right-pallidum
  \item RIGHT_PALLIDUM_VOLUME N mm^3 right-pallidum ROI size in mm^3
  \item RIGHT_PUTAMEN_SUVR N av1451 uptake - subregion right-putamen
  \item RIGHT_PUTAMEN_VOLUME N mm^3 right-putamen ROI size in mm^3
  \item RIGHT_THALAMUS_PROPER_SUVR N av1451 uptake - subregion right-thalamus-proper
  \item RIGHT_THALAMUS_PROPER_VOLUME N mm^3 right-thalamus-proper ROI size in mm^3
  \item RIGHT_VENTRALDC_SUVR N av1451 uptake - subregion right-ventraldc
  \item RIGHT_VENTRALDC_VOLUME N mm^3 right-ventraldc ROI size in mm^3
  \item RIGHT_VESSEL_SUVR N av1451 uptake - subregion right-vessel
  \item RIGHT_VESSEL_VOLUME N mm^3 right-vessel ROI size in mm^3
  \item WM_HYPOINTENSITIES_SUVR N av1451 uptake - subregion wm-hypointensities
  \item WM_HYPOINTENSITIES_VOLUME N mm^3 wm-hypointensities ROI size in mm^3
  \item CEREBELLUMGREYMATTER N av1451 uptake - reference region Cerebellum grey matter
  \item CEREBELLUMGREYMATTER_SIZE N mm^3 Cerebellum grey matter ROI size in mm^3
  \item BRAAK1 N av1451 uptake - reference region Braak region 1
  \item BRAAK1_SIZE N mm^3 Braak region 1 ROI size in mm^3
  \item BRAAK2 N av1451 uptake - subregion Braak region 2
  \item BRAAK2_SIZE N mm^3 Braak region 2 ROI size in mm^3
  \item BRAAK3 N av1451 uptake - subregion Braak region 3
  \item BRAAK3_SIZE N mm^3 Braak region 3 ROI size in mm^3
  \item BRAAK4 N av1451 uptake - subregion Braak region 4
  \item BRAAK4_SIZE N mm^3 Braak region 4 ROI size in mm^3
  \item BRAAK5 N av1451 uptake - subregion Braak region 5
  \item BRAAK5_SIZE N mm^3 Braak region 5 ROI size in mm^3
  \item BRAAK6 N av1451 uptake - subregion Braak region 6
  \item BRAAK6_SIZE N mm^3 Braak region 6 ROI size in mm^3
  \item ERODED_SUBCORTICALWM N av1451 uptake - subregion subcortical WM eroded away from cortex
  \item ERODED_SUBCORTICALWM_SIZE N mm^3 subcortical WM eroded away from cortex in mm^3
  \item BRAINSTEM N av1451 uptake - subregion brain-stem
  \item BRAINSTEM_SIZE N mm^3 brain-stem ROI size in mm^3
  \item VENTRICLE_3RD N av1451 uptake - subregion 3rd-ventricle
  \item VENTRICLE_3RD_SIZE N mm^3 3rd-ventricle ROI size in mm^3
  \item VENTRICLE_4TH N av1451 uptake - subregion 4th-ventricle
  \item VENTRICLE_4TH_SIZE N mm^3 4th-ventricle ROI size in mm^3
  \item VENTRICLE_5TH N av1451 uptake - subregion 5th-ventricle
  \item VENTRICLE_5TH_SIZE N mm^3 5th-ventricle ROI size in mm^3
  \item CC_ANTERIOR N av1451 uptake - subregion cc-anterior
  \item CC_ANTERIOR_SIZE N mm^3 cc-anterior ROI size in mm^3
  \item CC_CENTRAL N av1451 uptake - subregion cc-central
  \item CC_CENTRAL_SIZE N mm^3 cc-central ROI size in mm^3
  \item CC_MID_ANTERIOR N av1451 uptake - subregion cc-mid-anterior
  \item CC_MID_ANTERIOR_SIZE N mm^3 cc-mid-anterior ROI size in mm^3
  \item CC_MID_POSTERIOR N av1451 uptake - subregion cc-mid-posterior
  \item CC_MID_POSTERIOR_SIZE N mm^3 cc-mid-posterior ROI size in mm^3
  \item CC_POSTERIOR N av1451 uptake - subregion cc-posterior
  \item CC_POSTERIOR_SIZE N mm^3 cc-posterior ROI size in mm^3
  \item CSF N av1451 uptake - subregion csf
  \item CSF_SIZE N mm^3 csf ROI size in mm^3
  \item CTX_LH_BANKSSTS N av1451 uptake - subregion ctx-lh-bankssts
  \item CTX_LH_BANKSSTS_SIZE N mm^3 ctx-lh-bankssts ROI size in mm^3
  \item CTX_LH_CAUDALANTERIORCINGULATE N av1451 uptake - subregion ctx-lh-caudalanteriorcingulate
  \item CTX_LH_CAUDALANTERIORCINGULATE_SIZE N mm^3 ctx-lh-caudalanteriorcingulate ROI size in mm^3
  \item CTX_LH_CAUDALMIDDLEFRONTAL N av1451 uptake - subregion ctx-lh-caudalmiddlefrontal
  \item CTX_LH_CAUDALMIDDLEFRONTAL_SIZE N mm^3 ctx-lh-caudalmiddlefrontal ROI size in mm^3
  \item CTX_LH_CUNEUS N av1451 uptake - subregion ctx-lh-cuneus
  \item CTX_LH_CUNEUS_SIZE N mm^3 ctx-lh-cuneus ROI size in mm^3
  \item CTX_LH_ENTORHINAL N av1451 uptake - subregion ctx-lh-entorhinal
  \item CTX_LH_ENTORHINAL_SIZE N mm^3 ctx-lh-entorhinal ROI size in mm^3
  \item CTX_LH_FRONTALPOLE N av1451 uptake - subregion ctx-lh-frontalpole
  \item CTX_LH_FRONTALPOLE_SIZE N mm^3 ctx-lh-frontalpole ROI size in mm^3
  \item CTX_LH_FUSIFORM N av1451 uptake - subregion ctx-lh-fusiform
  \item CTX_LH_FUSIFORM_SIZE N mm^3 ctx-lh-fusiform ROI size in mm^3
  \item CTX_LH_INFERIORPARIETAL N av1451 uptake - subregion ctx-lh-inferiorparietal
  \item CTX_LH_INFERIORPARIETAL_SIZE N mm^3 ctx-lh-inferiorparietal ROI size in mm^3
  \item CTX_LH_INFERIORTEMPORAL N av1451 uptake - subregion ctx-lh-inferiortemporal
  \item CTX_LH_INFERIORTEMPORAL_SIZE N mm^3 ctx-lh-inferiortemporal ROI size in mm^3
  \item CTX_LH_INSULA N av1451 uptake - subregion ctx-lh-insula
  \item CTX_LH_INSULA_SIZE N mm^3 ctx-lh-insula ROI size in mm^3
  \item CTX_LH_ISTHMUSCINGULATE N av1451 uptake - subregion ctx-lh-isthmuscingulate
  \item CTX_LH_ISTHMUSCINGULATE_SIZE N mm^3 ctx-lh-isthmuscingulate ROI size in mm^3
  \item CTX_LH_LATERALOCCIPITAL N av1451 uptake - subregion ctx-lh-lateraloccipital
  \item CTX_LH_LATERALOCCIPITAL_SIZE N mm^3 ctx-lh-lateraloccipital ROI size in mm^3
  \item CTX_LH_LATERALORBITOFRONTAL N av1451 uptake - subregion ctx-lh-lateralorbitofrontal
  \item CTX_LH_LATERALORBITOFRONTAL_SIZE N mm^3 ctx-lh-lateralorbitofrontal ROI size in mm^3
  \item CTX_LH_LINGUAL N av1451 uptake - subregion ctx-lh-lingual
  \item CTX_LH_LINGUAL_SIZE N mm^3 ctx-lh-lingual ROI size in mm^3
  \item CTX_LH_MEDIALORBITOFRONTAL N av1451 uptake - subregion ctx-lh-medialorbitofrontal
  \item CTX_LH_MEDIALORBITOFRONTAL_SIZE N mm^3 ctx-lh-medialorbitofrontal ROI size in mm^3
  \item CTX_LH_MIDDLETEMPORAL N av1451 uptake - subregion ctx-lh-middletemporal
  \item CTX_LH_MIDDLETEMPORAL_SIZE N mm^3 ctx-lh-middletemporal ROI size in mm^3
  \item CTX_LH_PARACENTRAL N av1451 uptake - subregion ctx-lh-paracentral
  \item CTX_LH_PARACENTRAL_SIZE N mm^3 ctx-lh-paracentral ROI size in mm^3
  \item CTX_LH_PARAHIPPOCAMPAL N av1451 uptake - subregion ctx-lh-parahippocampal
  \item CTX_LH_PARAHIPPOCAMPAL_SIZE N mm^3 ctx-lh-parahippocampal ROI size in mm^3
  \item CTX_LH_PARSOPERCULARIS N av1451 uptake - subregion ctx-lh-parsopercularis
  \item CTX_LH_PARSOPERCULARIS_SIZE N mm^3 ctx-lh-parsopercularis ROI size in mm^3
  \item CTX_LH_PARSORBITALIS N av1451 uptake - subregion ctx-lh-parsorbitalis
  \item CTX_LH_PARSORBITALIS_SIZE N mm^3 ctx-lh-parsorbitalis ROI size in mm^3
  \item CTX_LH_PARSTRIANGULARIS N av1451 uptake - subregion ctx-lh-parstriangularis
  \item CTX_LH_PARSTRIANGULARIS_SIZE N mm^3 ctx-lh-parstriangularis ROI size in mm^3
  \item CTX_LH_PERICALCARINE N av1451 uptake - subregion ctx-lh-pericalcarine
  \item CTX_LH_PERICALCARINE_SIZE N mm^3 ctx-lh-pericalcarine ROI size in mm^3
  \item CTX_LH_POSTCENTRAL N av1451 uptake - subregion ctx-lh-postcentral
  \item CTX_LH_POSTCENTRAL_SIZE N mm^3 ctx-lh-postcentral ROI size in mm^3
  \item CTX_LH_POSTERIORCINGULATE N av1451 uptake - subregion ctx-lh-posteriorcingulate
  \item CTX_LH_POSTERIORCINGULATE_SIZE N mm^3 ctx-lh-posteriorcingulate ROI size in mm^3
  \item CTX_LH_PRECENTRAL N av1451 uptake - subregion ctx-lh-precentral
  \item CTX_LH_PRECENTRAL_SIZE N mm^3 ctx-lh-precentral ROI size in mm^3
  \item CTX_LH_PRECUNEUS N av1451 uptake - subregion ctx-lh-precuneus
  \item CTX_LH_PRECUNEUS_SIZE N mm^3 ctx-lh-precuneus ROI size in mm^3
  \item CTX_LH_ROSTRALANTERIORCINGULATE N av1451 uptake - subregion ctx-lh-rostralanteriorcingulate
  \item CTX_LH_ROSTRALANTERIORCINGULATE_SIZE N mm^3 ctx-lh-rostralanteriorcingulate ROI size in mm^3
  \item CTX_LH_ROSTRALMIDDLEFRONTAL N av1451 uptake - subregion ctx-lh-rostralmiddlefrontal
  \item CTX_LH_ROSTRALMIDDLEFRONTAL_SIZE N mm^3 ctx-lh-rostralmiddlefrontal ROI size in mm^3
  \item CTX_LH_SUPERIORFRONTAL N av1451 uptake - subregion ctx-lh-superiorfrontal
  \item CTX_LH_SUPERIORFRONTAL_SIZE N mm^3 ctx-lh-superiorfrontal ROI size in mm^3
  \item CTX_LH_SUPERIORPARIETAL N av1451 uptake - subregion ctx-lh-superiorparietal
  \item CTX_LH_SUPERIORPARIETAL_SIZE N mm^3 ctx-lh-superiorparietal ROI size in mm^3
  \item CTX_LH_SUPERIORTEMPORAL N av1451 uptake - subregion ctx-lh-superiortemporal
  \item CTX_LH_SUPERIORTEMPORAL_SIZE N mm^3 ctx-lh-superiortemporal ROI size in mm^3
  \item CTX_LH_SUPRAMARGINAL N av1451 uptake - subregion ctx-lh-supramarginal
  \item CTX_LH_SUPRAMARGINAL_SIZE N mm^3 ctx-lh-supramarginal ROI size in mm^3
  \item CTX_LH_TEMPORALPOLE N av1451 uptake - subregion ctx-lh-temporalpole
  \item CTX_LH_TEMPORALPOLE_SIZE N mm^3 ctx-lh-temporalpole ROI size in mm^3
  \item CTX_LH_TRANSVERSETEMPORAL N av1451 uptake - subregion ctx-lh-transversetemporal
  \item CTX_LH_TRANSVERSETEMPORAL_SIZE N mm^3 ctx-lh-transversetemporal ROI size in mm^3
  \item CTX_LH_UNKNOWN N av1451 uptake - subregion ctx-lh-unknown
  \item CTX_LH_UNKNOWN_SIZE N mm^3 ctx-lh-unknown ROI size in mm^3
  \item CTX_RH_BANKSSTS N av1451 uptake - subregion ctx-rh-bankssts
  \item CTX_RH_BANKSSTS_SIZE N mm^3 ctx-rh-bankssts ROI size in mm^3
  \item CTX_RH_CAUDALANTERIORCINGULATE N av1451 uptake - subregion ctx-rh-caudalanteriorcingulate
  \item CTX_RH_CAUDALANTERIORCINGULATE_SIZE N mm^3 ctx-rh-caudalanteriorcingulate ROI size in mm^3
  \item CTX_RH_CAUDALMIDDLEFRONTAL N av1451 uptake - subregion ctx-rh-caudalmiddlefrontal
  \item CTX_RH_CAUDALMIDDLEFRONTAL_SIZE N mm^3 ctx-rh-caudalmiddlefrontal ROI size in mm^3
  \item CTX_RH_CUNEUS N av1451 uptake - subregion ctx-rh-cuneus
  \item CTX_RH_CUNEUS_SIZE N mm^3 ctx-rh-cuneus ROI size in mm^3
  \item CTX_RH_ENTORHINAL N av1451 uptake - subregion ctx-rh-entorhinal
  \item CTX_RH_ENTORHINAL_SIZE N mm^3 ctx-rh-entorhinal ROI size in mm^3
  \item CTX_RH_FRONTALPOLE N av1451 uptake - subregion ctx-rh-frontalpole
  \item CTX_RH_FRONTALPOLE_SIZE N mm^3 ctx-rh-frontalpole ROI size in mm^3
  \item CTX_RH_FUSIFORM N av1451 uptake - subregion ctx-rh-fusiform
  \item CTX_RH_FUSIFORM_SIZE N mm^3 ctx-rh-fusiform ROI size in mm^3
  \item CTX_RH_INFERIORPARIETAL N av1451 uptake - subregion ctx-rh-inferiorparietal
  \item CTX_RH_INFERIORPARIETAL_SIZE N mm^3 ctx-rh-inferiorparietal ROI size in mm^3
  \item CTX_RH_INFERIORTEMPORAL N av1451 uptake - subregion ctx-rh-inferiortemporal
  \item CTX_RH_INFERIORTEMPORAL_SIZE N mm^3 ctx-rh-inferiortemporal ROI size in mm^3
  \item CTX_RH_INSULA N av1451 uptake - subregion ctx-rh-insula
  \item CTX_RH_INSULA_SIZE N mm^3 ctx-rh-insula ROI size in mm^3
  \item CTX_RH_ISTHMUSCINGULATE N av1451 uptake - subregion ctx-rh-isthmuscingulate
  \item CTX_RH_ISTHMUSCINGULATE_SIZE N mm^3 ctx-rh-isthmuscingulate ROI size in mm^3
  \item CTX_RH_LATERALOCCIPITAL N av1451 uptake - subregion ctx-rh-lateraloccipital
  \item CTX_RH_LATERALOCCIPITAL_SIZE N mm^3 ctx-rh-lateraloccipital ROI size in mm^3
  \item CTX_RH_LATERALORBITOFRONTAL N av1451 uptake - subregion ctx-rh-lateralorbitofrontal
  \item CTX_RH_LATERALORBITOFRONTAL_SIZE N mm^3 ctx-rh-lateralorbitofrontal ROI size in mm^3
  \item CTX_RH_LINGUAL N av1451 uptake - subregion ctx-rh-lingual
  \item CTX_RH_LINGUAL_SIZE N mm^3 ctx-rh-lingual ROI size in mm^3
  \item CTX_RH_MEDIALORBITOFRONTAL N av1451 uptake - subregion ctx-rh-medialorbitofrontal
  \item CTX_RH_MEDIALORBITOFRONTAL_SIZE N mm^3 ctx-rh-medialorbitofrontal ROI size in mm^3
  \item CTX_RH_MIDDLETEMPORAL N av1451 uptake - subregion ctx-rh-middletemporal
  \item CTX_RH_MIDDLETEMPORAL_SIZE N mm^3 ctx-rh-middletemporal ROI size in mm^3
  \item CTX_RH_PARACENTRAL N av1451 uptake - subregion ctx-rh-paracentral
  \item CTX_RH_PARACENTRAL_SIZE N mm^3 ctx-rh-paracentral ROI size in mm^3
  \item CTX_RH_PARAHIPPOCAMPAL N av1451 uptake - subregion ctx-rh-parahippocampal
  \item CTX_RH_PARAHIPPOCAMPAL_SIZE N mm^3 ctx-rh-parahippocampal ROI size in mm^3
  \item CTX_RH_PARSOPERCULARIS N av1451 uptake - subregion ctx-rh-parsopercularis
  \item CTX_RH_PARSOPERCULARIS_SIZE N mm^3 ctx-rh-parsopercularis ROI size in mm^3
  \item CTX_RH_PARSORBITALIS N av1451 uptake - subregion ctx-rh-parsorbitalis
  \item CTX_RH_PARSORBITALIS_SIZE N mm^3 ctx-rh-parsorbitalis ROI size in mm^3
  \item CTX_RH_PARSTRIANGULARIS N av1451 uptake - subregion ctx-rh-parstriangularis
  \item CTX_RH_PARSTRIANGULARIS_SIZE N mm^3 ctx-rh-parstriangularis ROI size in mm^3
  \item CTX_RH_PERICALCARINE N av1451 uptake - subregion ctx-rh-pericalcarine
  \item CTX_RH_PERICALCARINE_SIZE N mm^3 ctx-rh-pericalcarine ROI size in mm^3
  \item CTX_RH_POSTCENTRAL N av1451 uptake - subregion ctx-rh-postcentral
  \item CTX_RH_POSTCENTRAL_SIZE N mm^3 ctx-rh-postcentral ROI size in mm^3
  \item CTX_RH_POSTERIORCINGULATE N av1451 uptake - subregion ctx-rh-posteriorcingulate
  \item CTX_RH_POSTERIORCINGULATE_SIZE N mm^3 ctx-rh-posteriorcingulate ROI size in mm^3
  \item CTX_RH_PRECENTRAL N av1451 uptake - subregion ctx-rh-precentral
  \item CTX_RH_PRECENTRAL_SIZE N mm^3 ctx-rh-precentral ROI size in mm^3
  \item CTX_RH_PRECUNEUS N av1451 uptake - subregion ctx-rh-precuneus
  \item CTX_RH_PRECUNEUS_SIZE N mm^3 ctx-rh-precuneus ROI size in mm^3
  \item CTX_RH_ROSTRALANTERIORCINGULATE N av1451 uptake - subregion ctx-rh-rostralanteriorcingulate
  \item CTX_RH_ROSTRALANTERIORCINGULATE_SIZE N mm^3 ctx-rh-rostralanteriorcingulate ROI size in mm^3
  \item CTX_RH_ROSTRALMIDDLEFRONTAL N av1451 uptake - subregion ctx-rh-rostralmiddlefrontal
  \item CTX_RH_ROSTRALMIDDLEFRONTAL_SIZE N mm^3 ctx-rh-rostralmiddlefrontal ROI size in mm^3
  \item CTX_RH_SUPERIORFRONTAL N av1451 uptake - subregion ctx-rh-superiorfrontal
  \item CTX_RH_SUPERIORFRONTAL_SIZE N mm^3 ctx-rh-superiorfrontal ROI size in mm^3
  \item CTX_RH_SUPERIORPARIETAL N av1451 uptake - subregion ctx-rh-superiorparietal
  \item CTX_RH_SUPERIORPARIETAL_SIZE N mm^3 ctx-rh-superiorparietal ROI size in mm^3
  \item CTX_RH_SUPERIORTEMPORAL N av1451 uptake - subregion ctx-rh-superiortemporal
  \item CTX_RH_SUPERIORTEMPORAL_SIZE N mm^3 ctx-rh-superiortemporal ROI size in mm^3
  \item CTX_RH_SUPRAMARGINAL N av1451 uptake - subregion ctx-rh-supramarginal
  \item CTX_RH_SUPRAMARGINAL_SIZE N mm^3 ctx-rh-supramarginal ROI size in mm^3
  \item CTX_RH_TEMPORALPOLE N av1451 uptake - subregion ctx-rh-temporalpole
  \item CTX_RH_TEMPORALPOLE_SIZE N mm^3 ctx-rh-temporalpole ROI size in mm^3
  \item CTX_RH_TRANSVERSETEMPORAL N av1451 uptake - subregion ctx-rh-transversetemporal
  \item CTX_RH_TRANSVERSETEMPORAL_SIZE N mm^3 ctx-rh-transversetemporal ROI size in mm^3
  \item CTX_RH_UNKNOWN N av1451 uptake - subregion ctx-rh-unknown
  \item CTX_RH_UNKNOWN_SIZE N mm^3 ctx-rh-unknown ROI size in mm^3
  \item LEFT_ACCUMBENS_AREA N av1451 uptake - subregion left-accumbens-area
  \item LEFT_ACCUMBENS_AREA_SIZE N mm^3 left-accumbens-area ROI size in mm^3
  \item LEFT_AMYGDALA N av1451 uptake - subregion left-amygdala
  \item LEFT_AMYGDALA_SIZE N mm^3 left-amygdala ROI size in mm^3
  \item LEFT_CAUDATE N av1451 uptake - subregion left-caudate
  \item LEFT_CAUDATE_SIZE N mm^3 left-caudate ROI size in mm^3
  \item LEFT_CEREBELLUM_CORTEX N av1451 uptake - subregion left-cerebellum-cortex
  \item LEFT_CEREBELLUM_CORTEX_SIZE N mm^3 left-cerebellum-cortex ROI size in mm^3
  \item LEFT_CEREBELLUM_WHITE_MATTER N av1451 uptake - subregion left-cerebellum-white-matter
  \item LEFT_CEREBELLUM_WHITE_MATTER_SIZE N mm^3 left-cerebellum-white-matter ROI size in mm^3
  \item LEFT_CEREBRAL_WHITE_MATTER N av1451 uptake - subregion left-cerebral-white-matter
  \item LEFT_CEREBRAL_WHITE_MATTER_SIZE N mm^3 left-cerebral-white-matter ROI size in mm^3
  \item LEFT_CHOROID_PLEXUS N av1451 uptake - subregion left-choroid-plexus
  \item LEFT_CHOROID_PLEXUS_SIZE N mm^3 left-choroid-plexus ROI size in mm^3
  \item LEFT_HIPPOCAMPUS N av1451 uptake - subregion left-hippocampus
  \item LEFT_HIPPOCAMPUS_SIZE N mm^3 left-hippocampus ROI size in mm^3
  \item LEFT_INF_LAT_VENT N av1451 uptake - subregion left-inf-lat-vent
  \item LEFT_INF_LAT_VENT_SIZE N mm^3 left-inf-lat-vent ROI size in mm^3
  \item LEFT_LATERAL_VENTRICLE N av1451 uptake - subregion left-lateral-ventricle
  \item LEFT_LATERAL_VENTRICLE_SIZE N mm^3 left-lateral-ventricle ROI size in mm^3
  \item LEFT_PALLIDUM N av1451 uptake - subregion left-pallidum
  \item LEFT_PALLIDUM_SIZE N mm^3 left-pallidum ROI size in mm^3
  \item LEFT_PUTAMEN N av1451 uptake - subregion left-putamen
  \item LEFT_PUTAMEN_SIZE N mm^3 left-putamen ROI size in mm^3
  \item LEFT_THALAMUS_PROPER N av1451 uptake - subregion left-thalamus-proper
  \item LEFT_THALAMUS_PROPER_SIZE N mm^3 left-thalamus-proper ROI size in mm^3
  \item LEFT_VENTRALDC N av1451 uptake - subregion left-ventraldc
  \item LEFT_VENTRALDC_SIZE N mm^3 left-ventraldc ROI size in mm^3
  \item LEFT_VESSEL N av1451 uptake - subregion left-vessel
  \item LEFT_VESSEL_SIZE N mm^3 left-vessel ROI size in mm^3
  \item NON_WM_HYPOINTENSITIES N av1451 uptake - subregion non-wm-hypointensities
  \item NON_WM_HYPOINTENSITIES_SIZE N mm^3 non-wm-hypointensities ROI size in mm^3
  \item OPTIC_CHIASM N av1451 uptake - subregion optic-chiasm
  \item OPTIC_CHIASM_SIZE N mm^3 optic-chiasm ROI size in mm^3
  \item RIGHT_ACCUMBENS_AREA N av1451 uptake - subregion right-accumbens-area
  \item RIGHT_ACCUMBENS_AREA_SIZE N mm^3 right-accumbens-area ROI size in mm^3
  \item RIGHT_AMYGDALA N av1451 uptake - subregion right-amygdala
  \item RIGHT_AMYGDALA_SIZE N mm^3 right-amygdala ROI size in mm^3
  \item RIGHT_CAUDATE N av1451 uptake - subregion right-caudate
  \item RIGHT_CAUDATE_SIZE N mm^3 right-caudate ROI size in mm^3
  \item RIGHT_CEREBELLUM_CORTEX N av1451 uptake - subregion right-cerebellum-cortex
  \item RIGHT_CEREBELLUM_CORTEX_SIZE N mm^3 right-cerebellum-cortex ROI size in mm^3
  \item RIGHT_CEREBELLUM_WHITE_MATTER N av1451 uptake - subregion right-cerebellum-white-matter
  \item RIGHT_CEREBELLUM_WHITE_MATTER_SIZE N mm^3 right-cerebellum-white-matter ROI size in mm^3
  \item RIGHT_CEREBRAL_WHITE_MATTER N av1451 uptake - subregion right-cerebral-white-matter
  \item RIGHT_CEREBRAL_WHITE_MATTER_SIZE N mm^3 right-cerebral-white-matter ROI size in mm^3
  \item RIGHT_CHOROID_PLEXUS N av1451 uptake - subregion right-choroid-plexus
  \item RIGHT_CHOROID_PLEXUS_SIZE N mm^3 right-choroid-plexus ROI size in mm^3
  \item RIGHT_HIPPOCAMPUS N av1451 uptake - subregion right-hippocampus
  \item RIGHT_HIPPOCAMPUS_SIZE N mm^3 right-hippocampus ROI size in mm^3
  \item RIGHT_INF_LAT_VENT N av1451 uptake - subregion right-inf-lat-vent
  \item RIGHT_INF_LAT_VENT_SIZE N mm^3 right-inf-lat-vent ROI size in mm^3
  \item RIGHT_LATERAL_VENTRICLE N av1451 uptake - subregion right-lateral-ventricle
  \item RIGHT_LATERAL_VENTRICLE_SIZE N mm^3 right-lateral-ventricle ROI size in mm^3
  \item RIGHT_PALLIDUM N av1451 uptake - subregion right-pallidum
  \item RIGHT_PALLIDUM_SIZE N mm^3 right-pallidum ROI size in mm^3
  \item RIGHT_PUTAMEN N av1451 uptake - subregion right-putamen
  \item RIGHT_PUTAMEN_SIZE N mm^3 right-putamen ROI size in mm^3
  \item RIGHT_THALAMUS_PROPER N av1451 uptake - subregion right-thalamus-proper
  \item RIGHT_THALAMUS_PROPER_SIZE N mm^3 right-thalamus-proper ROI size in mm^3
  \item RIGHT_VENTRALDC N av1451 uptake - subregion right-ventraldc
  \item RIGHT_VENTRALDC_SIZE N mm^3 right-ventraldc ROI size in mm^3
  \item RIGHT_VESSEL N av1451 uptake - subregion right-vessel
  \item RIGHT_VESSEL_SIZE N mm^3 right-vessel ROI size in mm^3
  \item WM_HYPOINTENSITIES N av1451 uptake - subregion wm-hypointensities
  \item WM_HYPOINTENSITIES_SIZE N mm^3 wm-hypointensities ROI size in mm^3
}
}
\examples{
\dontrun{
describe(ucberkeleyav1451)
}
}
\seealso{
\url{https://adni.bitbucket.io/reference/docs/UCBERKELEYAV1451/UCBERKELEY_AV1451_Methods_03.02.16.pdf}, \url{https://adni.bitbucket.io/reference/docs/UCBERKELEYAV1451/UCBERKELEY_AV1451_Methods_04.25.2022.pdf}, \url{https://adni.bitbucket.io/reference/docs/UCBERKELEYAV1451/UCBERKELEY_AV1451_Methods_11.15.2021.pdf}, \url{https://adni.bitbucket.io/reference/docs/UCBERKELEYAV1451/UCBERKELEYAV1451_Methods_20171114.pdf}, \url{https://adni.bitbucket.io/reference/docs/UCBERKELEYAV1451/UCBERKELEY_AV1451_Methods_2021-01-14.pdf}, \url{https://adni.bitbucket.io/reference/docs/UCBERKELEYAV1451/UCBERKELEY_AV1451_Methods_Aug2018.pdf}, \url{https://adni.bitbucket.io/reference/docs/UCBERKELEYAV1451/UCBERKELEY_AV1451_Methods_Dec2019.pdf}
}
\keyword{datasets}
