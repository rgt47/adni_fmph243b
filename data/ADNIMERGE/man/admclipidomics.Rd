% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataDoc.R
\docType{data}
\name{admclipidomics}
\alias{admclipidomics}
\title{ADMC Lipidomics Measured Lipid Molecules}
\format{
A data frame with 910 rows and 525 variables
}
\usage{
data(admclipidomics)
}
\description{
ADMC Lipidomics Measured Lipid Molecules
}
\details{
A dataset with variables as follows:

\itemize{
  \item RID -4 -4 Participant roster ID
  \item RECNO -4 -4 Record number (2 indicates replicate sample)
  \item VARID -4 -4 -4
  \item UCD.Lipid.1 -4 -4 AC (8:0)
  \item UCD.Lipid.2 -4 -4 AC (8:1)
  \item UCD.Lipid.3 -4 -4 AC (10:1)
  \item UCD.Lipid.4 -4 -4 AC (12:0)
  \item UCD.Lipid.5 -4 -4 AC (16:0)
  \item UCD.Lipid.6 -4 -4 AC (18:0)
  \item UCD.Lipid.7 -4 -4 AC (18:1)
  \item UCD.Lipid.8 -4 -4 AC (18:2)
  \item UCD.Lipid.9 -4 -4 AC (18:3)
  \item UCD.Lipid.10 -4 -4 CE (16:1)
  \item UCD.Lipid.11 -4 -4 CE (18:1)
  \item UCD.Lipid.12 -4 -4 CE (18:2)
  \item UCD.Lipid.13 -4 -4 CE (18:3)
  \item UCD.Lipid.14 -4 -4 CE (20:3)
  \item UCD.Lipid.15 -4 -4 CE (20:4)
  \item UCD.Lipid.16 -4 -4 CE (20:5)
  \item UCD.Lipid.17 -4 -4 CE (22:6)
  \item UCD.Lipid.18 -4 -4 Cer (41:1)
  \item UCD.Lipid.19 -4 -4 Cer (42:2)
  \item UCD.Lipid.20 -4 -4 Cholesterol
  \item UCD.Lipid.21 -4 -4 CSH_posESI #002
  \item UCD.Lipid.22 -4 -4 CSH_posESI #003
  \item UCD.Lipid.23 -4 -4 CSH_posESI #004
  \item UCD.Lipid.24 -4 -4 CSH_posESI #006
  \item UCD.Lipid.25 -4 -4 CSH_posESI #015
  \item UCD.Lipid.26 -4 -4 CSH_posESI #017
  \item UCD.Lipid.27 -4 -4 CSH_posESI #022
  \item UCD.Lipid.28 -4 -4 CSH_posESI #025
  \item UCD.Lipid.29 -4 -4 CSH_posESI #027
  \item UCD.Lipid.30 -4 -4 CSH_posESI #033
  \item UCD.Lipid.31 -4 -4 CSH_posESI #034
  \item UCD.Lipid.32 -4 -4 CSH_posESI #036
  \item UCD.Lipid.33 -4 -4 CSH_posESI #037
  \item UCD.Lipid.34 -4 -4 CSH_posESI #040
  \item UCD.Lipid.35 -4 -4 CSH_posESI #041
  \item UCD.Lipid.36 -4 -4 CSH_posESI #042
  \item UCD.Lipid.37 -4 -4 CSH_posESI #043
  \item UCD.Lipid.38 -4 -4 CSH_posESI #056
  \item UCD.Lipid.39 -4 -4 CSH_posESI #057
  \item UCD.Lipid.40 -4 -4 CSH_posESI #058
  \item UCD.Lipid.41 -4 -4 CSH_posESI #061
  \item UCD.Lipid.42 -4 -4 CSH_posESI #079
  \item UCD.Lipid.43 -4 -4 CSH_posESI #081
  \item UCD.Lipid.44 -4 -4 CSH_posESI #091
  \item UCD.Lipid.45 -4 -4 CSH_posESI #098
  \item UCD.Lipid.46 -4 -4 CSH_posESI #103
  \item UCD.Lipid.47 -4 -4 CSH_posESI #104
  \item UCD.Lipid.48 -4 -4 CSH_posESI #106
  \item UCD.Lipid.49 -4 -4 CSH_posESI #109
  \item UCD.Lipid.50 -4 -4 CSH_posESI #116
  \item UCD.Lipid.51 -4 -4 CSH_posESI #117
  \item UCD.Lipid.52 -4 -4 CSH_posESI #121
  \item UCD.Lipid.53 -4 -4 CSH_posESI #122
  \item UCD.Lipid.54 -4 -4 CSH_posESI #124
  \item UCD.Lipid.55 -4 -4 CSH_posESI #125
  \item UCD.Lipid.56 -4 -4 CSH_posESI #127
  \item UCD.Lipid.57 -4 -4 CSH_posESI #130
  \item UCD.Lipid.58 -4 -4 CSH_posESI #134
  \item UCD.Lipid.59 -4 -4 CSH_posESI #137
  \item UCD.Lipid.60 -4 -4 CSH_posESI #141
  \item UCD.Lipid.61 -4 -4 CSH_posESI #142
  \item UCD.Lipid.62 -4 -4 CSH_posESI #155
  \item UCD.Lipid.63 -4 -4 CSH_posESI #156
  \item UCD.Lipid.64 -4 -4 CSH_posESI #159
  \item UCD.Lipid.65 -4 -4 CSH_posESI #164
  \item UCD.Lipid.66 -4 -4 CSH_posESI #165
  \item UCD.Lipid.67 -4 -4 CSH_posESI #168
  \item UCD.Lipid.68 -4 -4 CSH_posESI #170
  \item UCD.Lipid.69 -4 -4 CSH_posESI #176
  \item UCD.Lipid.70 -4 -4 CSH_posESI #187
  \item UCD.Lipid.71 -4 -4 CSH_posESI #192
  \item UCD.Lipid.72 -4 -4 CSH_posESI #196
  \item UCD.Lipid.73 -4 -4 CSH_posESI #197
  \item UCD.Lipid.74 -4 -4 CSH_posESI #201
  \item UCD.Lipid.75 -4 -4 CSH_posESI #203
  \item UCD.Lipid.76 -4 -4 CSH_posESI #205
  \item UCD.Lipid.77 -4 -4 CSH_posESI #211
  \item UCD.Lipid.78 -4 -4 CSH_posESI #214
  \item UCD.Lipid.79 -4 -4 CSH_posESI #215
  \item UCD.Lipid.80 -4 -4 CSH_posESI #230
  \item UCD.Lipid.81 -4 -4 CSH_posESI #231
  \item UCD.Lipid.82 -4 -4 CSH_posESI #232
  \item UCD.Lipid.83 -4 -4 CSH_posESI #270
  \item UCD.Lipid.84 -4 -4 CSH_posESI #275
  \item UCD.Lipid.85 -4 -4 CSH_posESI #282
  \item UCD.Lipid.86 -4 -4 CSH_posESI #285
  \item UCD.Lipid.87 -4 -4 CSH_posESI #286
  \item UCD.Lipid.88 -4 -4 CSH_posESI #288
  \item UCD.Lipid.89 -4 -4 CSH_posESI #290
  \item UCD.Lipid.90 -4 -4 CSH_posESI #293
  \item UCD.Lipid.91 -4 -4 CSH_posESI #296
  \item UCD.Lipid.92 -4 -4 CSH_posESI #298
  \item UCD.Lipid.93 -4 -4 CSH_posESI #301
  \item UCD.Lipid.94 -4 -4 CSH_posESI #303
  \item UCD.Lipid.95 -4 -4 CSH_posESI #304
  \item UCD.Lipid.96 -4 -4 CSH_posESI #307
  \item UCD.Lipid.97 -4 -4 CSH_posESI #308
  \item UCD.Lipid.98 -4 -4 CSH_posESI #310
  \item UCD.Lipid.99 -4 -4 CSH_posESI #336
  \item UCD.Lipid.100 -4 -4 CSH_posESI #359
  \item UCD.Lipid.101 -4 -4 CSH_posESI #362
  \item UCD.Lipid.102 -4 -4 CSH_posESI #363
  \item UCD.Lipid.103 -4 -4 CSH_posESI #364
  \item UCD.Lipid.104 -4 -4 DG (32:1)
  \item UCD.Lipid.105 -4 -4 DG (34:1)
  \item UCD.Lipid.106 -4 -4 DG (34:2)
  \item UCD.Lipid.107 -4 -4 DG (34:3)
  \item UCD.Lipid.108 -4 -4 DG (36:1)
  \item UCD.Lipid.109 -4 -4 DG (36:2)
  \item UCD.Lipid.110 -4 -4 DG (36:3)
  \item UCD.Lipid.111 -4 -4 DG (36:4)
  \item UCD.Lipid.112 -4 -4 DG (36:5)
  \item UCD.Lipid.113 -4 -4 DG (38:0)
  \item UCD.Lipid.114 -4 -4 DG (38:3)
  \item UCD.Lipid.115 -4 -4 DG (38:5)
  \item UCD.Lipid.116 -4 -4 DG (38:6)
  \item UCD.Lipid.117 -4 -4 Hex2Cer (34:1)
  \item UCD.Lipid.118 -4 -4 LacCer (42:2)
  \item UCD.Lipid.119 -4 -4 LPC (14:0)
  \item UCD.Lipid.120 -4 -4 LPC (15:0)
  \item UCD.Lipid.121 -4 -4 LPC (17:1)
  \item UCD.Lipid.122 -4 -4 LPC (18:0)
  \item UCD.Lipid.123 -4 -4 LPC (18:3)
  \item UCD.Lipid.124 -4 -4 LPC (20:0)
  \item UCD.Lipid.125 -4 -4 LPC (20:4)
  \item UCD.Lipid.126 -4 -4 LPC (20:5)
  \item UCD.Lipid.127 -4 -4 LPC (22:4)
  \item UCD.Lipid.128 -4 -4 LPC (22:6)
  \item UCD.Lipid.129 -4 -4 LPC (o-16:0)
  \item UCD.Lipid.130 -4 -4 LPC (p-16:0) or LPC (o-16:1)
  \item UCD.Lipid.131 -4 -4 LPC (p-18:0) or LPC (o-18:1)
  \item UCD.Lipid.132 -4 -4 PC (25:0(CHO))
  \item UCD.Lipid.133 -4 -4 PC (28:0)
  \item UCD.Lipid.134 -4 -4 PC (30:0)
  \item UCD.Lipid.135 -4 -4 PC (30:1)
  \item UCD.Lipid.136 -4 -4 PC (31:0)
  \item UCD.Lipid.137 -4 -4 PC (31:1)
  \item UCD.Lipid.138 -4 -4 PC (32:3)
  \item UCD.Lipid.139 -4 -4 PC (33:0)
  \item UCD.Lipid.140 -4 -4 PC (35:3)
  \item UCD.Lipid.141 -4 -4 PC (36:3)
  \item UCD.Lipid.142 -4 -4 PC (36:4)
  \item UCD.Lipid.143 -4 -4 PC (36:6)
  \item UCD.Lipid.144 -4 -4 PC (37:3)
  \item UCD.Lipid.145 -4 -4 PC (37:6)
  \item UCD.Lipid.146 -4 -4 PC (38:4)
  \item UCD.Lipid.147 -4 -4 PC (38:5)
  \item UCD.Lipid.148 -4 -4 PC (38:6)
  \item UCD.Lipid.149 -4 -4 PC (38:7)
  \item UCD.Lipid.150 -4 -4 PC (39:6)
  \item UCD.Lipid.151 -4 -4 PC (40:5)
  \item UCD.Lipid.152 -4 -4 PC (40:6)
  \item UCD.Lipid.153 -4 -4 PC (42:5)
  \item UCD.Lipid.154 -4 -4 PC (42:6)
  \item UCD.Lipid.155 -4 -4 PC (o-34:0)
  \item UCD.Lipid.156 -4 -4 PC (p-32:1) or PC (o-32:2)
  \item UCD.Lipid.157 -4 -4 PC (p-34:1) or PC (o-34:2)
  \item UCD.Lipid.158 -4 -4 PC (p-38:2) or PC (o-38:3)
  \item UCD.Lipid.159 -4 -4 PC (p-38:4) or PC (o-38:5)
  \item UCD.Lipid.160 -4 -4 PC (p-40:1) or PC (o-40:2)
  \item UCD.Lipid.161 -4 -4 PC (p-40:5) or PC (o-40:6)
  \item UCD.Lipid.162 -4 -4 PC (p-40:6) or PC (o-40:7) A
  \item UCD.Lipid.163 -4 -4 PC (p-40:6) or PC (o-40:7) B
  \item UCD.Lipid.164 -4 -4 PC (p-40:7) or PC (o-40:8)
  \item UCD.Lipid.165 -4 -4 PC (p-42:3) or PC (o-42:4)
  \item UCD.Lipid.166 -4 -4 PE (36:4)
  \item UCD.Lipid.167 -4 -4 PE (38:4)
  \item UCD.Lipid.168 -4 -4 PE (p-34:1) or PE (o-34:2)
  \item UCD.Lipid.169 -4 -4 SM (30:1)
  \item UCD.Lipid.170 -4 -4 SM (34:0)
  \item UCD.Lipid.171 -4 -4 SM (40:2)
  \item UCD.Lipid.172 -4 -4 SM (42:2)
  \item UCD.Lipid.173 -4 -4 TG (42:0)
  \item UCD.Lipid.174 -4 -4 TG (40:0)
  \item UCD.Lipid.175 -4 -4 TG (40:1)
  \item UCD.Lipid.176 -4 -4 TG (42:1)
  \item UCD.Lipid.177 -4 -4 TG (42:2)
  \item UCD.Lipid.178 -4 -4 TG (42:3)
  \item UCD.Lipid.179 -4 -4 TG (44:0)
  \item UCD.Lipid.180 -4 -4 TG (44:1)
  \item UCD.Lipid.181 -4 -4 TG (44:2)
  \item UCD.Lipid.182 -4 -4 TG (46:0)
  \item UCD.Lipid.183 -4 -4 TG (46:1)
  \item UCD.Lipid.184 -4 -4 TG (46:2)
  \item UCD.Lipid.185 -4 -4 TG (46:3)
  \item UCD.Lipid.186 -4 -4 TG (46:4)
  \item UCD.Lipid.187 -4 -4 TG (48:0)
  \item UCD.Lipid.188 -4 -4 TG (48:1)
  \item UCD.Lipid.189 -4 -4 TG (48:2)
  \item UCD.Lipid.190 -4 -4 TG (48:3)
  \item UCD.Lipid.191 -4 -4 TG (48:4)
  \item UCD.Lipid.192 -4 -4 TG (49:0)
  \item UCD.Lipid.193 -4 -4 TG (49:1)
  \item UCD.Lipid.194 -4 -4 TG (49:2)
  \item UCD.Lipid.195 -4 -4 TG (49:3)
  \item UCD.Lipid.196 -4 -4 TG (50:0)
  \item UCD.Lipid.197 -4 -4 TG (50:1)
  \item UCD.Lipid.198 -4 -4 TG (50:2)
  \item UCD.Lipid.199 -4 -4 TG (50:3)
  \item UCD.Lipid.200 -4 -4 TG (50:4)
  \item UCD.Lipid.201 -4 -4 TG (50:5)
  \item UCD.Lipid.202 -4 -4 TG (50:6)
  \item UCD.Lipid.203 -4 -4 TG (51:1)
  \item UCD.Lipid.204 -4 -4 TG (51:2)
  \item UCD.Lipid.205 -4 -4 TG (51:3)
  \item UCD.Lipid.206 -4 -4 TG (51:4)
  \item UCD.Lipid.207 -4 -4 TG (51:5)
  \item UCD.Lipid.208 -4 -4 TG (52:0)
  \item UCD.Lipid.209 -4 -4 TG (52:1)
  \item UCD.Lipid.210 -4 -4 TG (52:2)
  \item UCD.Lipid.211 -4 -4 TG (52:3)
  \item UCD.Lipid.212 -4 -4 TG (52:4)
  \item UCD.Lipid.213 -4 -4 TG (52:5)
  \item UCD.Lipid.214 -4 -4 TG (52:6)
  \item UCD.Lipid.215 -4 -4 TG (53:1)
  \item UCD.Lipid.216 -4 -4 TG (53:2)
  \item UCD.Lipid.217 -4 -4 TG (53:3)
  \item UCD.Lipid.218 -4 -4 TG (53:4)
  \item UCD.Lipid.219 -4 -4 TG (53:5)
  \item UCD.Lipid.220 -4 -4 TG (54:0)
  \item UCD.Lipid.221 -4 -4 TG (54:1)
  \item UCD.Lipid.222 -4 -4 TG (54:2)
  \item UCD.Lipid.223 -4 -4 TG (54:3)
  \item UCD.Lipid.224 -4 -4 TG (54:4)
  \item UCD.Lipid.225 -4 -4 TG (54:5)
  \item UCD.Lipid.226 -4 -4 TG (54:6)
  \item UCD.Lipid.227 -4 -4 TG (54:8)
  \item UCD.Lipid.228 -4 -4 TG (56:1)
  \item UCD.Lipid.229 -4 -4 TG (56:2)
  \item UCD.Lipid.230 -4 -4 TG (56:3)
  \item UCD.Lipid.231 -4 -4 TG (56:4)
  \item UCD.Lipid.232 -4 -4 TG (56:5)
  \item UCD.Lipid.233 -4 -4 TG (56:6)
  \item UCD.Lipid.234 -4 -4 TG (56:7)
  \item UCD.Lipid.235 -4 -4 TG (56:8)
  \item UCD.Lipid.236 -4 -4 TG (56:9)
  \item UCD.Lipid.237 -4 -4 TG (57:1)
  \item UCD.Lipid.238 -4 -4 TG (57:2)
  \item UCD.Lipid.239 -4 -4 TG (58:1)
  \item UCD.Lipid.240 -4 -4 TG (58:10)
  \item UCD.Lipid.241 -4 -4 TG (58:2)
  \item UCD.Lipid.242 -4 -4 TG (58:3)
  \item UCD.Lipid.243 -4 -4 TG (58:4)
  \item UCD.Lipid.244 -4 -4 TG (58:6)
  \item UCD.Lipid.245 -4 -4 TG (58:8)
  \item UCD.Lipid.246 -4 -4 TG (58:9)
  \item UCD.Lipid.247 -4 -4 TG (59:2)
  \item UCD.Lipid.248 -4 -4 TG (59:3)
  \item UCD.Lipid.249 -4 -4 TG (60:11)
  \item UCD.Lipid.250 -4 -4 TG (60:2)
  \item UCD.Lipid.251 -4 -4 TG (60:3)
  \item UCD.Lipid.252 -4 -4 TG (60:4)
  \item UCD.Lipid.253 -4 -4 TG (60:6)
  \item UCD.Lipid.254 -4 -4 TG (62:3)
  \item UCD.Lipid.255 -4 -4 TG (62:4)
  \item UCD.Lipid.256 -4 -4 TG (64:4)
  \item UCD.Lipid.257 -4 -4 Cer (32:1)
  \item UCD.Lipid.258 -4 -4 Cer (33:1)
  \item UCD.Lipid.259 -4 -4 Cer (34:0)
  \item UCD.Lipid.260 -4 -4 Cer (34:1)
  \item UCD.Lipid.261 -4 -4 Cer (34:2)
  \item UCD.Lipid.262 -4 -4 Cer (36:1)
  \item UCD.Lipid.263 -4 -4 Cer (38:1)
  \item UCD.Lipid.264 -4 -4 Cer (39:1)
  \item UCD.Lipid.265 -4 -4 Cer (40:0)
  \item UCD.Lipid.266 -4 -4 Cer (40:1)
  \item UCD.Lipid.267 -4 -4 Cer (40:2)
  \item UCD.Lipid.268 -4 -4 Cer (41:1)
  \item UCD.Lipid.269 -4 -4 Cer (42:0)
  \item UCD.Lipid.270 -4 -4 Cer (42:1)
  \item UCD.Lipid.271 -4 -4 Cer (42:2)
  \item UCD.Lipid.272 -4 -4 Cer (43:1)
  \item UCD.Lipid.273 -4 -4 Cer (44:1)
  \item UCD.Lipid.274 -4 -4 CSH_negESI #014
  \item UCD.Lipid.275 -4 -4 CSH_negESI #056
  \item UCD.Lipid.276 -4 -4 CSH_negESI #109
  \item UCD.Lipid.277 -4 -4 CSH_negESI #113
  \item UCD.Lipid.278 -4 -4 CSH_negESI #130
  \item UCD.Lipid.279 -4 -4 CSH_negESI #131
  \item UCD.Lipid.280 -4 -4 CSH_negESI #223
  \item UCD.Lipid.281 -4 -4 CSH_negESI #293
  \item UCD.Lipid.282 -4 -4 CSH_negESI #377
  \item UCD.Lipid.283 -4 -4 CSH_negESI #378
  \item UCD.Lipid.284 -4 -4 CSH_negESI #387
  \item UCD.Lipid.285 -4 -4 CSH_negESI #405
  \item UCD.Lipid.286 -4 -4 CSH_negESI #413
  \item UCD.Lipid.287 -4 -4 CSH_negESI #414
  \item UCD.Lipid.288 -4 -4 CSH_negESI #420
  \item UCD.Lipid.289 -4 -4 CSH_negESI #501
  \item UCD.Lipid.290 -4 -4 CSH_negESI #502
  \item UCD.Lipid.291 -4 -4 CSH_negESI #504
  \item UCD.Lipid.292 -4 -4 CSH_negESI #510
  \item UCD.Lipid.293 -4 -4 CSH_negESI #511
  \item UCD.Lipid.294 -4 -4 CSH_negESI #520
  \item UCD.Lipid.295 -4 -4 CSH_negESI #524
  \item UCD.Lipid.296 -4 -4 CSH_negESI #527
  \item UCD.Lipid.297 -4 -4 CSH_negESI #528
  \item UCD.Lipid.298 -4 -4 CSH_negESI #530
  \item UCD.Lipid.299 -4 -4 CSH_negESI #550
  \item UCD.Lipid.300 -4 -4 CSH_negESI #564
  \item UCD.Lipid.301 -4 -4 CSH_negESI #565
  \item UCD.Lipid.302 -4 -4 CSH_negESI #569
  \item UCD.Lipid.303 -4 -4 CSH_negESI #573
  \item UCD.Lipid.304 -4 -4 CSH_negESI #579
  \item UCD.Lipid.305 -4 -4 CSH_negESI #580
  \item UCD.Lipid.306 -4 -4 CSH_negESI #581
  \item UCD.Lipid.307 -4 -4 CSH_negESI #586
  \item UCD.Lipid.308 -4 -4 CSH_negESI #587
  \item UCD.Lipid.309 -4 -4 CSH_negESI #589
  \item UCD.Lipid.310 -4 -4 CSH_negESI #590
  \item UCD.Lipid.311 -4 -4 CSH_negESI #591
  \item UCD.Lipid.312 -4 -4 CSH_negESI #592
  \item UCD.Lipid.313 -4 -4 CSH_negESI #593
  \item UCD.Lipid.314 -4 -4 CSH_negESI #598
  \item UCD.Lipid.315 -4 -4 CSH_negESI #599
  \item UCD.Lipid.316 -4 -4 CSH_negESI #600
  \item UCD.Lipid.317 -4 -4 CSH_negESI #612
  \item UCD.Lipid.318 -4 -4 CSH_negESI #615
  \item UCD.Lipid.319 -4 -4 CSH_negESI #617
  \item UCD.Lipid.320 -4 -4 CSH_negESI #624
  \item UCD.Lipid.321 -4 -4 CSH_negESI #627
  \item UCD.Lipid.322 -4 -4 CSH_negESI #628
  \item UCD.Lipid.323 -4 -4 CSH_negESI #631
  \item UCD.Lipid.324 -4 -4 CSH_negESI #633
  \item UCD.Lipid.325 -4 -4 CSH_negESI #638
  \item UCD.Lipid.326 -4 -4 CSH_negESI #639
  \item UCD.Lipid.327 -4 -4 CSH_negESI #644
  \item UCD.Lipid.328 -4 -4 CSH_negESI #648
  \item UCD.Lipid.329 -4 -4 CSH_negESI #650
  \item UCD.Lipid.330 -4 -4 CSH_negESI #653
  \item UCD.Lipid.331 -4 -4 CSH_negESI #658
  \item UCD.Lipid.332 -4 -4 CSH_negESI #668
  \item UCD.Lipid.333 -4 -4 CSH_negESI #670
  \item UCD.Lipid.334 -4 -4 CSH_negESI #678
  \item UCD.Lipid.335 -4 -4 CSH_negESI #682
  \item UCD.Lipid.336 -4 -4 CSH_negESI #683
  \item UCD.Lipid.337 -4 -4 CSH_negESI #690
  \item UCD.Lipid.338 -4 -4 CSH_negESI #691
  \item UCD.Lipid.339 -4 -4 CSH_negESI #693
  \item UCD.Lipid.340 -4 -4 CSH_negESI #695
  \item UCD.Lipid.341 -4 -4 CSH_negESI #698
  \item UCD.Lipid.342 -4 -4 CSH_negESI #701
  \item UCD.Lipid.343 -4 -4 CSH_negESI #706
  \item UCD.Lipid.344 -4 -4 CSH_negESI #709
  \item UCD.Lipid.345 -4 -4 CSH_negESI #710
  \item UCD.Lipid.346 -4 -4 CSH_negESI #713
  \item UCD.Lipid.347 -4 -4 CSH_negESI #716
  \item UCD.Lipid.348 -4 -4 CSH_negESI #722
  \item UCD.Lipid.349 -4 -4 CSH_negESI #723
  \item UCD.Lipid.350 -4 -4 CSH_negESI #724
  \item UCD.Lipid.351 -4 -4 CSH_negESI #728
  \item UCD.Lipid.352 -4 -4 CSH_negESI #733
  \item UCD.Lipid.353 -4 -4 CSH_negESI #735
  \item UCD.Lipid.354 -4 -4 CSH_negESI #743
  \item UCD.Lipid.355 -4 -4 CSH_negESI #744
  \item UCD.Lipid.356 -4 -4 CSH_negESI #749
  \item UCD.Lipid.357 -4 -4 CSH_negESI #750
  \item UCD.Lipid.358 -4 -4 CSH_negESI #762
  \item UCD.Lipid.359 -4 -4 CSH_negESI #764
  \item UCD.Lipid.360 -4 -4 CSH_negESI #767
  \item UCD.Lipid.361 -4 -4 CSH_negESI #770
  \item UCD.Lipid.362 -4 -4 CSH_negESI #773
  \item UCD.Lipid.363 -4 -4 HexCer (38:1)
  \item UCD.Lipid.364 -4 -4 HexCer (40:1)
  \item UCD.Lipid.365 -4 -4 HexCer (41:1)
  \item UCD.Lipid.366 -4 -4 HexCer (42:1)
  \item UCD.Lipid.367 -4 -4 HexCer (42:2)
  \item UCD.Lipid.368 -4 -4 HexCer (34:1(2OH))
  \item UCD.Lipid.369 -4 -4 LPC (16:0)
  \item UCD.Lipid.370 -4 -4 LPC (16:1)
  \item UCD.Lipid.371 -4 -4 LPC (18:0)
  \item UCD.Lipid.372 -4 -4 LPC (18:1)
  \item UCD.Lipid.373 -4 -4 LPC (18:2)
  \item UCD.Lipid.374 -4 -4 LPC (20:1)
  \item UCD.Lipid.375 -4 -4 LPC (20:2)
  \item UCD.Lipid.376 -4 -4 LPC (20:3)
  \item UCD.Lipid.377 -4 -4 LPC (22:5)
  \item UCD.Lipid.378 -4 -4 LPE (16:0)
  \item UCD.Lipid.379 -4 -4 LPE (18:2)
  \item UCD.Lipid.380 -4 -4 LPE (20:4)
  \item UCD.Lipid.381 -4 -4 LPE (22:6)
  \item UCD.Lipid.382 -4 -4 PC (32:0)
  \item UCD.Lipid.383 -4 -4 PC (32:1)
  \item UCD.Lipid.384 -4 -4 PC (32:2)
  \item UCD.Lipid.385 -4 -4 PC (33:1)
  \item UCD.Lipid.386 -4 -4 PC (33:2)
  \item UCD.Lipid.387 -4 -4 PC (34:0)
  \item UCD.Lipid.388 -4 -4 PC (34:1)
  \item UCD.Lipid.389 -4 -4 PC (34:2)
  \item UCD.Lipid.390 -4 -4 PC (34:3)
  \item UCD.Lipid.391 -4 -4 PC (34:4)
  \item UCD.Lipid.392 -4 -4 PC (35:1)
  \item UCD.Lipid.393 -4 -4 PC (35:2)
  \item UCD.Lipid.394 -4 -4 PC (35:4)
  \item UCD.Lipid.395 -4 -4 PC (36:1)
  \item UCD.Lipid.396 -4 -4 PC (36:2)
  \item UCD.Lipid.397 -4 -4 PC (36:3)
  \item UCD.Lipid.398 -4 -4 PC (36:4)
  \item UCD.Lipid.399 -4 -4 PC (36:5)
  \item UCD.Lipid.400 -4 -4 PC (37:2)
  \item UCD.Lipid.401 -4 -4 PC (37:4)
  \item UCD.Lipid.402 -4 -4 PC (38:2)
  \item UCD.Lipid.403 -4 -4 PC (38:3)
  \item UCD.Lipid.404 -4 -4 PC (38:4)
  \item UCD.Lipid.405 -4 -4 PC (38:5)
  \item UCD.Lipid.406 -4 -4 PC (38:6)
  \item UCD.Lipid.407 -4 -4 PC (40:4)
  \item UCD.Lipid.408 -4 -4 PC (40:5)
  \item UCD.Lipid.409 -4 -4 PC (40:6)
  \item UCD.Lipid.410 -4 -4 PC (40:7)
  \item UCD.Lipid.411 -4 -4 PC (40:8)
  \item UCD.Lipid.412 -4 -4 PC (o-32:0)
  \item UCD.Lipid.413 -4 -4 PC (p-32:0) or PC (o-32:1)
  \item UCD.Lipid.414 -4 -4 PC (p-34:0) or PC (o-34:1)
  \item UCD.Lipid.415 -4 -4 PC (p-34:1) or PC (o-34:2)
  \item UCD.Lipid.416 -4 -4 PC (p-34:2) or PC (o-34:3)
  \item UCD.Lipid.417 -4 -4 PC (p-36:1) or PC (o-36:2)
  \item UCD.Lipid.418 -4 -4 PC (p-36:2) or PC (o-36:3)
  \item UCD.Lipid.419 -4 -4 PC (p-36:3) or PC (o-36:4)
  \item UCD.Lipid.420 -4 -4 PC (p-36:4) or PC (o-36:5)
  \item UCD.Lipid.421 -4 -4 PC (p-38:3) or PC (o-38:4)
  \item UCD.Lipid.422 -4 -4 PC (p-38:4) or PC (o-38:5)
  \item UCD.Lipid.423 -4 -4 PC (p-38:5) or PC (o-38:6)
  \item UCD.Lipid.424 -4 -4 PC (p-40:3) or PC (o-40:4)
  \item UCD.Lipid.425 -4 -4 PC (p-40:4) or PC (o-40:5)
  \item UCD.Lipid.426 -4 -4 PC (p-42:4) or PC (o-42:5)
  \item UCD.Lipid.427 -4 -4 PC (p-42:5) or PC (o-42:6)
  \item UCD.Lipid.428 -4 -4 PC (p-44:4) or PC (o-44:5)
  \item UCD.Lipid.429 -4 -4 PE (34:1)
  \item UCD.Lipid.430 -4 -4 PE (34:2)
  \item UCD.Lipid.431 -4 -4 PE (36:1)
  \item UCD.Lipid.432 -4 -4 PE (36:2)
  \item UCD.Lipid.433 -4 -4 PE (36:3)
  \item UCD.Lipid.434 -4 -4 PE (38:2)
  \item UCD.Lipid.435 -4 -4 PE (38:4)
  \item UCD.Lipid.436 -4 -4 PE (38:6)
  \item UCD.Lipid.437 -4 -4 PE (40:6)
  \item UCD.Lipid.438 -4 -4 PE (p-34:2) or PE (o-34:3)
  \item UCD.Lipid.439 -4 -4 PE (p-36:1) or PE (o-36:2)
  \item UCD.Lipid.440 -4 -4 PE (p-36:2) or PE (o-36:3)
  \item UCD.Lipid.441 -4 -4 PE (p-36:4) or PE (o-36:5)
  \item UCD.Lipid.442 -4 -4 PE (p-36:5) or PE (o-36:6)
  \item UCD.Lipid.443 -4 -4 PE (p-38:2) or PE (o-38:3)
  \item UCD.Lipid.444 -4 -4 PE (p-38:3) or PE (o-38:4)
  \item UCD.Lipid.445 -4 -4 PE (p-38:4) or PE (o-38:5)
  \item UCD.Lipid.446 -4 -4 PE (p-38:5) or PE (o-38:6)
  \item UCD.Lipid.447 -4 -4 PE (p-38:6) or PE (o-38:7)
  \item UCD.Lipid.448 -4 -4 PE (p-40:4) or PE (o-40:5)
  \item UCD.Lipid.449 -4 -4 PE (p-40:5) or PE (o-40:6)
  \item UCD.Lipid.450 -4 -4 PE (p-40:6) or PE (o-40:7)
  \item UCD.Lipid.451 -4 -4 PE (p-40:7) or PE (o-40:8)
  \item UCD.Lipid.452 -4 -4 PI (32:1)
  \item UCD.Lipid.453 -4 -4 PI (34:1)
  \item UCD.Lipid.454 -4 -4 PI (34:2)
  \item UCD.Lipid.455 -4 -4 PI (36:1)
  \item UCD.Lipid.456 -4 -4 PI (36:2)
  \item UCD.Lipid.457 -4 -4 PI (36:3)
  \item UCD.Lipid.458 -4 -4 PI (36:4)
  \item UCD.Lipid.459 -4 -4 PI (38:3)
  \item UCD.Lipid.460 -4 -4 PI (38:4)
  \item UCD.Lipid.461 -4 -4 PI (38:5)
  \item UCD.Lipid.462 -4 -4 PI (40:6)
  \item UCD.Lipid.463 -4 -4 SM (32:0)
  \item UCD.Lipid.464 -4 -4 SM (32:1)
  \item UCD.Lipid.465 -4 -4 SM (32:2)
  \item UCD.Lipid.466 -4 -4 SM (33:1)
  \item UCD.Lipid.467 -4 -4 SM (34:1)
  \item UCD.Lipid.468 -4 -4 SM (34:2)
  \item UCD.Lipid.469 -4 -4 SM (36:0)
  \item UCD.Lipid.470 -4 -4 SM (36:1)
  \item UCD.Lipid.471 -4 -4 SM (36:2)
  \item UCD.Lipid.472 -4 -4 SM (36:3)
  \item UCD.Lipid.473 -4 -4 SM (37:1)
  \item UCD.Lipid.474 -4 -4 SM (38:0)
  \item UCD.Lipid.475 -4 -4 SM (38:1)
  \item UCD.Lipid.476 -4 -4 SM (38:2)
  \item UCD.Lipid.477 -4 -4 SM (39:1)
  \item UCD.Lipid.478 -4 -4 SM (39:2)
  \item UCD.Lipid.479 -4 -4 SM (40:0)
  \item UCD.Lipid.480 -4 -4 SM (40:1)
  \item UCD.Lipid.481 -4 -4 SM (40:2) A
  \item UCD.Lipid.482 -4 -4 SM (40:2) B
  \item UCD.Lipid.483 -4 -4 SM (40:3)
  \item UCD.Lipid.484 -4 -4 SM (41:1)
  \item UCD.Lipid.485 -4 -4 SM (41:2)
  \item UCD.Lipid.486 -4 -4 SM (42:0)
  \item UCD.Lipid.487 -4 -4 SM (42:1)
  \item UCD.Lipid.488 -4 -4 SM (42:2)
  \item UCD.Lipid.489 -4 -4 SM (42:3)
  \item UCD.Lipid.490 -4 -4 SM (43:1)
  \item UCD.Lipid.491 -4 -4 SM (43:2)
  \item UCD.Lipid.492 -4 -4 SM (44:2)
  \item UCD.Lipid.493 -4 -4 FA (11:0)
  \item UCD.Lipid.494 -4 -4 FA (12:0)
  \item UCD.Lipid.495 -4 -4 FA (13:0)
  \item UCD.Lipid.496 -4 -4 FA (14:0)
  \item UCD.Lipid.497 -4 -4 FA (14:1)
  \item UCD.Lipid.498 -4 -4 FA (15:0)
  \item UCD.Lipid.499 -4 -4 FA (15:1)
  \item UCD.Lipid.500 -4 -4 FA (16:0)
  \item UCD.Lipid.501 -4 -4 FA (16:1)
  \item UCD.Lipid.502 -4 -4 FA (17:0)
  \item UCD.Lipid.503 -4 -4 FA (17:1)
  \item UCD.Lipid.504 -4 -4 FA (18:0)
  \item UCD.Lipid.505 -4 -4 FA (18:1)
  \item UCD.Lipid.506 -4 -4 FA (18:2)
  \item UCD.Lipid.507 -4 -4 FA (18:3)
  \item UCD.Lipid.508 -4 -4 FA (20:0)
  \item UCD.Lipid.509 -4 -4 FA (20:1)
  \item UCD.Lipid.510 -4 -4 FA (20:2)
  \item UCD.Lipid.511 -4 -4 FA (20:3)
  \item UCD.Lipid.512 -4 -4 FA (20:4)
  \item UCD.Lipid.513 -4 -4 FA (20:5)
  \item UCD.Lipid.514 -4 -4 FA (22:0)
  \item UCD.Lipid.515 -4 -4 FA (22:1)
  \item UCD.Lipid.516 -4 -4 FA (22:2)
  \item UCD.Lipid.517 -4 -4 FA (22:6)
  \item UCD.Lipid.518 -4 -4 FA (24:0)
  \item UCD.Lipid.519 -4 -4 FA (24:1)
  \item UCD.Lipid.520 -4 -4 FA (26:0)
  \item UCD.Lipid.521 -4 -4 FA (28:0)
}
}
\examples{
\dontrun{
describe(admclipidomics)
}
}
\seealso{
\url{https://adni.bitbucket.io/reference/docs/ADMCLIPIDOMICS/ADMC\%20Lipidomics\%20methods.pdf}
}
\keyword{datasets}
