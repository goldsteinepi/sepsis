#################
# NICU sepsis analysis
# Citation: Goldstein ND, Eppes SC, Ingraham BC, Paul DA. Characteristics of late-onset sepsis in the NICU: Does occupancy impact risk of infection? J Perinatol. 2016 Sep;36(9):753-7.
# 4/7/15 -- Neal Goldstein
#################


### FUNCTIONS ###

library(psych) #describe, describeBy
library(gmodels) #CrossTable
library(epitools) #rate ratio
library(survival) #cox regression
library(RColorBrewer) #color palette

season = function(event_date)
{
  season=NA
  
  if (!is.na(event_date))
  {
    month = as.numeric(format(as.Date(event_date), format = "%m"))
    day = as.numeric(format(as.Date(event_date), format = "%d"))
    
    if ((month==12 && day>=21) || (month==1) || (month==2) || (month==3 && day<=20)) {
      season="Winter"
    } else if ((month==3 && day>=21) || (month==4) || (month==5) || (month==6 && day<=20)) {
      season="Spring"
    } else if ((month==6 && day>=21) || (month==7) || (month==8) || (month==9 && day<=20)) {
      season="Summer"
    } else if ((month==9 && day>=21) || (month==10) || (month==11) || (month==12 && day<=20)) {
      season="Fall"
    }
  }

  return(season)
}


### READ DATA ###

load("NICU.2015-12-28.RData")


### RECODE ###

#black or hispanic
NICU$Race_ethnicity = ifelse(NICU$Mom_race_ethnicity==1, 1, ifelse(NICU$Mom_race_ethnicity==2, 2, 0))

#preterm
NICU$Preterm = ifelse(NICU$Gestational_age<37, 1, 0)

#preterm by other gestations
NICU$Preterm_25 = ifelse(NICU$Gestational_age<25, 1, 0)
NICU$Preterm_28 = ifelse(NICU$Gestational_age<28, 1, 0)
NICU$Preterm_32 = ifelse(NICU$Gestational_age<32, 1, 0)

#vlbw, 1500g
NICU$Birthweight_1500 = ifelse(NICU$Birthweight<1500, 1, 0)

#LOS quantile
NICU$LOS_quantile = ifelse(NICU$LOS<quantile(NICU$LOS, na.rm=T, probs=seq(0,1,0.25))[2], 0, ifelse(NICU$LOS<quantile(NICU$LOS, na.rm=T, probs=seq(0,1,0.25))[3], 1, ifelse(NICU$LOS<quantile(NICU$LOS, na.rm=T, probs=seq(0,1,0.25))[4], 2, 3)))

#LOS log transform
NICU$LOS_log = log(NICU$LOS)

#month of sepsis
NICU$Sepsis_onset_month = as.numeric(format(as.Date(NICU$Sepsis_onset_date), format = "%m"))

#year of sepsis
NICU$Sepsis_onset_year = as.numeric(format(as.Date(NICU$Sepsis_onset_date), format = "%Y"))

#clustering unit for random effects model
NICU$Cluster = as.numeric(format(as.Date(NICU$Date_admission), format = "%m%Y"))

#manually code sepsis organisms that were not in NIS
NICU$Sepsis_gram_pos_cocci[NICU$ID==3] = 1
NICU$Sepsis_staph_aureus[NICU$ID==157] = 1
NICU$Sepsis_strep_alpha[NICU$ID==440] = 1
NICU$Sepsis_coag_neg_staph[NICU$ID==489] = 1
NICU$Sepsis_GBS[NICU$ID==603] = 1
NICU$Sepsis_gram_pos_cocci[NICU$ID==860] = 1
NICU$Sepsis_GBS[NICU$ID==1569] = 1
NICU$Sepsis_coag_neg_staph[NICU$ID==27500] = 1
NICU$Sepsis_pseudomonas[NICU$ID==27500] = 1
NICU$Sepsis_GBS[NICU$ID==10306] = 1
NICU$Sepsis_coag_neg_staph[NICU$ID==11658] = 1
NICU$Sepsis_e_coli[NICU$ID==12854] = 1
NICU$Sepsis_coag_neg_staph[NICU$ID==14277] = 1
NICU$Sepsis_coag_neg_staph[NICU$ID==14560] = 1
NICU$Sepsis_coag_neg_staph[NICU$ID==15206] = 1
NICU$Sepsis_e_coli[NICU$ID==15525] = 1
NICU$Sepsis_gram_pos_cocci[NICU$ID==16137] = 1
NICU$Sepsis_coag_neg_staph[NICU$ID==16473] = 1
NICU$Sepsis_coag_neg_staph[NICU$ID==18385] = 1
NICU$Sepsis_coag_neg_staph[NICU$ID==19450] = 1
NICU$Sepsis_GBS[NICU$ID==20789] = 1
NICU$Sepsis_coag_neg_staph[NICU$ID==21347] = 1
NICU$Sepsis_coag_neg_staph[NICU$ID==21408] = 1
NICU$Sepsis_staph_aureus[NICU$ID==23160] = 1
NICU$Sepsis_klebsiella[NICU$ID==23160] = 1
NICU$Sepsis_coag_neg_staph[NICU$ID==23215] = 1
NICU$Sepsis_enterobacter[NICU$ID==23215] = 1
NICU$Sepsis_strep_alpha[NICU$ID==23927] = 1
NICU$Sepsis_staph_aureus[NICU$ID==23939] = 1
NICU$Sepsis_coag_neg_staph[NICU$ID==24381] = 1
NICU$Sepsis_GBS[NICU$ID==24907] = 1
NICU$Sepsis_coag_neg_staph[NICU$ID==25064] = 1
NICU$Sepsis_coag_neg_staph[NICU$ID==25077] = 1
NICU$Sepsis_staph_aureus[NICU$ID==25383] = 1
NICU$Sepsis_klebsiella[NICU$ID==25383] = 1
NICU$Sepsis_GBS[NICU$ID==25390] = 1
NICU$Sepsis_klebsiella[NICU$ID==25390] = 1
NICU$Sepsis_GBS[NICU$ID==25432] = 1
NICU$Sepsis_e_coli[NICU$ID==25768] = 1
NICU$Sepsis_coag_neg_staph[NICU$ID==25933] = 1
NICU$Sepsis_GBS[NICU$ID==26201] = 1
NICU$Sepsis_gram_pos_cocci[NICU$ID==11804] = 1
NICU$Sepsis_staph_aureus[NICU$ID==11804] = 1
NICU$Sepsis_staph_aureus[NICU$ID==13115] = 1
NICU$Sepsis_gram_pos_cocci[NICU$ID==14687] = 1
NICU$Sepsis_fungal[NICU$ID==14687] = 1
NICU$Sepsis_gram_pos_cocci[NICU$ID==19361] = 1
NICU$Sepsis_staph_aureus[NICU$ID==19361] = 1
NICU$Sepsis_coag_neg_staph[NICU$ID==9867] = 1
NICU$Sepsis_e_coli[NICU$ID==11137] = 1
NICU$Sepsis_gram_pos_cocci[NICU$ID==11248] = 1
NICU$Sepsis_gram_pos_cocci[NICU$ID==11892] = 1
NICU$Sepsis_coag_neg_staph[NICU$ID==13367] = 1
NICU$Sepsis_GBS[NICU$ID==23985] = 1
NICU$Sepsis_staph_aureus[NICU$ID==24172] = 0
NICU$Sepsis_coag_neg_staph[NICU$ID==24172] = 1
NICU$Sepsis_gram_pos_cocci[NICU$ID==14837] = 1

#remove sepsis dx that were not confirmed by chart review
NICU$Sepsis[NICU$ID==1053] = 0
NICU$Sepsis[NICU$ID==1114] = 0
NICU$Sepsis[NICU$ID==1130] = 0
NICU$Sepsis[NICU$ID==1349] = 0
NICU$Sepsis[NICU$ID==1386] = 0
NICU$Sepsis[NICU$ID==1608] = 0
NICU$Sepsis[NICU$ID==1622] = 0
NICU$Sepsis[NICU$ID==1866] = 0
NICU$Sepsis[NICU$ID==9359] = 0
NICU$Sepsis[NICU$ID==9874] = 0
NICU$Sepsis[NICU$ID==10088] = 0
NICU$Sepsis[NICU$ID==10703] = 0
NICU$Sepsis[NICU$ID==13943] = 0
NICU$Sepsis[NICU$ID==14638] = 0
NICU$Sepsis[NICU$ID==14837] = 0
NICU$Sepsis[NICU$ID==16301] = 0
NICU$Sepsis[NICU$ID==23910] = 0
NICU$Sepsis[NICU$ID==23969] = 0
NICU$Sepsis[NICU$ID==24829] = 0
NICU$Sepsis[NICU$ID==25132] = 0
NICU$Sepsis[NICU$ID==25235] = 0
NICU$Sepsis[NICU$ID==25322] = 0
NICU$Sepsis[NICU$ID==25412] = 0
NICU$Sepsis[NICU$ID==25664] = 0
NICU$Sepsis[NICU$ID==25879] = 0
NICU$Sepsis[NICU$ID==26061] = 0
NICU$Sepsis[NICU$ID==26761] = 0
NICU$Sepsis[NICU$ID==27020] = 0
NICU$Sepsis[NICU$ID==27054] = 0
NICU$Sepsis[NICU$ID==27062] = 0
NICU$Sepsis[NICU$ID==27240] = 0
NICU$Sepsis[NICU$ID==27300] = 0
NICU$Sepsis[NICU$ID==27324] = 0
NICU$Sepsis[NICU$ID==1349] = 0
NICU$Sepsis[NICU$ID==1386] = 0
NICU$Sepsis[NICU$ID==9359] = 0
NICU$Sepsis[NICU$ID==10088] = 0
NICU$Sepsis[NICU$ID==23910] = 0
NICU$Sepsis[NICU$ID==24829] = 0
NICU$Sepsis[NICU$ID==25132] = 0
NICU$Sepsis[NICU$ID==25235] = 0
NICU$Sepsis[NICU$ID==26061] = 0

#collapse organisms to categories, give priority to S. aureus, GBS, CoNS, E. coli
#numbers force order in barplot
NICU$Sepsis_cause = ifelse(NICU$Sepsis_staph_aureus==1, "8 Gram+, S. aureus", ifelse(NICU$Sepsis_GBS==1, "7 Gram+, GBS", ifelse(NICU$Sepsis_coag_neg_staph==1, "6 Gram+, CoNS", ifelse(NICU$Sepsis_strep_alpha==1 | NICU$Sepsis_GAS==1 | NICU$Sepsis_gram_pos_cocci==1 | NICU$Sepsis_listeria==1, "5 Gram+, other", ifelse(NICU$Sepsis_e_coli==1, "4 Gram-, E.coli", ifelse(NICU$Sepsis_acinetobacter==1 | NICU$Sepsis_enterobacter==1 | NICU$Sepsis_klebsiella==1 | NICU$Sepsis_pseudomonas==1 | NICU$Sepsis_salmonella==1 | NICU$Sepsis_serratia==1, "3 Gram-, other", ifelse(NICU$Sepsis_fungal==1, "2 Fungal", ifelse(NICU$Sepsis==1, "1 Other", NA))))))))

#manually code sepsis categories for uncommon organisms
NICU$Sepsis_cause[NICU$ID==15574] =  "3 Gram-, other"
NICU$Sepsis_cause[NICU$ID==16704] =  "3 Gram-, other"
NICU$Sepsis_cause[NICU$ID==17189] =  "3 Gram-, other"
NICU$Sepsis_cause[NICU$ID==20203] =  "3 Gram-, other"
NICU$Sepsis_cause[NICU$ID==23901] =  "5 Gram+, other"
NICU$Sepsis_cause[NICU$ID==26132] =  "3 Gram-, other"
NICU$Sepsis_cause[NICU$ID==27287] =  "3 Gram-, other"
NICU$Sepsis_cause[NICU$ID==10040] =  "3 Gram-, other"
NICU$Sepsis_cause[NICU$ID==13574] =  "3 Gram-, other"
NICU$Sepsis_cause[NICU$ID==14315] =  "3 Gram-, other"
NICU$Sepsis_cause[NICU$ID==10703] =  "3 Gram-, other"
NICU$Sepsis_cause[NICU$ID==13943] =  "3 Gram-, other"

#MSSA/MRSA
NICU$Sepsis_MRSA = NA
NICU$Sepsis_MRSA[NICU$ID==11804] = "2 MRSA"
NICU$Sepsis_MRSA[NICU$ID==13115] = "2 MRSA"
NICU$Sepsis_MRSA[NICU$ID==468] = "2 MRSA"
NICU$Sepsis_MRSA[NICU$ID==549] = "2 MRSA"
NICU$Sepsis_MRSA[NICU$ID==11018] = "2 MRSA"
NICU$Sepsis_MRSA[NICU$ID==12129] = "2 MRSA"
NICU$Sepsis_MRSA[NICU$ID==12285] = "2 MRSA"
NICU$Sepsis_MRSA[NICU$ID==15335] = "2 MRSA"
NICU$Sepsis_MRSA[NICU$ID==18302] = "2 MRSA"
NICU$Sepsis_MRSA[NICU$ID==21561] = "2 MRSA"
NICU$Sepsis_MRSA[NICU$ID==24048] = "2 MRSA"
NICU$Sepsis_MRSA[NICU$ID==24056] = "2 MRSA"
NICU$Sepsis_MRSA[NICU$ID==24191] = "2 MRSA"
NICU$Sepsis_MRSA[NICU$ID==25407] = "2 MRSA"
NICU$Sepsis_MRSA[NICU$ID==21233] = "2 MRSA"
NICU$Sepsis_MRSA[NICU$ID==11439] = "2 MRSA"
NICU$Sepsis_MRSA[NICU$ID==11876] = "2 MRSA"
NICU$Sepsis_MRSA[NICU$ID==12300] = "2 MRSA"
NICU$Sepsis_MRSA[NICU$ID==12655] = "2 MRSA"
NICU$Sepsis_MRSA[NICU$ID==157] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==23160] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==23939] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==25383] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==19361] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==9721] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==9817] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==10201] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==10497] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==10576] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==10717] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==11874] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==11985] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==12041] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==14169] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==14246] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==14377] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==14681] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==15014] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==17114] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==18163] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==18466] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==20899] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==21307] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==21822] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==23017] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==23429] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==24025] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==24234] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==24669] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==24830] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==25765] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==26002] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==27177] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==20510] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==11331] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==13819] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==14234] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==12202] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==14832] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==15016] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==23967] = "3 MSSA"
NICU$Sepsis_MRSA[NICU$ID==24380] = "3 MSSA"
NICU$Sepsis_MRSA = ifelse(NICU$Sepsis_staph_aureus==1 & is.na(NICU$Sepsis_MRSA), "1 Unknown", NICU$Sepsis_MRSA)

#late onset sepsis
NICU$Sepsis_late = ifelse(NICU$Sepsis==1 & NICU$Sepsis_onset==1, 1, 0)

#sepsis dx during intial admission
#NICU$Sepsis_initial_admission = ifelse(NICU$Sepsis==1 & NICU$Sepsis_onset_date<=NICU$Date_discharge_initial, 1, ifelse(NICU$Sepsis==1 & NICU$Sepsis_onset_date>NICU$Date_discharge_initial, 0, NA))
#NICU$Sepsis_initial_admission = ifelse(NICU$Sepsis==1 & NICU$Admission_n==1, 1, NICU$Sepsis_initial_admission)

#add variables that require individual handling
NICU$Admission_season = NA
NICU$Sepsis_onset_season = NA
NICU$Sepsis_onset_seasonyear = NA
NICU$Cohort_time = NA
NICU$Sepsis_stay = NA
NICU$Census_45_n = NA
NICU$Census_50_n = NA
NICU$Census_55_n = NA
NICU$Census_60_n = NA
NICU$Census_65_n = NA
NICU$Census_preterm25_cumulative = NA
NICU$Census_preterm28_cumulative = NA
NICU$Census_preterm32_cumulative = NA
NICU$Census_cumulative = NA

#census for each day, up to max(LOS)
# NICU$Census_daily1 = NA
# NICU$Census_daily2 = NA
# NICU$Census_daily3 = NA
# NICU$Census_daily4 = NA
# NICU$Census_daily5 = NA
# NICU$Census_daily6 = NA
# NICU$Census_daily7 = NA
# NICU$Census_daily8 = NA
# NICU$Census_daily9 = NA
# NICU$Census_daily10 = NA
# NICU$Census_daily11 = NA
# NICU$Census_daily12 = NA
# NICU$Census_daily13 = NA
# NICU$Census_daily14 = NA
# NICU$Census_daily15 = NA
# NICU$Census_daily16 = NA
# NICU$Census_daily17 = NA
# NICU$Census_daily18 = NA
# NICU$Census_daily19 = NA
# NICU$Census_daily20 = NA
# NICU$Census_daily21 = NA
# NICU$Census_daily22 = NA
# NICU$Census_daily23 = NA
# NICU$Census_daily24 = NA
# NICU$Census_daily25 = NA
# NICU$Census_daily26 = NA
# NICU$Census_daily27 = NA
# NICU$Census_daily28 = NA
# NICU$Census_daily29 = NA
# NICU$Census_daily30 = NA
# NICU$Census_daily31 = NA
# NICU$Census_daily32 = NA
# NICU$Census_daily33 = NA
# NICU$Census_daily34 = NA
# NICU$Census_daily35 = NA
# NICU$Census_daily36 = NA
# NICU$Census_daily37 = NA
# NICU$Census_daily38 = NA
# NICU$Census_daily39 = NA
# NICU$Census_daily40 = NA
# NICU$Census_daily41 = NA
# NICU$Census_daily42 = NA
# NICU$Census_daily43 = NA
# NICU$Census_daily44 = NA
# NICU$Census_daily45 = NA
# NICU$Census_daily46 = NA
# NICU$Census_daily47 = NA
# NICU$Census_daily48 = NA
# NICU$Census_daily49 = NA
# NICU$Census_daily50 = NA
# NICU$Census_daily51 = NA
# NICU$Census_daily52 = NA
# NICU$Census_daily53 = NA
# NICU$Census_daily54 = NA
# NICU$Census_daily55 = NA
# NICU$Census_daily56 = NA
# NICU$Census_daily57 = NA
# NICU$Census_daily58 = NA
# NICU$Census_daily59 = NA
# NICU$Census_daily60 = NA
# NICU$Census_daily61 = NA
# NICU$Census_daily62 = NA
# NICU$Census_daily63 = NA
# NICU$Census_daily64 = NA
# NICU$Census_daily65 = NA
# NICU$Census_daily66 = NA
# NICU$Census_daily67 = NA
# NICU$Census_daily68 = NA
# NICU$Census_daily69 = NA
# NICU$Census_daily70 = NA
# NICU$Census_daily71 = NA
# NICU$Census_daily72 = NA
# NICU$Census_daily73 = NA
# NICU$Census_daily74 = NA
# NICU$Census_daily75 = NA
# NICU$Census_daily76 = NA
# NICU$Census_daily77 = NA
# NICU$Census_daily78 = NA
# NICU$Census_daily79 = NA
# NICU$Census_daily80 = NA
# NICU$Census_daily81 = NA
# NICU$Census_daily82 = NA
# NICU$Census_daily83 = NA
# NICU$Census_daily84 = NA
# NICU$Census_daily85 = NA
# NICU$Census_daily86 = NA
# NICU$Census_daily87 = NA
# NICU$Census_daily88 = NA
# NICU$Census_daily89 = NA
# NICU$Census_daily90 = NA
# NICU$Census_daily91 = NA
# NICU$Census_daily92 = NA
# NICU$Census_daily93 = NA
# NICU$Census_daily94 = NA
# NICU$Census_daily95 = NA
# NICU$Census_daily96 = NA
# NICU$Census_daily97 = NA
# NICU$Census_daily98 = NA
# NICU$Census_daily99 = NA
# NICU$Census_daily100 = NA
# NICU$Census_daily101 = NA
# NICU$Census_daily102 = NA
# NICU$Census_daily103 = NA
# NICU$Census_daily104 = NA
# NICU$Census_daily105 = NA
# NICU$Census_daily106 = NA
# NICU$Census_daily107 = NA
# NICU$Census_daily108 = NA
# NICU$Census_daily109 = NA
# NICU$Census_daily110 = NA
# NICU$Census_daily111 = NA
# NICU$Census_daily112 = NA
# NICU$Census_daily113 = NA
# NICU$Census_daily114 = NA
# NICU$Census_daily115 = NA
# NICU$Census_daily116 = NA
# NICU$Census_daily117 = NA
# NICU$Census_daily118 = NA
# NICU$Census_daily119 = NA
# NICU$Census_daily120 = NA
# NICU$Census_daily121 = NA
# NICU$Census_daily122 = NA
# NICU$Census_daily123 = NA
# NICU$Census_daily124 = NA
# NICU$Census_daily125 = NA
# NICU$Census_daily126 = NA
# NICU$Census_daily127 = NA
# NICU$Census_daily128 = NA
# NICU$Census_daily129 = NA
# NICU$Census_daily130 = NA
# NICU$Census_daily131 = NA
# NICU$Census_daily132 = NA
# NICU$Census_daily133 = NA
# NICU$Census_daily134 = NA
# NICU$Census_daily135 = NA
# NICU$Census_daily136 = NA
# NICU$Census_daily137 = NA
# NICU$Census_daily138 = NA
# NICU$Census_daily139 = NA
# NICU$Census_daily140 = NA
# NICU$Census_daily141 = NA
# NICU$Census_daily142 = NA
# NICU$Census_daily143 = NA
# NICU$Census_daily144 = NA
# NICU$Census_daily145 = NA
# NICU$Census_daily146 = NA
# NICU$Census_daily147 = NA
# NICU$Census_daily148 = NA
# NICU$Census_daily149 = NA
# NICU$Census_daily150 = NA
# NICU$Census_daily151 = NA
# NICU$Census_daily152 = NA
# NICU$Census_daily153 = NA
# NICU$Census_daily154 = NA
# NICU$Census_daily155 = NA
# NICU$Census_daily156 = NA
# NICU$Census_daily157 = NA
# NICU$Census_daily158 = NA
# NICU$Census_daily159 = NA
# NICU$Census_daily160 = NA
# NICU$Census_daily161 = NA
# NICU$Census_daily162 = NA
# NICU$Census_daily163 = NA
# NICU$Census_daily164 = NA
# NICU$Census_daily165 = NA
# NICU$Census_daily166 = NA
# NICU$Census_daily167 = NA
# NICU$Census_daily168 = NA
# NICU$Census_daily169 = NA
# NICU$Census_daily170 = NA
# NICU$Census_daily171 = NA
# NICU$Census_daily172 = NA
# NICU$Census_daily173 = NA
# NICU$Census_daily174 = NA
# NICU$Census_daily175 = NA
# NICU$Census_daily176 = NA
# NICU$Census_daily177 = NA
# NICU$Census_daily178 = NA
# NICU$Census_daily179 = NA
# NICU$Census_daily180 = NA
# NICU$Census_daily181 = NA
# NICU$Census_daily182 = NA
# NICU$Census_daily183 = NA
# NICU$Census_daily184 = NA
# NICU$Census_daily185 = NA
# NICU$Census_daily186 = NA
# NICU$Census_daily187 = NA
# NICU$Census_daily188 = NA
# NICU$Census_daily189 = NA
# NICU$Census_daily190 = NA
# NICU$Census_daily191 = NA
# NICU$Census_daily192 = NA
# NICU$Census_daily193 = NA
# NICU$Census_daily194 = NA
# NICU$Census_daily195 = NA
# NICU$Census_daily196 = NA
# NICU$Census_daily197 = NA
# NICU$Census_daily198 = NA
# NICU$Census_daily199 = NA
# NICU$Census_daily200 = NA
# NICU$Census_daily201 = NA
# NICU$Census_daily202 = NA
# NICU$Census_daily203 = NA
# NICU$Census_daily204 = NA
# NICU$Census_daily205 = NA
# NICU$Census_daily206 = NA
# NICU$Census_daily207 = NA
# NICU$Census_daily208 = NA
# NICU$Census_daily209 = NA
# NICU$Census_daily210 = NA
# NICU$Census_daily211 = NA
# NICU$Census_daily212 = NA
# NICU$Census_daily213 = NA
# NICU$Census_daily214 = NA
# NICU$Census_daily215 = NA
# NICU$Census_daily216 = NA
# NICU$Census_daily217 = NA
# NICU$Census_daily218 = NA
# NICU$Census_daily219 = NA
# NICU$Census_daily220 = NA
# NICU$Census_daily221 = NA
# NICU$Census_daily222 = NA
# NICU$Census_daily223 = NA
# NICU$Census_daily224 = NA
# NICU$Census_daily225 = NA
# NICU$Census_daily226 = NA
# NICU$Census_daily227 = NA
# NICU$Census_daily228 = NA
# NICU$Census_daily229 = NA
# NICU$Census_daily230 = NA
# NICU$Census_daily231 = NA
# NICU$Census_daily232 = NA
# NICU$Census_daily233 = NA
# NICU$Census_daily234 = NA
# NICU$Census_daily235 = NA
# NICU$Census_daily236 = NA
# NICU$Census_daily237 = NA
# NICU$Census_daily238 = NA
# NICU$Census_daily239 = NA
# NICU$Census_daily240 = NA
# NICU$Census_daily241 = NA
# NICU$Census_daily242 = NA
# NICU$Census_daily243 = NA
# NICU$Census_daily244 = NA
# NICU$Census_daily245 = NA
# NICU$Census_daily246 = NA
# NICU$Census_daily247 = NA
# NICU$Census_daily248 = NA
# NICU$Census_daily249 = NA
# NICU$Census_daily250 = NA
# NICU$Census_daily251 = NA
# NICU$Census_daily252 = NA
# NICU$Census_daily253 = NA
# NICU$Census_daily254 = NA
# NICU$Census_daily255 = NA
# NICU$Census_daily256 = NA
# NICU$Census_daily257 = NA
# NICU$Census_daily258 = NA
# NICU$Census_daily259 = NA
# NICU$Census_daily260 = NA
# NICU$Census_daily261 = NA
# NICU$Census_daily262 = NA
# NICU$Census_daily263 = NA
# NICU$Census_daily264 = NA
# NICU$Census_daily265 = NA
# NICU$Census_daily266 = NA
# NICU$Census_daily267 = NA
# NICU$Census_daily268 = NA
# NICU$Census_daily269 = NA
# NICU$Census_daily270 = NA
# NICU$Census_daily271 = NA
# NICU$Census_daily272 = NA
# NICU$Census_daily273 = NA
# NICU$Census_daily274 = NA
# NICU$Census_daily275 = NA
# NICU$Census_daily276 = NA
# NICU$Census_daily277 = NA
# NICU$Census_daily278 = NA
# NICU$Census_daily279 = NA
# NICU$Census_daily280 = NA
# NICU$Census_daily281 = NA
# NICU$Census_daily282 = NA
# NICU$Census_daily283 = NA
# NICU$Census_daily284 = NA
# NICU$Census_daily285 = NA
# NICU$Census_daily286 = NA
# NICU$Census_daily287 = NA
# NICU$Census_daily288 = NA
# NICU$Census_daily289 = NA
# NICU$Census_daily290 = NA
# NICU$Census_daily291 = NA
# NICU$Census_daily292 = NA
# NICU$Census_daily293 = NA
# NICU$Census_daily294 = NA
# NICU$Census_daily295 = NA
# NICU$Census_daily296 = NA
# NICU$Census_daily297 = NA
# NICU$Census_daily298 = NA
# NICU$Census_daily299 = NA
# NICU$Census_daily300 = NA
# NICU$Census_daily301 = NA
# NICU$Census_daily302 = NA
# NICU$Census_daily303 = NA
# NICU$Census_daily304 = NA
# NICU$Census_daily305 = NA
# NICU$Census_daily306 = NA
# NICU$Census_daily307 = NA
# NICU$Census_daily308 = NA
# NICU$Census_daily309 = NA
# NICU$Census_daily310 = NA
# NICU$Census_daily311 = NA
# NICU$Census_daily312 = NA
# NICU$Census_daily313 = NA
# NICU$Census_daily314 = NA
# NICU$Census_daily315 = NA
# NICU$Census_daily316 = NA
# NICU$Census_daily317 = NA
# NICU$Census_daily318 = NA
# NICU$Census_daily319 = NA
# NICU$Census_daily320 = NA
# NICU$Census_daily321 = NA
# NICU$Census_daily322 = NA
# NICU$Census_daily323 = NA
# NICU$Census_daily324 = NA
# NICU$Census_daily325 = NA
# NICU$Census_daily326 = NA
# NICU$Census_daily327 = NA
# NICU$Census_daily328 = NA
# NICU$Census_daily329 = NA
# NICU$Census_daily330 = NA
# NICU$Census_daily331 = NA
# NICU$Census_daily332 = NA
# NICU$Census_daily333 = NA
# NICU$Census_daily334 = NA
# NICU$Census_daily335 = NA
# NICU$Census_daily336 = NA
# NICU$Census_daily337 = NA
# NICU$Census_daily338 = NA
# NICU$Census_daily339 = NA
# NICU$Census_daily340 = NA
# NICU$Census_daily341 = NA
# NICU$Census_daily342 = NA
# NICU$Census_daily343 = NA
# NICU$Census_daily344 = NA
# NICU$Census_daily345 = NA
# NICU$Census_daily346 = NA
# NICU$Census_daily347 = NA
# NICU$Census_daily348 = NA
# NICU$Census_daily349 = NA
# NICU$Census_daily350 = NA
# NICU$Census_daily351 = NA
# NICU$Census_daily352 = NA
# NICU$Census_daily353 = NA
# NICU$Census_daily354 = NA
# NICU$Census_daily355 = NA
# NICU$Census_daily356 = NA
# NICU$Census_daily357 = NA
# NICU$Census_daily358 = NA
# NICU$Census_daily359 = NA
# NICU$Census_daily360 = NA
# NICU$Census_daily361 = NA
# NICU$Census_daily362 = NA
# NICU$Census_daily363 = NA
# NICU$Census_daily364 = NA
# NICU$Census_daily365 = NA
# NICU$Census_daily366 = NA
# NICU$Census_daily367 = NA
# NICU$Census_daily368 = NA
# NICU$Census_daily369 = NA
# NICU$Census_daily370 = NA
# NICU$Census_daily371 = NA
# NICU$Census_daily372 = NA
# NICU$Census_daily373 = NA
# NICU$Census_daily374 = NA
# NICU$Census_daily375 = NA
# NICU$Census_daily376 = NA
# NICU$Census_daily377 = NA
# NICU$Census_daily378 = NA
# NICU$Census_daily379 = NA
# NICU$Census_daily380 = NA
# NICU$Census_daily381 = NA
# NICU$Census_daily382 = NA
# NICU$Census_daily383 = NA
# NICU$Census_daily384 = NA
# NICU$Census_daily385 = NA
# NICU$Census_daily386 = NA
# NICU$Census_daily387 = NA
# NICU$Census_daily388 = NA
# NICU$Census_daily389 = NA
# NICU$Census_daily390 = NA
# NICU$Census_daily391 = NA
# NICU$Census_daily392 = NA
# NICU$Census_daily393 = NA
# NICU$Census_daily394 = NA
# NICU$Census_daily395 = NA
# NICU$Census_daily396 = NA
# NICU$Census_daily397 = NA
# NICU$Census_daily398 = NA
# NICU$Census_daily399 = NA
# NICU$Census_daily400 = NA
# NICU$Census_daily401 = NA
# NICU$Census_daily402 = NA
# NICU$Census_daily403 = NA
# NICU$Census_daily404 = NA
# NICU$Census_daily405 = NA
# NICU$Census_daily406 = NA
# NICU$Census_daily407 = NA
# NICU$Census_daily408 = NA

for (i in 1:nrow(NICU))
  {
  cat("\n\n************** ","Observation: ",i," **************\n",sep="")
  
  #seasonality of admission
  NICU$Admission_season[i] = season(NICU$Date_admission[i])

  #seasonality of sepsis
  NICU$Sepsis_onset_season[i] = season(NICU$Sepsis_onset_date[i])
  NICU$Sepsis_onset_seasonyear[i] = NICU$Sepsis_onset_year[i] + ifelse(NICU$Sepsis_onset_season[i]=="Spring", 0.25, ifelse(NICU$Sepsis_onset_season[i]=="Summer", 0.5, ifelse(NICU$Sepsis_onset_season[i]=="Fall", 0.75, 0)))
  
  #date of outcome or censoring
  if (NICU$Sepsis_late[i]==1)
  {
    followup = NICU$Sepsis_onset_date[i]
  } else {
    followup = NICU$Date_discharge_initial[i]
  }
  
  #cohort time to event: outcome and censoring
  NICU$Cohort_time[i] = as.numeric(followup - NICU$Date_admission[i]) + 1

  #sepsis during stay (subtract current case, if there is one)
  if (NICU$Sepsis[i]==0) {
    NICU$Sepsis_stay[i] = ifelse(sum((NICU$Sepsis_onset_date >= NICU$Date_admission[i]) & (NICU$Sepsis_onset_date <= followup), na.rm=T)>0, 1, 0)
  } else {
    NICU$Sepsis_stay[i] = ifelse(sum((NICU$Sepsis_onset_date >= NICU$Date_admission[i]) & (NICU$Sepsis_onset_date <= followup), na.rm=T)>1, 1, 0)
  }

  #operationalize census variables
  if ((!is.na(NICU$Date_admission[i])) && (!is.na(followup)))
  {
    #reset daily census vars based on max LOS
#     for (j in 1:max(NICU$LOS))
#       assign(paste("Census_daily",j,sep=""),NA)
#     rm(j)
    
    count_census = 0
    count_45 = 0
    count_50 = 0
    count_55 = 0
    count_60 = 0
    count_65 = 0
    count_preterm25 = 0
    count_preterm28 = 0
    count_preterm32 = 0
    
    for (j in 0:as.numeric(followup-NICU$Date_admission[i]))
    {
      #count daily census
      daily_census = as.numeric(sum((NICU$Date_admission <= (NICU$Date_admission[i]+j)) & (NICU$Date_discharge_initial >= (NICU$Date_admission[i]+j)), na.rm=T))
      
      #count daily preterms
      daily_preterms25 = as.numeric(sum(NICU$Preterm_25[which((NICU$Date_admission <= (NICU$Date_admission[i]+j)) & (NICU$Date_discharge_initial >= (NICU$Date_admission[i]+j)))], na.rm=T))
      daily_preterms28 = as.numeric(sum(NICU$Preterm_28[which((NICU$Date_admission <= (NICU$Date_admission[i]+j)) & (NICU$Date_discharge_initial >= (NICU$Date_admission[i]+j)))], na.rm=T))
      daily_preterms32 = as.numeric(sum(NICU$Preterm_32[which((NICU$Date_admission <= (NICU$Date_admission[i]+j)) & (NICU$Date_discharge_initial >= (NICU$Date_admission[i]+j)))], na.rm=T))
      
      #cumulative total of census per day
      count_census = count_census + daily_census
      
      #cumulative total of preterm per day
      count_preterm25 = count_preterm25 + daily_preterms25
      count_preterm28 = count_preterm28 + daily_preterms28
      count_preterm32 = count_preterm32 + daily_preterms32
      
      #check for number of days exceeding thresholds
      count_45 = count_45 + ifelse(daily_census>=45, 1, 0)
      count_50 = count_50 + ifelse(daily_census>=50, 1, 0)
      count_55 = count_55 + ifelse(daily_census>=55, 1, 0)
      count_60 = count_60 + ifelse(daily_census>=60, 1, 0)
      count_65 = count_65 + ifelse(daily_census>=65, 1, 0)

      #compute proportion of days in NICU with high census
#       assign(paste("Census_daily",j+1,sep=""),(count_55/(j+1)))
    }
    
    #store census counts
    NICU$Census_45_n[i] = count_45
    NICU$Census_50_n[i] = count_50
    NICU$Census_55_n[i] = count_55
    NICU$Census_60_n[i] = count_60
    NICU$Census_65_n[i] = count_65
    NICU$Census_preterm25_cumulative[i] = count_preterm25
    NICU$Census_preterm28_cumulative[i] = count_preterm28
    NICU$Census_preterm32_cumulative[i] = count_preterm32
    NICU$Census_cumulative[i] = count_census
    
    #census for each day, up to max(LOS)
#     NICU$Census_daily1[i] = Census_daily1
#     NICU$Census_daily2[i] = Census_daily2
#     NICU$Census_daily3[i] = Census_daily3
#     NICU$Census_daily4[i] = Census_daily4
#     NICU$Census_daily5[i] = Census_daily5
#     NICU$Census_daily6[i] = Census_daily6
#     NICU$Census_daily7[i] = Census_daily7
#     NICU$Census_daily8[i] = Census_daily8
#     NICU$Census_daily9[i] = Census_daily9
#     NICU$Census_daily10[i] = Census_daily10
#     NICU$Census_daily11[i] = Census_daily11
#     NICU$Census_daily12[i] = Census_daily12
#     NICU$Census_daily13[i] = Census_daily13
#     NICU$Census_daily14[i] = Census_daily14
#     NICU$Census_daily15[i] = Census_daily15
#     NICU$Census_daily16[i] = Census_daily16
#     NICU$Census_daily17[i] = Census_daily17
#     NICU$Census_daily18[i] = Census_daily18
#     NICU$Census_daily19[i] = Census_daily19
#     NICU$Census_daily20[i] = Census_daily20
#     NICU$Census_daily21[i] = Census_daily21
#     NICU$Census_daily22[i] = Census_daily22
#     NICU$Census_daily23[i] = Census_daily23
#     NICU$Census_daily24[i] = Census_daily24
#     NICU$Census_daily25[i] = Census_daily25
#     NICU$Census_daily26[i] = Census_daily26
#     NICU$Census_daily27[i] = Census_daily27
#     NICU$Census_daily28[i] = Census_daily28
#     NICU$Census_daily29[i] = Census_daily29
#     NICU$Census_daily30[i] = Census_daily30
#     NICU$Census_daily31[i] = Census_daily31
#     NICU$Census_daily32[i] = Census_daily32
#     NICU$Census_daily33[i] = Census_daily33
#     NICU$Census_daily34[i] = Census_daily34
#     NICU$Census_daily35[i] = Census_daily35
#     NICU$Census_daily36[i] = Census_daily36
#     NICU$Census_daily37[i] = Census_daily37
#     NICU$Census_daily38[i] = Census_daily38
#     NICU$Census_daily39[i] = Census_daily39
#     NICU$Census_daily40[i] = Census_daily40
#     NICU$Census_daily41[i] = Census_daily41
#     NICU$Census_daily42[i] = Census_daily42
#     NICU$Census_daily43[i] = Census_daily43
#     NICU$Census_daily44[i] = Census_daily44
#     NICU$Census_daily45[i] = Census_daily45
#     NICU$Census_daily46[i] = Census_daily46
#     NICU$Census_daily47[i] = Census_daily47
#     NICU$Census_daily48[i] = Census_daily48
#     NICU$Census_daily49[i] = Census_daily49
#     NICU$Census_daily50[i] = Census_daily50
#     NICU$Census_daily51[i] = Census_daily51
#     NICU$Census_daily52[i] = Census_daily52
#     NICU$Census_daily53[i] = Census_daily53
#     NICU$Census_daily54[i] = Census_daily54
#     NICU$Census_daily55[i] = Census_daily55
#     NICU$Census_daily56[i] = Census_daily56
#     NICU$Census_daily57[i] = Census_daily57
#     NICU$Census_daily58[i] = Census_daily58
#     NICU$Census_daily59[i] = Census_daily59
#     NICU$Census_daily60[i] = Census_daily60
#     NICU$Census_daily61[i] = Census_daily61
#     NICU$Census_daily62[i] = Census_daily62
#     NICU$Census_daily63[i] = Census_daily63
#     NICU$Census_daily64[i] = Census_daily64
#     NICU$Census_daily65[i] = Census_daily65
#     NICU$Census_daily66[i] = Census_daily66
#     NICU$Census_daily67[i] = Census_daily67
#     NICU$Census_daily68[i] = Census_daily68
#     NICU$Census_daily69[i] = Census_daily69
#     NICU$Census_daily70[i] = Census_daily70
#     NICU$Census_daily71[i] = Census_daily71
#     NICU$Census_daily72[i] = Census_daily72
#     NICU$Census_daily73[i] = Census_daily73
#     NICU$Census_daily74[i] = Census_daily74
#     NICU$Census_daily75[i] = Census_daily75
#     NICU$Census_daily76[i] = Census_daily76
#     NICU$Census_daily77[i] = Census_daily77
#     NICU$Census_daily78[i] = Census_daily78
#     NICU$Census_daily79[i] = Census_daily79
#     NICU$Census_daily80[i] = Census_daily80
#     NICU$Census_daily81[i] = Census_daily81
#     NICU$Census_daily82[i] = Census_daily82
#     NICU$Census_daily83[i] = Census_daily83
#     NICU$Census_daily84[i] = Census_daily84
#     NICU$Census_daily85[i] = Census_daily85
#     NICU$Census_daily86[i] = Census_daily86
#     NICU$Census_daily87[i] = Census_daily87
#     NICU$Census_daily88[i] = Census_daily88
#     NICU$Census_daily89[i] = Census_daily89
#     NICU$Census_daily90[i] = Census_daily90
#     NICU$Census_daily91[i] = Census_daily91
#     NICU$Census_daily92[i] = Census_daily92
#     NICU$Census_daily93[i] = Census_daily93
#     NICU$Census_daily94[i] = Census_daily94
#     NICU$Census_daily95[i] = Census_daily95
#     NICU$Census_daily96[i] = Census_daily96
#     NICU$Census_daily97[i] = Census_daily97
#     NICU$Census_daily98[i] = Census_daily98
#     NICU$Census_daily99[i] = Census_daily99
#     NICU$Census_daily100[i] = Census_daily100
#     NICU$Census_daily101[i] = Census_daily101
#     NICU$Census_daily102[i] = Census_daily102
#     NICU$Census_daily103[i] = Census_daily103
#     NICU$Census_daily104[i] = Census_daily104
#     NICU$Census_daily105[i] = Census_daily105
#     NICU$Census_daily106[i] = Census_daily106
#     NICU$Census_daily107[i] = Census_daily107
#     NICU$Census_daily108[i] = Census_daily108
#     NICU$Census_daily109[i] = Census_daily109
#     NICU$Census_daily110[i] = Census_daily110
#     NICU$Census_daily111[i] = Census_daily111
#     NICU$Census_daily112[i] = Census_daily112
#     NICU$Census_daily113[i] = Census_daily113
#     NICU$Census_daily114[i] = Census_daily114
#     NICU$Census_daily115[i] = Census_daily115
#     NICU$Census_daily116[i] = Census_daily116
#     NICU$Census_daily117[i] = Census_daily117
#     NICU$Census_daily118[i] = Census_daily118
#     NICU$Census_daily119[i] = Census_daily119
#     NICU$Census_daily120[i] = Census_daily120
#     NICU$Census_daily121[i] = Census_daily121
#     NICU$Census_daily122[i] = Census_daily122
#     NICU$Census_daily123[i] = Census_daily123
#     NICU$Census_daily124[i] = Census_daily124
#     NICU$Census_daily125[i] = Census_daily125
#     NICU$Census_daily126[i] = Census_daily126
#     NICU$Census_daily127[i] = Census_daily127
#     NICU$Census_daily128[i] = Census_daily128
#     NICU$Census_daily129[i] = Census_daily129
#     NICU$Census_daily130[i] = Census_daily130
#     NICU$Census_daily131[i] = Census_daily131
#     NICU$Census_daily132[i] = Census_daily132
#     NICU$Census_daily133[i] = Census_daily133
#     NICU$Census_daily134[i] = Census_daily134
#     NICU$Census_daily135[i] = Census_daily135
#     NICU$Census_daily136[i] = Census_daily136
#     NICU$Census_daily137[i] = Census_daily137
#     NICU$Census_daily138[i] = Census_daily138
#     NICU$Census_daily139[i] = Census_daily139
#     NICU$Census_daily140[i] = Census_daily140
#     NICU$Census_daily141[i] = Census_daily141
#     NICU$Census_daily142[i] = Census_daily142
#     NICU$Census_daily143[i] = Census_daily143
#     NICU$Census_daily144[i] = Census_daily144
#     NICU$Census_daily145[i] = Census_daily145
#     NICU$Census_daily146[i] = Census_daily146
#     NICU$Census_daily147[i] = Census_daily147
#     NICU$Census_daily148[i] = Census_daily148
#     NICU$Census_daily149[i] = Census_daily149
#     NICU$Census_daily150[i] = Census_daily150
#     NICU$Census_daily151[i] = Census_daily151
#     NICU$Census_daily152[i] = Census_daily152
#     NICU$Census_daily153[i] = Census_daily153
#     NICU$Census_daily154[i] = Census_daily154
#     NICU$Census_daily155[i] = Census_daily155
#     NICU$Census_daily156[i] = Census_daily156
#     NICU$Census_daily157[i] = Census_daily157
#     NICU$Census_daily158[i] = Census_daily158
#     NICU$Census_daily159[i] = Census_daily159
#     NICU$Census_daily160[i] = Census_daily160
#     NICU$Census_daily161[i] = Census_daily161
#     NICU$Census_daily162[i] = Census_daily162
#     NICU$Census_daily163[i] = Census_daily163
#     NICU$Census_daily164[i] = Census_daily164
#     NICU$Census_daily165[i] = Census_daily165
#     NICU$Census_daily166[i] = Census_daily166
#     NICU$Census_daily167[i] = Census_daily167
#     NICU$Census_daily168[i] = Census_daily168
#     NICU$Census_daily169[i] = Census_daily169
#     NICU$Census_daily170[i] = Census_daily170
#     NICU$Census_daily171[i] = Census_daily171
#     NICU$Census_daily172[i] = Census_daily172
#     NICU$Census_daily173[i] = Census_daily173
#     NICU$Census_daily174[i] = Census_daily174
#     NICU$Census_daily175[i] = Census_daily175
#     NICU$Census_daily176[i] = Census_daily176
#     NICU$Census_daily177[i] = Census_daily177
#     NICU$Census_daily178[i] = Census_daily178
#     NICU$Census_daily179[i] = Census_daily179
#     NICU$Census_daily180[i] = Census_daily180
#     NICU$Census_daily181[i] = Census_daily181
#     NICU$Census_daily182[i] = Census_daily182
#     NICU$Census_daily183[i] = Census_daily183
#     NICU$Census_daily184[i] = Census_daily184
#     NICU$Census_daily185[i] = Census_daily185
#     NICU$Census_daily186[i] = Census_daily186
#     NICU$Census_daily187[i] = Census_daily187
#     NICU$Census_daily188[i] = Census_daily188
#     NICU$Census_daily189[i] = Census_daily189
#     NICU$Census_daily190[i] = Census_daily190
#     NICU$Census_daily191[i] = Census_daily191
#     NICU$Census_daily192[i] = Census_daily192
#     NICU$Census_daily193[i] = Census_daily193
#     NICU$Census_daily194[i] = Census_daily194
#     NICU$Census_daily195[i] = Census_daily195
#     NICU$Census_daily196[i] = Census_daily196
#     NICU$Census_daily197[i] = Census_daily197
#     NICU$Census_daily198[i] = Census_daily198
#     NICU$Census_daily199[i] = Census_daily199
#     NICU$Census_daily200[i] = Census_daily200
#     NICU$Census_daily201[i] = Census_daily201
#     NICU$Census_daily202[i] = Census_daily202
#     NICU$Census_daily203[i] = Census_daily203
#     NICU$Census_daily204[i] = Census_daily204
#     NICU$Census_daily205[i] = Census_daily205
#     NICU$Census_daily206[i] = Census_daily206
#     NICU$Census_daily207[i] = Census_daily207
#     NICU$Census_daily208[i] = Census_daily208
#     NICU$Census_daily209[i] = Census_daily209
#     NICU$Census_daily210[i] = Census_daily210
#     NICU$Census_daily211[i] = Census_daily211
#     NICU$Census_daily212[i] = Census_daily212
#     NICU$Census_daily213[i] = Census_daily213
#     NICU$Census_daily214[i] = Census_daily214
#     NICU$Census_daily215[i] = Census_daily215
#     NICU$Census_daily216[i] = Census_daily216
#     NICU$Census_daily217[i] = Census_daily217
#     NICU$Census_daily218[i] = Census_daily218
#     NICU$Census_daily219[i] = Census_daily219
#     NICU$Census_daily220[i] = Census_daily220
#     NICU$Census_daily221[i] = Census_daily221
#     NICU$Census_daily222[i] = Census_daily222
#     NICU$Census_daily223[i] = Census_daily223
#     NICU$Census_daily224[i] = Census_daily224
#     NICU$Census_daily225[i] = Census_daily225
#     NICU$Census_daily226[i] = Census_daily226
#     NICU$Census_daily227[i] = Census_daily227
#     NICU$Census_daily228[i] = Census_daily228
#     NICU$Census_daily229[i] = Census_daily229
#     NICU$Census_daily230[i] = Census_daily230
#     NICU$Census_daily231[i] = Census_daily231
#     NICU$Census_daily232[i] = Census_daily232
#     NICU$Census_daily233[i] = Census_daily233
#     NICU$Census_daily234[i] = Census_daily234
#     NICU$Census_daily235[i] = Census_daily235
#     NICU$Census_daily236[i] = Census_daily236
#     NICU$Census_daily237[i] = Census_daily237
#     NICU$Census_daily238[i] = Census_daily238
#     NICU$Census_daily239[i] = Census_daily239
#     NICU$Census_daily240[i] = Census_daily240
#     NICU$Census_daily241[i] = Census_daily241
#     NICU$Census_daily242[i] = Census_daily242
#     NICU$Census_daily243[i] = Census_daily243
#     NICU$Census_daily244[i] = Census_daily244
#     NICU$Census_daily245[i] = Census_daily245
#     NICU$Census_daily246[i] = Census_daily246
#     NICU$Census_daily247[i] = Census_daily247
#     NICU$Census_daily248[i] = Census_daily248
#     NICU$Census_daily249[i] = Census_daily249
#     NICU$Census_daily250[i] = Census_daily250
#     NICU$Census_daily251[i] = Census_daily251
#     NICU$Census_daily252[i] = Census_daily252
#     NICU$Census_daily253[i] = Census_daily253
#     NICU$Census_daily254[i] = Census_daily254
#     NICU$Census_daily255[i] = Census_daily255
#     NICU$Census_daily256[i] = Census_daily256
#     NICU$Census_daily257[i] = Census_daily257
#     NICU$Census_daily258[i] = Census_daily258
#     NICU$Census_daily259[i] = Census_daily259
#     NICU$Census_daily260[i] = Census_daily260
#     NICU$Census_daily261[i] = Census_daily261
#     NICU$Census_daily262[i] = Census_daily262
#     NICU$Census_daily263[i] = Census_daily263
#     NICU$Census_daily264[i] = Census_daily264
#     NICU$Census_daily265[i] = Census_daily265
#     NICU$Census_daily266[i] = Census_daily266
#     NICU$Census_daily267[i] = Census_daily267
#     NICU$Census_daily268[i] = Census_daily268
#     NICU$Census_daily269[i] = Census_daily269
#     NICU$Census_daily270[i] = Census_daily270
#     NICU$Census_daily271[i] = Census_daily271
#     NICU$Census_daily272[i] = Census_daily272
#     NICU$Census_daily273[i] = Census_daily273
#     NICU$Census_daily274[i] = Census_daily274
#     NICU$Census_daily275[i] = Census_daily275
#     NICU$Census_daily276[i] = Census_daily276
#     NICU$Census_daily277[i] = Census_daily277
#     NICU$Census_daily278[i] = Census_daily278
#     NICU$Census_daily279[i] = Census_daily279
#     NICU$Census_daily280[i] = Census_daily280
#     NICU$Census_daily281[i] = Census_daily281
#     NICU$Census_daily282[i] = Census_daily282
#     NICU$Census_daily283[i] = Census_daily283
#     NICU$Census_daily284[i] = Census_daily284
#     NICU$Census_daily285[i] = Census_daily285
#     NICU$Census_daily286[i] = Census_daily286
#     NICU$Census_daily287[i] = Census_daily287
#     NICU$Census_daily288[i] = Census_daily288
#     NICU$Census_daily289[i] = Census_daily289
#     NICU$Census_daily290[i] = Census_daily290
#     NICU$Census_daily291[i] = Census_daily291
#     NICU$Census_daily292[i] = Census_daily292
#     NICU$Census_daily293[i] = Census_daily293
#     NICU$Census_daily294[i] = Census_daily294
#     NICU$Census_daily295[i] = Census_daily295
#     NICU$Census_daily296[i] = Census_daily296
#     NICU$Census_daily297[i] = Census_daily297
#     NICU$Census_daily298[i] = Census_daily298
#     NICU$Census_daily299[i] = Census_daily299
#     NICU$Census_daily300[i] = Census_daily300
#     NICU$Census_daily301[i] = Census_daily301
#     NICU$Census_daily302[i] = Census_daily302
#     NICU$Census_daily303[i] = Census_daily303
#     NICU$Census_daily304[i] = Census_daily304
#     NICU$Census_daily305[i] = Census_daily305
#     NICU$Census_daily306[i] = Census_daily306
#     NICU$Census_daily307[i] = Census_daily307
#     NICU$Census_daily308[i] = Census_daily308
#     NICU$Census_daily309[i] = Census_daily309
#     NICU$Census_daily310[i] = Census_daily310
#     NICU$Census_daily311[i] = Census_daily311
#     NICU$Census_daily312[i] = Census_daily312
#     NICU$Census_daily313[i] = Census_daily313
#     NICU$Census_daily314[i] = Census_daily314
#     NICU$Census_daily315[i] = Census_daily315
#     NICU$Census_daily316[i] = Census_daily316
#     NICU$Census_daily317[i] = Census_daily317
#     NICU$Census_daily318[i] = Census_daily318
#     NICU$Census_daily319[i] = Census_daily319
#     NICU$Census_daily320[i] = Census_daily320
#     NICU$Census_daily321[i] = Census_daily321
#     NICU$Census_daily322[i] = Census_daily322
#     NICU$Census_daily323[i] = Census_daily323
#     NICU$Census_daily324[i] = Census_daily324
#     NICU$Census_daily325[i] = Census_daily325
#     NICU$Census_daily326[i] = Census_daily326
#     NICU$Census_daily327[i] = Census_daily327
#     NICU$Census_daily328[i] = Census_daily328
#     NICU$Census_daily329[i] = Census_daily329
#     NICU$Census_daily330[i] = Census_daily330
#     NICU$Census_daily331[i] = Census_daily331
#     NICU$Census_daily332[i] = Census_daily332
#     NICU$Census_daily333[i] = Census_daily333
#     NICU$Census_daily334[i] = Census_daily334
#     NICU$Census_daily335[i] = Census_daily335
#     NICU$Census_daily336[i] = Census_daily336
#     NICU$Census_daily337[i] = Census_daily337
#     NICU$Census_daily338[i] = Census_daily338
#     NICU$Census_daily339[i] = Census_daily339
#     NICU$Census_daily340[i] = Census_daily340
#     NICU$Census_daily341[i] = Census_daily341
#     NICU$Census_daily342[i] = Census_daily342
#     NICU$Census_daily343[i] = Census_daily343
#     NICU$Census_daily344[i] = Census_daily344
#     NICU$Census_daily345[i] = Census_daily345
#     NICU$Census_daily346[i] = Census_daily346
#     NICU$Census_daily347[i] = Census_daily347
#     NICU$Census_daily348[i] = Census_daily348
#     NICU$Census_daily349[i] = Census_daily349
#     NICU$Census_daily350[i] = Census_daily350
#     NICU$Census_daily351[i] = Census_daily351
#     NICU$Census_daily352[i] = Census_daily352
#     NICU$Census_daily353[i] = Census_daily353
#     NICU$Census_daily354[i] = Census_daily354
#     NICU$Census_daily355[i] = Census_daily355
#     NICU$Census_daily356[i] = Census_daily356
#     NICU$Census_daily357[i] = Census_daily357
#     NICU$Census_daily358[i] = Census_daily358
#     NICU$Census_daily359[i] = Census_daily359
#     NICU$Census_daily360[i] = Census_daily360
#     NICU$Census_daily361[i] = Census_daily361
#     NICU$Census_daily362[i] = Census_daily362
#     NICU$Census_daily363[i] = Census_daily363
#     NICU$Census_daily364[i] = Census_daily364
#     NICU$Census_daily365[i] = Census_daily365
#     NICU$Census_daily366[i] = Census_daily366
#     NICU$Census_daily367[i] = Census_daily367
#     NICU$Census_daily368[i] = Census_daily368
#     NICU$Census_daily369[i] = Census_daily369
#     NICU$Census_daily370[i] = Census_daily370
#     NICU$Census_daily371[i] = Census_daily371
#     NICU$Census_daily372[i] = Census_daily372
#     NICU$Census_daily373[i] = Census_daily373
#     NICU$Census_daily374[i] = Census_daily374
#     NICU$Census_daily375[i] = Census_daily375
#     NICU$Census_daily376[i] = Census_daily376
#     NICU$Census_daily377[i] = Census_daily377
#     NICU$Census_daily378[i] = Census_daily378
#     NICU$Census_daily379[i] = Census_daily379
#     NICU$Census_daily380[i] = Census_daily380
#     NICU$Census_daily381[i] = Census_daily381
#     NICU$Census_daily382[i] = Census_daily382
#     NICU$Census_daily383[i] = Census_daily383
#     NICU$Census_daily384[i] = Census_daily384
#     NICU$Census_daily385[i] = Census_daily385
#     NICU$Census_daily386[i] = Census_daily386
#     NICU$Census_daily387[i] = Census_daily387
#     NICU$Census_daily388[i] = Census_daily388
#     NICU$Census_daily389[i] = Census_daily389
#     NICU$Census_daily390[i] = Census_daily390
#     NICU$Census_daily391[i] = Census_daily391
#     NICU$Census_daily392[i] = Census_daily392
#     NICU$Census_daily393[i] = Census_daily393
#     NICU$Census_daily394[i] = Census_daily394
#     NICU$Census_daily395[i] = Census_daily395
#     NICU$Census_daily396[i] = Census_daily396
#     NICU$Census_daily397[i] = Census_daily397
#     NICU$Census_daily398[i] = Census_daily398
#     NICU$Census_daily399[i] = Census_daily399
#     NICU$Census_daily400[i] = Census_daily400
#     NICU$Census_daily401[i] = Census_daily401
#     NICU$Census_daily402[i] = Census_daily402
#     NICU$Census_daily403[i] = Census_daily403
#     NICU$Census_daily404[i] = Census_daily404
#     NICU$Census_daily405[i] = Census_daily405
#     NICU$Census_daily406[i] = Census_daily406
#     NICU$Census_daily407[i] = Census_daily407
#     NICU$Census_daily408[i] = Census_daily408
  }
}
rm(i,j,count_45,count_50,count_55,count_60,count_65,count_census,daily_census,followup,count_preterm25,count_preterm28,count_preterm32,daily_preterms25,daily_preterms28,daily_preterms32)
#rm(Census_daily1,Census_daily2,Census_daily3,Census_daily4,Census_daily5,Census_daily6,Census_daily7,Census_daily8,Census_daily9,Census_daily10,Census_daily11,Census_daily12,Census_daily13,Census_daily14,Census_daily15,Census_daily16,Census_daily17,Census_daily18,Census_daily19,Census_daily20,Census_daily21,Census_daily22,Census_daily23,Census_daily24,Census_daily25,Census_daily26,Census_daily27,Census_daily28,Census_daily29,Census_daily30,Census_daily31,Census_daily32,Census_daily33,Census_daily34,Census_daily35,Census_daily36,Census_daily37,Census_daily38,Census_daily39,Census_daily40,Census_daily41,Census_daily42,Census_daily43,Census_daily44,Census_daily45,Census_daily46,Census_daily47,Census_daily48,Census_daily49,Census_daily50,Census_daily51,Census_daily52,Census_daily53,Census_daily54,Census_daily55,Census_daily56,Census_daily57,Census_daily58,Census_daily59,Census_daily60,Census_daily61,Census_daily62,Census_daily63,Census_daily64,Census_daily65,Census_daily66,Census_daily67,Census_daily68,Census_daily69,Census_daily70,Census_daily71,Census_daily72,Census_daily73,Census_daily74,Census_daily75,Census_daily76,Census_daily77,Census_daily78,Census_daily79,Census_daily80,Census_daily81,Census_daily82,Census_daily83,Census_daily84,Census_daily85,Census_daily86,Census_daily87,Census_daily88,Census_daily89,Census_daily90,Census_daily91,Census_daily92,Census_daily93,Census_daily94,Census_daily95,Census_daily96,Census_daily97,Census_daily98,Census_daily99,Census_daily100,Census_daily101,Census_daily102,Census_daily103,Census_daily104,Census_daily105,Census_daily106,Census_daily107,Census_daily108,Census_daily109,Census_daily110,Census_daily111,Census_daily112,Census_daily113,Census_daily114,Census_daily115,Census_daily116,Census_daily117,Census_daily118,Census_daily119,Census_daily120,Census_daily121,Census_daily122,Census_daily123,Census_daily124,Census_daily125,Census_daily126,Census_daily127,Census_daily128,Census_daily129,Census_daily130,Census_daily131,Census_daily132,Census_daily133,Census_daily134,Census_daily135,Census_daily136,Census_daily137,Census_daily138,Census_daily139,Census_daily140,Census_daily141,Census_daily142,Census_daily143,Census_daily144,Census_daily145,Census_daily146,Census_daily147,Census_daily148,Census_daily149,Census_daily150,Census_daily151,Census_daily152,Census_daily153,Census_daily154,Census_daily155,Census_daily156,Census_daily157,Census_daily158,Census_daily159,Census_daily160,Census_daily161,Census_daily162,Census_daily163,Census_daily164,Census_daily165,Census_daily166,Census_daily167,Census_daily168,Census_daily169,Census_daily170,Census_daily171,Census_daily172,Census_daily173,Census_daily174,Census_daily175,Census_daily176,Census_daily177,Census_daily178,Census_daily179,Census_daily180,Census_daily181,Census_daily182,Census_daily183,Census_daily184,Census_daily185,Census_daily186,Census_daily187,Census_daily188,Census_daily189,Census_daily190,Census_daily191,Census_daily192,Census_daily193,Census_daily194,Census_daily195,Census_daily196,Census_daily197,Census_daily198,Census_daily199,
#   Census_daily200,Census_daily201,Census_daily202,Census_daily203,Census_daily204,Census_daily205,Census_daily206,Census_daily207,Census_daily208,Census_daily209,Census_daily210,Census_daily211,Census_daily212,Census_daily213,Census_daily214,Census_daily215,Census_daily216,Census_daily217,Census_daily218,Census_daily219,Census_daily220,Census_daily221,Census_daily222,Census_daily223,Census_daily224,Census_daily225,Census_daily226,Census_daily227,Census_daily228,Census_daily229,Census_daily230,Census_daily231,Census_daily232,Census_daily233,Census_daily234,Census_daily235,Census_daily236,Census_daily237,Census_daily238,Census_daily239,Census_daily240,Census_daily241,Census_daily242,Census_daily243,Census_daily244,Census_daily245,Census_daily246,Census_daily247,Census_daily248,Census_daily249,Census_daily250,Census_daily251,Census_daily252,Census_daily253,Census_daily254,Census_daily255,Census_daily256,Census_daily257,Census_daily258,Census_daily259,Census_daily260,Census_daily261,Census_daily262,Census_daily263,Census_daily264,Census_daily265,Census_daily266,Census_daily267,Census_daily268,Census_daily269,Census_daily270,Census_daily271,Census_daily272,Census_daily273,Census_daily274,Census_daily275,Census_daily276,Census_daily277,Census_daily278,Census_daily279,Census_daily280,Census_daily281,Census_daily282,Census_daily283,Census_daily284,Census_daily285,Census_daily286,Census_daily287,Census_daily288,Census_daily289,Census_daily290,Census_daily291,Census_daily292,Census_daily293,Census_daily294,Census_daily295,Census_daily296,Census_daily297,Census_daily298,Census_daily299,Census_daily300,Census_daily301,Census_daily302,Census_daily303,Census_daily304,Census_daily305,Census_daily306,Census_daily307,Census_daily308,Census_daily309,Census_daily310,Census_daily311,Census_daily312,Census_daily313,Census_daily314,Census_daily315,Census_daily316,Census_daily317,Census_daily318,Census_daily319,Census_daily320,Census_daily321,Census_daily322,Census_daily323,Census_daily324,Census_daily325,Census_daily326,Census_daily327,Census_daily328,Census_daily329,Census_daily330,Census_daily331,Census_daily332,Census_daily333,Census_daily334,Census_daily335,Census_daily336,Census_daily337,Census_daily338,Census_daily339,Census_daily340,Census_daily341,Census_daily342,Census_daily343,Census_daily344,Census_daily345,Census_daily346,Census_daily347,Census_daily348,Census_daily349,Census_daily350,Census_daily351,Census_daily352,Census_daily353,Census_daily354,Census_daily355,Census_daily356,Census_daily357,Census_daily358,Census_daily359,Census_daily360,Census_daily361,Census_daily362,Census_daily363,Census_daily364,Census_daily365,Census_daily366,Census_daily367,Census_daily368,Census_daily369,Census_daily370,Census_daily371,Census_daily372,Census_daily373,Census_daily374,Census_daily375,Census_daily376,Census_daily377,Census_daily378,Census_daily379,Census_daily380,Census_daily381,Census_daily382,Census_daily383,Census_daily384,Census_daily385,Census_daily386,Census_daily387,Census_daily388,Census_daily389,Census_daily390,Census_daily391,
#   Census_daily392,Census_daily393,Census_daily394,Census_daily395,Census_daily396,Census_daily397,Census_daily398,Census_daily399,Census_daily400,Census_daily401,Census_daily402,Census_daily403,Census_daily404,Census_daily405,Census_daily406,Census_daily407,Census_daily408)

#force season as factor to change reference groups to Winter
NICU$Admission_season = factor(NICU$Admission_season, levels=c("Winter","Spring","Summer","Fall"))

#census threshhold indicators, from RL: High census from a staffing perspective is above 55. When above 65, even if we have bed capacity into low 70’s, we usually bump kids to other units if possible.
NICU$Census_45 = ifelse(NICU$Census_45_n>=1, 1, 0)
NICU$Census_50 = ifelse(NICU$Census_50_n>=1, 1, 0)
NICU$Census_55 = ifelse(NICU$Census_55_n>=1, 1, 0)
NICU$Census_60 = ifelse(NICU$Census_60_n>=1, 1, 0)
NICU$Census_65 = ifelse(NICU$Census_65_n>=1, 1, 0)

#cumulative percent of very preterm infants in unit
NICU$Census_preterm25_percent = round(NICU$Census_preterm25_cumulative / NICU$Census_cumulative, digits=2)*100
NICU$Census_preterm28_percent = round(NICU$Census_preterm28_cumulative / NICU$Census_cumulative, digits=2)*100
NICU$Census_preterm32_percent = round(NICU$Census_preterm32_cumulative / NICU$Census_cumulative, digits=2)*100

#average census during at-risk period for the cohort
NICU$Census_average_cohort = round(NICU$Census_cumulative / NICU$Cohort_time)

#NICU$Census_55 = ifelse(NICU$Census_average>=55, 1, 0)
#NICU$Census_65 = ifelse(NICU$Census_average>=65, 1, 0)


### SUBSET and CREATE COHORT ###

#no birthweight means not admitted to NICU
NICU = NICU[!is.na(NICU$Birthweight), ]

#limit data to 1997 to 2014, 1997 is the current Christiana Hospital NICU design
NICU = NICU[NICU$Admission_year>=1997 & NICU$Admission_year<=2014, ]

#limit to single admissions only
NICU = NICU[NICU$Admission_n==1, ]

#drop all unnecessary census vars
# NICU$Census_daily194 = NULL
# NICU$Census_daily195 = NULL
# NICU$Census_daily196 = NULL
# NICU$Census_daily197 = NULL
# NICU$Census_daily198 = NULL
# NICU$Census_daily199 = NULL
# NICU$Census_daily200 = NULL
# NICU$Census_daily201 = NULL
# NICU$Census_daily202 = NULL
# NICU$Census_daily203 = NULL
# NICU$Census_daily204 = NULL
# NICU$Census_daily205 = NULL
# NICU$Census_daily206 = NULL
# NICU$Census_daily207 = NULL
# NICU$Census_daily208 = NULL
# NICU$Census_daily209 = NULL
# NICU$Census_daily210 = NULL
# NICU$Census_daily211 = NULL
# NICU$Census_daily212 = NULL
# NICU$Census_daily213 = NULL
# NICU$Census_daily214 = NULL
# NICU$Census_daily215 = NULL
# NICU$Census_daily216 = NULL
# NICU$Census_daily217 = NULL
# NICU$Census_daily218 = NULL
# NICU$Census_daily219 = NULL
# NICU$Census_daily220 = NULL
# NICU$Census_daily221 = NULL
# NICU$Census_daily222 = NULL
# NICU$Census_daily223 = NULL
# NICU$Census_daily224 = NULL
# NICU$Census_daily225 = NULL
# NICU$Census_daily226 = NULL
# NICU$Census_daily227 = NULL
# NICU$Census_daily228 = NULL
# NICU$Census_daily229 = NULL
# NICU$Census_daily230 = NULL
# NICU$Census_daily231 = NULL
# NICU$Census_daily232 = NULL
# NICU$Census_daily233 = NULL
# NICU$Census_daily234 = NULL
# NICU$Census_daily235 = NULL
# NICU$Census_daily236 = NULL
# NICU$Census_daily237 = NULL
# NICU$Census_daily238 = NULL
# NICU$Census_daily239 = NULL
# NICU$Census_daily240 = NULL
# NICU$Census_daily241 = NULL
# NICU$Census_daily242 = NULL
# NICU$Census_daily243 = NULL
# NICU$Census_daily244 = NULL
# NICU$Census_daily245 = NULL
# NICU$Census_daily246 = NULL
# NICU$Census_daily247 = NULL
# NICU$Census_daily248 = NULL
# NICU$Census_daily249 = NULL
# NICU$Census_daily250 = NULL
# NICU$Census_daily251 = NULL
# NICU$Census_daily252 = NULL
# NICU$Census_daily253 = NULL
# NICU$Census_daily254 = NULL
# NICU$Census_daily255 = NULL
# NICU$Census_daily256 = NULL
# NICU$Census_daily257 = NULL
# NICU$Census_daily258 = NULL
# NICU$Census_daily259 = NULL
# NICU$Census_daily260 = NULL
# NICU$Census_daily261 = NULL
# NICU$Census_daily262 = NULL
# NICU$Census_daily263 = NULL
# NICU$Census_daily264 = NULL
# NICU$Census_daily265 = NULL
# NICU$Census_daily266 = NULL
# NICU$Census_daily267 = NULL
# NICU$Census_daily268 = NULL
# NICU$Census_daily269 = NULL
# NICU$Census_daily270 = NULL
# NICU$Census_daily271 = NULL
# NICU$Census_daily272 = NULL
# NICU$Census_daily273 = NULL
# NICU$Census_daily274 = NULL
# NICU$Census_daily275 = NULL
# NICU$Census_daily276 = NULL
# NICU$Census_daily277 = NULL
# NICU$Census_daily278 = NULL
# NICU$Census_daily279 = NULL
# NICU$Census_daily280 = NULL
# NICU$Census_daily281 = NULL
# NICU$Census_daily282 = NULL
# NICU$Census_daily283 = NULL
# NICU$Census_daily284 = NULL
# NICU$Census_daily285 = NULL
# NICU$Census_daily286 = NULL
# NICU$Census_daily287 = NULL
# NICU$Census_daily288 = NULL
# NICU$Census_daily289 = NULL
# NICU$Census_daily290 = NULL
# NICU$Census_daily291 = NULL
# NICU$Census_daily292 = NULL
# NICU$Census_daily293 = NULL
# NICU$Census_daily294 = NULL
# NICU$Census_daily295 = NULL
# NICU$Census_daily296 = NULL
# NICU$Census_daily297 = NULL
# NICU$Census_daily298 = NULL
# NICU$Census_daily299 = NULL
# NICU$Census_daily300 = NULL
# NICU$Census_daily301 = NULL
# NICU$Census_daily302 = NULL
# NICU$Census_daily303 = NULL
# NICU$Census_daily304 = NULL
# NICU$Census_daily305 = NULL
# NICU$Census_daily306 = NULL
# NICU$Census_daily307 = NULL
# NICU$Census_daily308 = NULL
# NICU$Census_daily309 = NULL
# NICU$Census_daily310 = NULL
# NICU$Census_daily311 = NULL
# NICU$Census_daily312 = NULL
# NICU$Census_daily313 = NULL
# NICU$Census_daily314 = NULL
# NICU$Census_daily315 = NULL
# NICU$Census_daily316 = NULL
# NICU$Census_daily317 = NULL
# NICU$Census_daily318 = NULL
# NICU$Census_daily319 = NULL
# NICU$Census_daily320 = NULL
# NICU$Census_daily321 = NULL
# NICU$Census_daily322 = NULL
# NICU$Census_daily323 = NULL
# NICU$Census_daily324 = NULL
# NICU$Census_daily325 = NULL
# NICU$Census_daily326 = NULL
# NICU$Census_daily327 = NULL
# NICU$Census_daily328 = NULL
# NICU$Census_daily329 = NULL
# NICU$Census_daily330 = NULL
# NICU$Census_daily331 = NULL
# NICU$Census_daily332 = NULL
# NICU$Census_daily333 = NULL
# NICU$Census_daily334 = NULL
# NICU$Census_daily335 = NULL
# NICU$Census_daily336 = NULL
# NICU$Census_daily337 = NULL
# NICU$Census_daily338 = NULL
# NICU$Census_daily339 = NULL
# NICU$Census_daily340 = NULL
# NICU$Census_daily341 = NULL
# NICU$Census_daily342 = NULL
# NICU$Census_daily343 = NULL
# NICU$Census_daily344 = NULL
# NICU$Census_daily345 = NULL
# NICU$Census_daily346 = NULL
# NICU$Census_daily347 = NULL
# NICU$Census_daily348 = NULL
# NICU$Census_daily349 = NULL
# NICU$Census_daily350 = NULL
# NICU$Census_daily351 = NULL
# NICU$Census_daily352 = NULL
# NICU$Census_daily353 = NULL
# NICU$Census_daily354 = NULL
# NICU$Census_daily355 = NULL
# NICU$Census_daily356 = NULL
# NICU$Census_daily357 = NULL
# NICU$Census_daily358 = NULL
# NICU$Census_daily359 = NULL
# NICU$Census_daily360 = NULL
# NICU$Census_daily361 = NULL
# NICU$Census_daily362 = NULL
# NICU$Census_daily363 = NULL
# NICU$Census_daily364 = NULL
# NICU$Census_daily365 = NULL
# NICU$Census_daily366 = NULL
# NICU$Census_daily367 = NULL
# NICU$Census_daily368 = NULL
# NICU$Census_daily369 = NULL
# NICU$Census_daily370 = NULL
# NICU$Census_daily371 = NULL
# NICU$Census_daily372 = NULL
# NICU$Census_daily373 = NULL
# NICU$Census_daily374 = NULL
# NICU$Census_daily375 = NULL
# NICU$Census_daily376 = NULL
# NICU$Census_daily377 = NULL
# NICU$Census_daily378 = NULL
# NICU$Census_daily379 = NULL
# NICU$Census_daily380 = NULL
# NICU$Census_daily381 = NULL
# NICU$Census_daily382 = NULL
# NICU$Census_daily383 = NULL
# NICU$Census_daily384 = NULL
# NICU$Census_daily385 = NULL
# NICU$Census_daily386 = NULL
# NICU$Census_daily387 = NULL
# NICU$Census_daily388 = NULL
# NICU$Census_daily389 = NULL
# NICU$Census_daily390 = NULL
# NICU$Census_daily391 = NULL
# NICU$Census_daily392 = NULL
# NICU$Census_daily393 = NULL
# NICU$Census_daily394 = NULL
# NICU$Census_daily395 = NULL
# NICU$Census_daily396 = NULL
# NICU$Census_daily397 = NULL
# NICU$Census_daily398 = NULL
# NICU$Census_daily399 = NULL
# NICU$Census_daily400 = NULL
# NICU$Census_daily401 = NULL
# NICU$Census_daily402 = NULL
# NICU$Census_daily403 = NULL
# NICU$Census_daily404 = NULL
# NICU$Census_daily405 = NULL
# NICU$Census_daily406 = NULL
# NICU$Census_daily407 = NULL
# NICU$Census_daily408 = NULL

# #assemble census cohort, sampled on exposure: census_55
# census55 = NICU[NICU$Census_55==1, ]
# 
# #randomly select unexposed 1:1
# set.seed(27)
# census55 = rbind(census55, NICU[sample(which(NICU$Census_55==0), nrow(census55), replace=F), ])
# row.names(census55) = NULL

#load("/Users/ndgoldstein/Documents/CCHS/NICU/Sepsis/NICU_sepsis.RData")


### DESCRIPTIVES for ENTIRE COHORT ###

CrossTable(NICU$Sepsis_late)
describe(NICU$Census_preterm32_percent)
describe(NICU$Census_average_cohort)
describeBy(NICU$Census_preterm32_percent, NICU$Sepsis_late)
describeBy(NICU$Census_average_cohort, NICU$Sepsis_late)

describe(NICU$Mom_age)
CrossTable(NICU$Race_ethnicity)
CrossTable(NICU$Delivery_method)
describe(NICU$LOS); quantile(NICU$LOS,na.rm=T)

describeBy(NICU$Mom_age, NICU$Sepsis_late); t.test(NICU$Mom_age ~ NICU$Sepsis_late)
CrossTable(NICU$Race_ethnicity, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Delivery_method, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Chorio_clinical, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$PROM, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(NICU$Census_average_cohort, NICU$Sepsis_late); t.test(NICU$Mom_age ~ NICU$Sepsis_late)
describeBy(NICU$Census_preterm32_percent, NICU$Sepsis_late); t.test(NICU$Mom_age ~ NICU$Sepsis_late)
CrossTable(NICU$Preterm, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(NICU$LOS, NICU$Sepsis_late); t.test(NICU$LOS_log ~ NICU$Sepsis_late)
CrossTable(NICU$Vent, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
#CrossTable(NICU$Sepsis_stay, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Admission_season, NICU$Sepsis_late, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(NICU$Admission_year, NICU$Sepsis_late); t.test(NICU$Admission_year ~ NICU$Sepsis_late)


### DISTRIBUTION of late-onset SEPSIS ###

#plot distribution of late onset sepsis
barplot(tapply(NICU$Sepsis_late, INDEX=NICU$Sepsis_onset_year, FUN=sum), ylim=c(0,60), main="Distribution of cases of sepsis in cohort", xlab="Year", ylab="Cases of sepsis")

#plot distribution of sepsis by month
barplot(tapply(NICU$Sepsis_late, INDEX=NICU$Sepsis_onset_month, FUN=sum), ylim=c(0,70), main="Distribution of cases of sepsis in cohort", xlab="Month", ylab="Cases of sepsis")

#plot distribution of sepsis by season
barplot(tapply(NICU$Sepsis_late, INDEX=NICU$Sepsis_onset_season, FUN=sum), ylim=c(0,200), main="Distribution of cases of sepsis in cohort", xlab="Season", ylab="Cases of sepsis")

#global test for difference by season
prop.test(tapply(NICU$Sepsis_late, INDEX=NICU$Sepsis_onset_season, FUN=sum), n=rep(sum(NICU$Sepsis_late,na.rm=T),4))

#plot distribution of sepsis causes by year & season
midpoints = barplot(tapply(NICU$Sepsis_late, INDEX=NICU$Sepsis_onset_seasonyear, FUN=sum), ylim=c(0,22), main="Confirmed cases of late-onset sepsis:\nseasonal variation 1997-2014", xlab="Year", ylab="Cases of sepsis", xaxt="n")
axis(1, at=midpoints[seq(1,72,4)], labels=seq(1997,2014), las=2, cex.axis=0.8)

#fit a time-series object
sepsis_ts = ts(as.numeric(tapply(NICU$Sepsis_late, INDEX=NICU$Sepsis_onset_seasonyear, FUN=sum)), frequency=4, start=1997)

#plot the time series
plot.ts(sepsis_ts, main="Seasonality of Disease: 1997-2014", ylab="Outcome")

#decomposition plot; using multiplicative since season effect varies over time
plot(decompose(sepsis_ts, type="multiplicative"))

#check significance of trend and season
#summary(lm(sepsis_ts ~ 0+time(sepsis_ts)+factor(cycle(sepsis_ts)), na.action=NULL))


### DISTRIBUTION of late-onset SEPSIS ORGANISMS ###

CrossTable(NICU$Sepsis_cause[NICU$Sepsis_late==1])
CrossTable(NICU$Sepsis_MRSA[NICU$Sepsis_late==1])

#output to file
tiff("Figure1.tif",height=4,width=6,units='in',res=1200)

#plot distribution of sepsis organisms by year
#patterns: http://r.789695.n4.nabble.com/How-to-fill-bar-plot-with-textile-rather-than-color-td847003.html
organisms = table(NICU$Sepsis_cause[NICU$Sepsis_late==1], NICU$Sepsis_onset_year[NICU$Sepsis_late==1])
n = dim(organisms)[1]
par(xpd=TRUE)
barcolors=rev(brewer.pal(7,"Greys")); barcolors[7] = "#FFFFFF"
#barcolors=rev(brewer.pal(7,"YlOrRd"))
#barcolors[2] = barcolors[1]
#barcolors[1] = "tan4"
barplot(organisms, col=barcolors, ylim=c(0,50), xlab="Year", ylab="Cases of sepsis")
legend(x=15, y=65, legend=sub("[0-9] ", "", rownames(organisms)[n:1]), fill=rev(barcolors))
rm(barcolors)

dev.off()

#plot distribution of S. aureus by year
organisms = table(NICU$Sepsis_MRSA[NICU$Sepsis_late==1], NICU$Sepsis_onset_year[NICU$Sepsis_late==1])
n = dim(organisms)[1]
par(xpd=TRUE)
barplot(organisms, col=gray(0:(n-1)/(n-1)), ylim=c(0,15), main="Methicillin profile of late-onset sepsis\n due to S. aureus infection", xlab="Year", ylab="Cases of S. aureus")
legend(x=16, y=12, legend=sub("[0-9] ", "", rownames(organisms)[n:1]), fill=gray(0:(n-1)/(n-1))[n:1])


### INCIDENCE of LATE ONSET SEPSIS ###

#all cause IR (missing cohort time means loss to follow up) per 1000 NICU days
outcome = sum(NICU$Sepsis_late[!is.na(NICU$Cohort_time)]) #numerator: cases
persontime = sum(NICU$Cohort_time[!is.na(NICU$Cohort_time)]) #denominator: cohort time
outcome/persontime*1000

#all cause IR, start of cohort
outcome = sum(NICU$Sepsis_late[!is.na(NICU$Cohort_time) & NICU$Admission_year==1997]) #numerator: cases
persontime = sum(NICU$Cohort_time[!is.na(NICU$Cohort_time) & NICU$Admission_year==1997]) #denominator: cohort time
outcome/persontime*1000

#all cause IR, end of cohort
outcome = sum(NICU$Sepsis_late[!is.na(NICU$Cohort_time) & NICU$Admission_year==2014]) #numerator: cases
persontime = sum(NICU$Cohort_time[!is.na(NICU$Cohort_time) & NICU$Admission_year==2014]) #denominator: cohort time
outcome/persontime*1000

#exposed (missing cohort time means loss to follow up) per 1000 NICU days
outcome_exposed = sum(NICU$Sepsis_late[NICU$Census_55==1 & !is.na(NICU$Cohort_time)]) #numerator: cases
persontime_exposed = sum(NICU$Cohort_time[NICU$Census_55==1 & !is.na(NICU$Cohort_time)]) #denominator: cohort time
outcome_exposed/persontime_exposed*1000

#unexposed (missing cohort time means loss to follow up) per 1000 NICU days
outcome_unexposed = sum(NICU$Sepsis_late[NICU$Census_55==0 & !is.na(NICU$Cohort_time)]) #numerator: cases
persontime_unexposed = sum(NICU$Cohort_time[NICU$Census_55==0 & !is.na(NICU$Cohort_time)]) #denominator: cohort time
outcome_unexposed/persontime_unexposed*1000

#incidence rate ratio
rateratio(c(outcome_unexposed,outcome_exposed,persontime_unexposed,persontime_exposed))

# #exposed
# outcome_exposed = sum(NICU$Sepsis_late[NICU$Census_55==1], na.rm=T) #numerator: cases
# no_outcome_exposed = sum(!NICU$Sepsis_late[NICU$Census_55==1], na.rm=T) #numerator: not cases
# 
# #unexposed
# outcome_unexposed = sum(NICU$Sepsis_late[NICU$Census_55==0], na.rm=T) #numerator: cases
# no_outcome_unexposed = sum(!NICU$Sepsis_late[NICU$Census_55==0], na.rm=T) #numerator: not cases
# 
# #cumulative incidence relative risk
# riskratio(c(no_outcome_unexposed,outcome_unexposed,no_outcome_exposed,outcome_exposed))


### ADJUSTED MODELS ###

#potential confounders
t.test(NICU$Mom_age ~ NICU$Sepsis_late)
chisq.test(NICU$Race_ethnicity, NICU$Sepsis_late)
chisq.test(NICU$Delivery_method, NICU$Sepsis_late)
chisq.test(NICU$Chorio_clinical, NICU$Sepsis_late)
chisq.test(NICU$PROM, NICU$Sepsis_late)
chisq.test(NICU$Preterm, NICU$Sepsis_late)
t.test(NICU$LOS_log ~ NICU$Sepsis_late)
chisq.test(NICU$Vent, NICU$Sepsis_late)
chisq.test(NICU$Sepsis_stay, NICU$Sepsis_late)
chisq.test(NICU$Admission_season, NICU$Sepsis_late)
t.test(NICU$Admission_year ~ NICU$Sepsis_late)

t.test(NICU$Mom_age ~ NICU$Census_55)
chisq.test(NICU$Race_ethnicity, NICU$Census_55)
chisq.test(NICU$Delivery_method, NICU$Census_55)
chisq.test(NICU$Chorio_clinical, NICU$Census_55)
chisq.test(NICU$PROM, NICU$Census_55)
chisq.test(NICU$Preterm, NICU$Census_55)
t.test(NICU$LOS_log ~ NICU$Census_55)
chisq.test(NICU$Vent, NICU$Census_55)
chisq.test(NICU$Sepsis_stay, NICU$Census_55)
chisq.test(NICU$Admission_season, NICU$Census_55)
t.test(NICU$Admission_year ~ NICU$Census_55)

#full model, backward selection removal p<0.10 and change in effect <10%
summary(coxph(Surv(Cohort_time, Sepsis_late)~as.factor(Census_55)+Mom_age+as.factor(Race_ethnicity)+as.factor(Delivery_method)+as.factor(Chorio_clinical)+as.factor(PROM)+as.factor(Preterm)+log(LOS)+as.factor(Vent)+as.factor(Sepsis_stay)+as.factor(Admission_season)+Admission_year, data=NICU))
summary(coxph(Surv(Cohort_time, Sepsis_late)~as.factor(Census_55)+as.factor(Race_ethnicity)+as.factor(Delivery_method)+as.factor(Chorio_clinical)+as.factor(PROM)+as.factor(Preterm)+LOS_log+as.factor(Vent)+as.factor(Sepsis_stay)+as.factor(Admission_season)+Admission_year, data=NICU))
summary(coxph(Surv(Cohort_time, Sepsis_late)~as.factor(Census_55)+as.factor(Race_ethnicity)+as.factor(Chorio_clinical)+as.factor(PROM)+as.factor(Preterm)+LOS_log+as.factor(Vent)+as.factor(Sepsis_stay)+as.factor(Admission_season)+Admission_year, data=NICU))
summary(coxph(Surv(Cohort_time, Sepsis_late)~as.factor(Census_55)+as.factor(Race_ethnicity)+as.factor(PROM)+as.factor(Preterm)+LOS_log+as.factor(Vent)+as.factor(Sepsis_stay)+as.factor(Admission_season)+Admission_year, data=NICU))
summary(coxph(Surv(Cohort_time, Sepsis_late)~as.factor(Census_55)+as.factor(Race_ethnicity)+as.factor(Preterm)+LOS_log+as.factor(Vent)+as.factor(Sepsis_stay)+as.factor(Admission_season)+Admission_year, data=NICU))
summary(coxph(Surv(Cohort_time, Sepsis_late)~as.factor(Census_55)+as.factor(Race_ethnicity)+as.factor(Preterm)+LOS_log+as.factor(Vent)+as.factor(Admission_season)+Admission_year, data=NICU))

#cox regression, no time dependent covariates
model = coxph(Surv(Cohort_time, Sepsis_late)~as.factor(Census_55)+as.factor(Race_ethnicity)+as.factor(Preterm)+LOS_log+as.factor(Vent)+as.factor(Admission_season)+Admission_year, data=NICU) 
summary(model) #summary of the model 
round(exp(coef(model)),2) #coefficient estimates: hazard ratios 
round(exp(confint(model)),2) #confidence intervals 
plot(survfit(model)) #survival plot 

#check proportional hazard assumption -- violated for each covariate p<0.05
cox.zph(model)
plot(cox.zph(model))

#cox regression, time dependent covariates
model = coxph(Surv(Cohort_time, Sepsis_late)~as.factor(Census_55)+as.factor(Race_ethnicity)+as.factor(Admission_season)+Admission_year+tt(Admission_year)+LOS+as.factor(Vent), data=NICU, tt=function(x, t, ...) x*log(t))

#parametric survival analysis using weibull distribution
#note: coefficient signs will be opposite, as modeling survival instead of event
model_exponential = survreg(Surv(Cohort_time, Sepsis_late)~as.factor(Census_55)+as.factor(Race_ethnicity)+as.factor(Preterm)+LOS+as.factor(Vent)+as.factor(Sepsis_stay)+as.factor(Admission_season)+Admission_year, data=NICU, dist="exponential")
model_weibull = survreg(Surv(Cohort_time, Sepsis_late)~as.factor(Census_55)+as.factor(Race_ethnicity)+as.factor(Preterm)+LOS+as.factor(Vent)+as.factor(Sepsis_stay)+as.factor(Admission_season)+Admission_year, data=NICU, dist="weibull")
anova(model_exponential,model_weibull)

#mixed effects cox regression
model_me = coxme(Surv(Cohort_time, Sepsis_late)~(1|Cluster)+as.factor(Census_55)+as.factor(Race_ethnicity)+as.factor(Preterm)+LOS+as.factor(Vent)+as.factor(Sepsis_stay)+as.factor(Admission_season)+Admission_year, data=NICU) 

#export analysis dataset for BI
# write.csv(NICU[,c("ID","Cohort_time","Sepsis_late","Cluster","Census_55","Race_ethnicity","Preterm", "LOS", "Vent", "Sepsis_stay", "Admission_season", "Admission_year","Census_daily1","Census_daily2","Census_daily3","Census_daily4","Census_daily5","Census_daily6","Census_daily7","Census_daily8","Census_daily9","Census_daily10","Census_daily11","Census_daily12","Census_daily13","Census_daily14","Census_daily15","Census_daily16","Census_daily17","Census_daily18","Census_daily19","Census_daily20","Census_daily21","Census_daily22","Census_daily23","Census_daily24","Census_daily25","Census_daily26","Census_daily27","Census_daily28","Census_daily29","Census_daily30","Census_daily31","Census_daily32","Census_daily33","Census_daily34","Census_daily35","Census_daily36","Census_daily37","Census_daily38","Census_daily39","Census_daily40","Census_daily41","Census_daily42","Census_daily43","Census_daily44","Census_daily45","Census_daily46","Census_daily47","Census_daily48","Census_daily49","Census_daily50","Census_daily51","Census_daily52","Census_daily53","Census_daily54","Census_daily55","Census_daily56","Census_daily57","Census_daily58","Census_daily59","Census_daily60","Census_daily61","Census_daily62","Census_daily63","Census_daily64","Census_daily65","Census_daily66","Census_daily67","Census_daily68","Census_daily69","Census_daily70","Census_daily71","Census_daily72","Census_daily73","Census_daily74","Census_daily75","Census_daily76","Census_daily77","Census_daily78","Census_daily79","Census_daily80","Census_daily81","Census_daily82","Census_daily83","Census_daily84","Census_daily85","Census_daily86","Census_daily87","Census_daily88","Census_daily89","Census_daily90","Census_daily91","Census_daily92","Census_daily93","Census_daily94","Census_daily95","Census_daily96","Census_daily97","Census_daily98","Census_daily99",
#                   "Census_daily100","Census_daily101","Census_daily102","Census_daily103","Census_daily104","Census_daily105","Census_daily106","Census_daily107","Census_daily108","Census_daily109","Census_daily110","Census_daily111","Census_daily112","Census_daily113","Census_daily114","Census_daily115","Census_daily116","Census_daily117","Census_daily118","Census_daily119","Census_daily120","Census_daily121","Census_daily122","Census_daily123","Census_daily124","Census_daily125","Census_daily126","Census_daily127","Census_daily128","Census_daily129","Census_daily130","Census_daily131","Census_daily132","Census_daily133","Census_daily134","Census_daily135","Census_daily136","Census_daily137","Census_daily138","Census_daily139","Census_daily140","Census_daily141","Census_daily142","Census_daily143","Census_daily144","Census_daily145","Census_daily146","Census_daily147","Census_daily148","Census_daily149","Census_daily150","Census_daily151","Census_daily152","Census_daily153","Census_daily154","Census_daily155","Census_daily156","Census_daily157","Census_daily158","Census_daily159","Census_daily160","Census_daily161","Census_daily162","Census_daily163","Census_daily164","Census_daily165","Census_daily166","Census_daily167","Census_daily168","Census_daily169","Census_daily170","Census_daily171","Census_daily172","Census_daily173","Census_daily174","Census_daily175","Census_daily176","Census_daily177","Census_daily178","Census_daily179","Census_daily180","Census_daily181","Census_daily182","Census_daily183","Census_daily184","Census_daily185","Census_daily186","Census_daily187","Census_daily188","Census_daily189","Census_daily190","Census_daily191","Census_daily192","Census_daily193")], file="Sepsis_BI.csv", na="", row.names=F)
write.csv(NICU[,c("ID","Cohort_time","Sepsis_late","Census_preterm32_percent","Race_ethnicity","Preterm","Birthweight_1500","LOS", "Vent", "Sepsis_stay", "Admission_season", "Admission_year")], file="Sepsis_BI.csv", na="", row.names=F)


### FINAL MODEL for PAPER ###

model = coxph(Surv(Cohort_time, Sepsis_late)~Census_average_cohort+Census_preterm32_percent+as.factor(Race_ethnicity)+as.factor(Preterm)+as.factor(Vent)+as.factor(Admission_season)+Admission_year, data=NICU) 
model_interaction_census = coxph(Surv(Cohort_time, Sepsis_late)~scale(Census_average_cohort)*scale(Census_preterm32_percent)+as.factor(Race_ethnicity)+as.factor(Preterm)+as.factor(Vent)+as.factor(Admission_season)+Admission_year, data=NICU) 
model_interaction_vlbw = coxph(Surv(Cohort_time, Sepsis_late)~scale(Census_preterm32_percent)*as.factor(Birthweight_1500)+as.factor(Race_ethnicity)+as.factor(Preterm)+as.factor(Vent)+as.factor(Admission_season)+Admission_year, data=NICU) 

summary(model) #summary of the model 
round(exp(coef(model)),2) #coefficient estimates: hazard ratios 
round(exp(confint(model)),2) #confidence intervals 

#sensitivity analysis
model = coxph(Surv(Cohort_time, Sepsis_late)~Census_average_cohort+Census_preterm28_percent+as.factor(Race_ethnicity)+as.factor(Preterm)+as.factor(Vent)+as.factor(Admission_season)+Admission_year, data=NICU) 
model = coxph(Surv(Cohort_time, Sepsis_late)~Census_average_cohort+Census_preterm25_percent+as.factor(Race_ethnicity)+as.factor(Preterm)+as.factor(Vent)+as.factor(Admission_season)+Admission_year, data=NICU) 
summary(model) #summary of the model 


### PLOTS ###

#estimated survival from BI's SAS output
SAS_survival = data.frame(Cohort_time=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,45,46,48,53,55,56,57,62,63,65,66,67,82,87,91,99,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,45,46,48,53,55,56,57,62,63,65,66,67,82,87,91,99), Survival=c(1,0.999943896,0.999936192,0.99992777,0.999882655,0.999834079,0.999689631,0.999462775,0.999115698,0.998816456,0.998335897,0.997810601,0.997310699,0.996955494,0.996686008,0.996161252,0.995827468,0.995514898,0.995285429,0.995067428,0.9947499,0.994630642,0.99443012,0.994377738,0.994296311,0.994154337,0.993948015,0.99382499,0.993601343,0.993502638,0.993366486,0.993224727,0.993114937,0.993077289,0.993000038,0.992881275,0.992798277,0.992670193,0.992626092,0.992535862,0.992442028,0.992393215,0.992291296,0.992237347,0.992119254,0.99199772,0.991932076,0.991852701,0.991683522,0.991507954,0.991415527,0.991188675,0.990955352,0.990706471,0.990577264,0.990443633,0.989972076,0.989694562,0.989355988,0.988860565,1,0.999901788,0.999888302,0.99987356,0.999794588,0.99970956,0.999456737,0.999059733,0.998452472,0.99792903,0.997088671,0.996170428,0.995296912,0.994676435,0.994205804,0.993289642,0.992707083,0.992161684,0.991761367,0.991381121,0.990827389,0.99061945,0.990269862,0.990178548,0.990036611,0.98978915,0.989429581,0.989215206,0.988825542,0.988653587,0.988416418,0.988169508,0.987978296,0.987912733,0.987778206,0.987571403,0.987426891,0.987203894,0.987127118,0.986970043,0.986806706,0.986721741,0.98654435,0.986450456,0.98624494,0.986033455,0.985919233,0.985781128,0.985486797,0.985181391,0.985020629,0.984626101,0.984220391,0.983787706,0.98356311,0.983330846,0.982511421,0.982029322,0.981441286,0.980581109), stringsAsFactors=F)

SAS_survival$Hazard = 1 - SAS_survival$Survival           
SAS_survival$Census_55 = c(rep(0,60),rep(1,60))

#output to file
tiff("Figure2.tif",height=4,width=6,units='in',res=1200)

#hazard plot
plot(x=SAS_survival$Cohort_time[SAS_survival$Census_55==0], y=SAS_survival$Hazard[SAS_survival$Census_55==0], type="s", xlab="Days in NICU", ylab="Hazard", lty=1, lwd=2, ylim=c(0,0.02))
lines(x=SAS_survival$Cohort_time[SAS_survival$Census_55==1], y=SAS_survival$Hazard[SAS_survival$Census_55==1], type="s", lty=2, lwd=2)
axis(4)

#add legend
legend(10,0.026,lty=c(1,2),c("High census","Non-high census"), horiz=T, xpd=T, lwd=2)

#add HR
text(20,0.019,"HR 0.57, 95% CI: 0.46, 0.70",cex=0.8)

dev.off()

#adjusted plot needs covariate data (see Kleinbaum pg658)
#complete_case = na.omit(NICU[,c("Cohort_time","Sepsis_late","Census_55","Race_ethnicity","Preterm","Vent","Sepsis_stay","Admission_season","Admission_year")])
pattern = data.frame(Census_55=0,Race_ethnicity=mean(NICU$Race_ethnicity,na.rm=T),Preterm=mean(NICU$Preterm,na.rm=T),Vent=mean(NICU$Vent,na.rm=T),Sepsis_stay=mean(NICU$Sepsis_stay,na.rm=T),Admission_season="Winter",Admission_year=mean(NICU$Admission_year,na.rm=T), stringsAsFactors=F)

plot(survfit(model,pattern),ylim=c(0.979,1.00),xlim=c(0,100))
