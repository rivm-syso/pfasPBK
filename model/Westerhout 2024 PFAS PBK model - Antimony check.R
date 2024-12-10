# Generic PFAS PBK model based on Westerhout et al., 2024 with oral and dermal exposure routes
# To check the antimony model output
# Parameter values are obtained from [EFSA 2020] unless specified otherwise
# Remark: the doses are given in ug, the volumes in L, so concentrations are in ug/L
# From ug/L to uM is by dividing by MW

rm(list=ls()) # to clear out the global environment

setwd("C:/Users/westerj/OneDrive - rivm.nl/PFAS/PFAS model development")

library(deSolve)
library(ggplot2)
library(tidyverse)
library(gridExtra)

#### Settings ----
TSTART <- 0
TSTOP <- 365*80 # simulation duration (days); 50 years
exposure_duration = 365*50 # Duration of exposure (d); 50 years
DT <- 1 # time step (days)
TIME <- seq(TSTART,TSTOP,by=DT)
AGE <- 0 # to use as a starting point, e.g. 20 years of age instead of birth
age_start_derm_exp <- 14 # dermal and inhalation exposure starting age (e.g. occupational, or personal care products; exposure may start at a later age than 0)
age_stop_derm_exp <- 50 # dermal and inhalation exposure stopping age (e.g. occupational, or personal care products; exposure may stop at an earlier age than 80)

Variables_df <- as.data.frame(list(TIME = TIME))
Variables_df <- Variables_df %>%
  mutate(day = TIME+1) %>%
  mutate(dayoftheweek = rep(seq(1:7),length.out = n())) %>%
  mutate(month = ceiling(day/(365/12))-(((ceiling(day/365))-1)*12)) %>%
  mutate(year = ceiling(day/(365))) %>%
  mutate(age = AGE + (TIME/(365)))

#### Physiological parameters ----
# Physiological parameters (from Brown, et al. 1997)
#fractional blood flows
QCC = 300 # 12.5*24   # Cardiac blood output (L/d/kg^0.75)
QFC = 0.052  # Fraction cardiac output going to fat
QLC = 0.069   # Fraction cardiac output going to liver, through hepatic artery
QKC = 0.175  # Fraction cardiac output going to kidney
QSkC = 0.058  # Fraction cardiac output going to skin
QGC = 0.181  # Fraction of cardiac output going to gut and in the liver via portal artery
# Not used #QfilC = 0.035  # Fraction cardiac output to the filtrate compartment (20% of kidney blood flow)

#fractional tissue volumes
VLC = 0.026  # Fraction liver volume
VFC = 0.214  # Fraction fat volume
VKC = 0.004  # Fraction kidney volume
VfilC = 0.0004  # Fraction filtrate compartment volume (10% of kidney volume) 
VGC = 0.0171  # Fraction gut volume 
VlunC = 0.007 # Fraction lung volume (ICRP 89)
VPlasC = 0.0428 # Fraction plasma volume (58% of blood) 
Htc = 0.44 # hematocrit 

# Mother-child parameters
BWbirth=3.68 # Body weight at birth in kg new add opinion 2020 

# Skin parameters 
# Dermal exposure 
# Dermexpo = 0 # 0 = NO, 1 = YES
Dermconc = 10.0 # Dermal concentration (ug/L) 
Dermal_duration = 12 # Dermal exposure duration (e.g. skin lotion = 16 h, whereas shampoo is 0.25 h)
Cdermal_day = Dermconc*Dermal_duration/24 # Correction of the daily dermal dose beased on exposure duration
Skinthickness = 0.1 # Skin thickness (cm) 

fss <- 0.005 # fraction of the total body surface area exposed. 0.005 = palm of hands (Sheridan et al., 1995, Rhodes et al., 2013)
hsurf <- 0.01 # applied layer thickness in cm, value is assumed to be 0.1 mm
hsc <- 0.0015 # Stratum corneum thickness (cm) # [Krüse 2007]
hve <- 0.0100 # Viable epidermis thickness (cm) # [Krüse 2007]
hcell <- 0.00005 # Blood vessel wall thickness (cm) # [Burton 1954] = 0.5 um capillary wall thickness

ffatsc <- 0.05 # Fraction of fat in stratum corneum [Polak 2012]
ffatve <- 0.02 # Fraction of fat in viable epidermis [Polak 2012]
ffatepi <- 0.02 # Fraction of fat in epidermis [Polak 2012]
ffatbl <- 0.007 # Fraction of fat in blood [Polak 2012]

Kpcell <- 0.93 # cm/h; Krüse model is normalized to 1 cm^2 # Permeability coefficient from arterial wall into the blood compartment [Krüse 2007]

# Physicochemical properties - PFOA
MW = 414.07
logP = 4.81
VP = 2.34 # 10^0.37 # Pa [Zhang 2021]; 1 mmHg = 133.322368 Pa; 0.15 mmHg (https://haz-map.com/Agents/6596) = 20 Pa

Kscve = ((1-ffatve)+(ffatve*(10^logP)))/((1-ffatsc)+(ffatsc*(10^logP)))
Kver = ((1-ffatbl)+(ffatbl*(10^logP)))/((1-ffatepi)+(ffatepi*(10^logP)))
Kpsc = 0.000088 # [Franko 2011]
Kpve <- 0.000088 # [Franko 2011]

# Oral exposure via food
Oralexpo = 0.000187 # µg/kg/day [EFSA 2020; page 143]

# Drinking water exposure 
Drinkrate = 13 # Drinking water rate (mL/kg/day)
Drinkconc = 0 # Drinking water concentration (ug/L or ppb) 

# Chemical-specific parameters (PFOA) 
Tmc = 144000 # 6000*24   # Maximum resorption rate 
Kt = 55   # Resorption affinity# same as monkey
kurinec = 0.0072 # 0.0003*24 # urinary elimination rate constant (/d/kg^-0.25)# estimated 
Free = 0.02  # Free fraction of PFOA in plasma# same as monkey 
# Partition coefficients from Harada, et al 2005 
PL = 2.2   # Liver/plasma partition coefficient of PFOA
PF = 0.04   # Fat/plasma partition coefficient of PFOA
PK = 1.05   # Kidney/plasma partition coefficient of PFOA
PSk = 0.1   # Skin/plasma partition coefficient of PFOA
PR = 0.12   # Rest of the body/plasma partition coefficient of PFOA
PG = 0.05   # Gut/blood plasma coefficient of PFOA  
PLun = 1.27   # Lung/plasma partition coefficient of PFOA; [Fabrega 2014]

# Free fraction of chemical in tissues 
FreeL = Free/PL                      #liver 
FreeF = Free/PF                     #fat 
FreeK = Free/PK                     #kidney 
FreeSk = Free/PSk                #skin 
FreeR = Free/PR                   #rest of tissues 
FreeG = Free/PG                    #gut 
FreeLun = Free/PLun                    #lung 

# Mother-child parameters
# PFOA
Maternal_conc = 2.0 # maternal concentration ng/mL at delivery [EFSA opinion 2020, p368]
PT= 0.74 # placenta transfer for PFOA, new add opinion 2020 
Ratio= 0.03 # milk concentration/maternal serum concentration for PFOA during breastfeeding, new add opinion 2020 
DECLINE = 0.077 # decline of PFOA in milk was 7.7% per month, new add opinion 2020 

# Oral exposure via the mother
Milkconcentration = Maternal_conc*Ratio # initial Milk concentration at birth ng/mL 
Milkconsumption = 0.8 # milkconsumption L per day 
Intakemilka = Milkconcentration*Milkconsumption# initial intake via breastfeeding first month
Intakemilkb = Intakemilka*(1-DECLINE)#intakeviabreastfeedingsecondmonth 
Intakemilkc = Intakemilkb*(1-DECLINE)#intakeviabreastfeeding3month
Intakemilkd = Intakemilkc*(1-DECLINE)#intakeviabreastfeeding4month 
Intakemilke = Intakemilkd*(1-DECLINE)#intakeviabreastfeeding5month
Intakemilkf = Intakemilke*(1-DECLINE)#intakeviabreastfeeding6month 
Intakemilkg = Intakemilkf*(1-DECLINE)#intakeviabreastfeeding7month 
Intakemilkh = Intakemilkg*(1-DECLINE)#intakeviabreastfeeding8month 
Intakemilki = Intakemilkh*(1-DECLINE)#intakeviabreastfeeding9month 
Intakemilkj = Intakemilki*(1-DECLINE)#intakeviabreastfeeding10month 
Intakemilkk = Intakemilkj*(1-DECLINE)#intakeviabreastfeeding11month 
Intakemilkl = Intakemilkk*(1-DECLINE)#intakeviabreastfeeding12month 

CONCbirth = Maternal_conc*PT # PFOA concentration at birth, new add in 2020 
APlasbirth = CONCbirth*Free*VPlasC*BWbirth # initial amount in plasma at birth 
Aartbirth = 0.39*APlasbirth # initial amount in arterial plasma at birth
Avenbirth = 0.61*APlasbirth # initial amount in venous plasma at birth
AGbirth = APlasbirth*PG # 0.99929853 # initial amount in gut at birth
ALbirth = APlasbirth*PL # 66.8360764 # initial amount in liver at birth
AFbirth = APlasbirth*PF # 10.000151 # initial amount in fat at birth
AKbirth = APlasbirth*PK # 5.88575746 # initial amount in kidney at birth
ALunbirth = APlasbirth*PLun # initial amount in lung tissue at birth
ASkbirth = APlasbirth*PSk # initial amount in skin tissue at birth
ARbirth = APlasbirth*PR # 74.9108792 # initial amount in rest of the body at birth
Afilbirth = APlasbirth*0.0000209606 # initial amount filtrate at birth
Adelaybirth = APlasbirth*42.27910714 # initial amount in storage at birth
Aurinebirth = APlasbirth*769.1347494 # initial amount in urine at birth

#### Time-dependent variables ----
Variables_df <- Variables_df %>%
  mutate(BW = BWbirth + 4.47*age - 0.093*age^2 + 0.00061*age^3) %>%
  mutate(SkinTarea = 9.1*((BW*1000)**0.666)) %>% # Total area of skin (cm^2) 
  mutate(Skinarea = fss*SkinTarea) %>% # Exposed area of skin (cm^2)
  mutate(kurine = kurinec*BW**(-0.25)) %>%
  # Oral exposure
  mutate(Oraldose = if_else(age < 0.083, Intakemilka/BW,
                            if_else(age >= 0.083 & age < 0.167, Intakemilkb/BW,
                                    if_else(age >= 0.167 & age < 0.250, Intakemilkc/BW,
                                            if_else(age >= 0.250 & age < 0.333, Intakemilkd/BW,
                                                    if_else(age >= 0.333 & age < 0.417, Intakemilke/BW,
                                                            if_else(age >= 0.417 & age < 0.500, Intakemilkf/BW,
                                                                    if_else(age >= 0.500 & age < 0.583, Intakemilkg/BW,
                                                                            if_else(age >= 0.583 & age < 0.667, Intakemilkh/BW,
                                                                                    if_else(age >= 0.667 & age < 0.750, Intakemilki/BW,
                                                                                            if_else(age >= 0.750 & age < 0.833, Intakemilkj/BW,
                                                                                                    if_else(age >= 0.833 & age < 0.917, Intakemilkk/BW,
                                                                                                            if_else(age >= 0.917 & age < 1, Intakemilkl/BW,
                                                                                                                    Oralexpo))))))))))))) %>% # (ug/day)
  mutate(Oraldose = if_else(day > exposure_duration,0,Oraldose*BW)) %>% # stop exposure after exposure_duration
  mutate(Drinkdose = (Drinkconc*Drinkrate/1000)*BW) %>%   # (ug/day) 
  mutate(Drinkdose = if_else(day > exposure_duration,0,Drinkdose)) %>% # stop exposure after exposure_duration
  
  # Blood flows
  mutate(QC = QCC*BW**0.75) %>% # Cardiac output (L/d) 
  mutate(QCP = QC*(1-Htc)) %>% # adjust for plasma flow 
  mutate(QL = QLC*QCP) %>% # Plasma flow to liver (L/d) 
  mutate(QF = QFC*QCP) %>% # Plasma flow to fat (L/d)
  mutate(QK = QKC*QCP) %>% # Plasma flow to kidney (L/d) 
  mutate(Qfil = 0.2*QK) %>% # Plasma flow to filtrate compartment (L/d)# 20% of QK 
  mutate(QG = QGC*QCP) %>% # Plasma flow to gut (L/d) 
  mutate(QSk = QSkC*QCP) %>% #*(Skinarea/SkinTarea)*Dermexpo) %>% #plasma flow to skin 
  mutate(QR = QCP - QL - QF - QK - QG - QSk) %>% # Plasma flow to rest of the body (L/d)
  mutate(Qbal = QCP - (QR+QL+QF+QK+QG+QSk)) %>% # balance check--better be 0 
  
  # Organ volumes
  mutate(VL = VLC*BW) %>% # Liver volume (L) 
  mutate(VF = VFC*BW) %>% # Fat volume (L) 
  mutate(VK = VKC*BW) %>% # Kidney volume (L) 
  mutate(Vfil = VfilC*BW) %>% # Filtrate compartment volume (L) 
  mutate(VG = VGC*BW) %>% # Gut volume (L) 
  mutate(VPlas = VPlasC*BW) %>% # Plasma volume (L) 
  mutate(Vart_plas = 0.39*VPlas) %>% # [Brown 1997]
  mutate(Vven_plas = 0.61*VPlas) %>% # [Brown 1997]
  #mutate(VSk = (SkinTarea*Skinthickness)/1000) %>% # Total skin volume (L) 
  mutate(VSk = (Skinarea*Skinthickness)/1000) %>% # Exposed skin volume (L) 
  mutate(Vlun = VlunC*BW) %>% # [ICRP 89]
  
  mutate(VR = 0.84*BW - VL - VF - VK - Vfil - VG - VPlas - VSk - Vlun) %>% # Rest of the body volume (L) 
  mutate(Vbal = (0.84*BW)-(VL+VF+VK+Vfil+VG+VPlas+VSk+Vlun)) %>% # Balance check--better be 0 
  
  mutate(Vsurf = fss*SkinTarea*hsurf/1000) %>% # volume of dermal layer on exposed skin surface (L) 1 m^2 = 10,000 cm^2 * cm = cm^3 = ml / 1000 = L
  
  # Dermal exposure
  mutate(Csurf = Cdermal_day) %>%
  
  # In case of age-dependent dermal exposure, set the value to 0 for the ages that are not included
  mutate(Csurf = if_else(age >= age_start_derm_exp & age <= age_stop_derm_exp,Csurf,0)) %>%
  mutate(Dermaldose = Csurf*Vsurf) %>%
  
  mutate(Vsc = fss*SkinTarea*hsc/1000) %>% # volume of exposed stratum corneum (L) cm^2 * cm = cm^3 = ml / 1000 = L
  mutate(Vve = fss*SkinTarea*hve/1000) %>% # volume of exposed viable epidermis (L) cm^2 * cm = cm^3 = ml / 1000 = L
  
  mutate(CLsc = 24*Kpsc*fss*SkinTarea/1000) %>% # stratum corneum clearance (L/h) cm^2 * cm/h = cm^3/h = ml/h / 1000 = L/h
  mutate(CLve = 24*Kpve*fss*SkinTarea/1000) %>% # viable epidermis clearance (L/h) cm^2 * cm/h = cm^3/h = ml/h / 1000 = L/h
  mutate(CLcell = 24*Kpcell*fss*SkinTarea/1000) %>% # clearance from arterial wall into the blood compartment (L/h) cm^2 * cm/h = cm^3/h = ml/h / 1000 = L/h
  
  mutate(Tm = Tmc*BW**0.75) #transporter maximum 

varkurine <- approxfun(Variables_df$TIME, Variables_df$kurine, rule = 2)

varQCP <- approxfun(Variables_df$TIME, Variables_df$QCP, rule = 2)
varQG <- approxfun(Variables_df$TIME, Variables_df$QG, rule = 2)
varQL <- approxfun(Variables_df$TIME, Variables_df$QL, rule = 2)
varQF <- approxfun(Variables_df$TIME, Variables_df$QF, rule = 2)
varQK <- approxfun(Variables_df$TIME, Variables_df$QK, rule = 2)
varQfil <- approxfun(Variables_df$TIME, Variables_df$Qfil, rule = 2)
varQSk <- approxfun(Variables_df$TIME, Variables_df$QSk, rule = 2)
varQR <- approxfun(Variables_df$TIME, Variables_df$QR, rule = 2)

varVPlas <- approxfun(Variables_df$TIME, Variables_df$VPlas, rule = 2)
varVart_plas <- approxfun(Variables_df$TIME, Variables_df$Vart_plas, rule = 2)
varVven_plas <- approxfun(Variables_df$TIME, Variables_df$Vven_plas, rule = 2)
varVG <- approxfun(Variables_df$TIME, Variables_df$VG, rule = 2)
varVL <- approxfun(Variables_df$TIME, Variables_df$VL, rule = 2)
varVF <- approxfun(Variables_df$TIME, Variables_df$VF, rule = 2)
varVK <- approxfun(Variables_df$TIME, Variables_df$VK, rule = 2)
varVfil <- approxfun(Variables_df$TIME, Variables_df$Vfil, rule = 2)
varVSk <- approxfun(Variables_df$TIME, Variables_df$VSk, rule = 2)
varVsurf <- approxfun(Variables_df$TIME, Variables_df$Vsurf, rule = 2)
varVsc <- approxfun(Variables_df$TIME, Variables_df$Vsc, rule = 2)
varVve <- approxfun(Variables_df$TIME, Variables_df$Vve, rule = 2)
varCLsc <- approxfun(Variables_df$TIME, Variables_df$CLsc, rule = 2)
varCLve <- approxfun(Variables_df$TIME, Variables_df$CLve, rule = 2)
varCLcell <- approxfun(Variables_df$TIME, Variables_df$CLcell, rule = 2)

varVlun <- approxfun(Variables_df$TIME, Variables_df$Vlun, rule = 2)
varVR <- approxfun(Variables_df$TIME, Variables_df$VR, rule = 2)

varOraldose <- approxfun(Variables_df$TIME, Variables_df$Oraldose, rule = 2)
varDrinkdose <- approxfun(Variables_df$TIME, Variables_df$Drinkdose, rule = 2)
varCsurf <- approxfun(Variables_df$TIME, Variables_df$Csurf, rule = 2)

varTm <- approxfun(Variables_df$TIME, Variables_df$Tm, rule = 2)

#### MODELS ----
# manual adaptation of Westerhout 2024 model, removing cholesterol
#### PFOA and PFOS PBK model ----
PFAS_generic <- function(t, A, parms) {
  with(as.list(c(A, parms)), {
    ## PFOA kinetics --
    kurine <- varkurine(t)
    QCP <- varQCP(t)
    QG <- varQG(t)
    QL <- varQL(t)
    QF <- varQF(t)
    QK <- varQK(t)
    Qfil <- varQfil(t)
    QSk <- varQSk(t)
    QR <- varQR(t)
    Vart_plas <- varVart_plas(t)
    Vven_plas <- varVven_plas(t)
    VG <- varVG(t)
    VL <- varVL(t)
    VF <- varVF(t)
    VK <- varVK(t)
    Vfil <- varVfil(t)
    VSk <- varVSk(t)
    Vsc <- varVsc(t)
    Vve <- varVve(t)
    CLsc <- varCLsc(t)
    CLve <- varCLve(t)
    CLcell <- varCLcell(t)
    Vlun <- varVlun(t)
    VR <- varVR(t)
    Oraldose <- varOraldose(t) # Dose expressed in ug/kg/day
    Drinkdose <- varDrinkdose(t) # Dose expressed in ug/kg/day
    Csurf <- varCsurf(t) # Dose expressed in ug/kg/day
    Tm <- varTm(t)
    
    # Concentrations
    # PFOA
    #CAFree <- APlas/VPlas # free concentration of chemical in plasma µg/L (ng/mL) 
    #CA <- CAFree/Free # total concentration of chemical in plasma 
    CG <- AG/VG # Concentration in gut (µg/L) 
    CVG <- CG/PG # Concentration leaving gut (µg/L) 
    CL = AL/VL # Concentration in liver (µg/L) 
    CVL = CL/PL # Concentration leaving liver (µg/L) 
    CF = AF/VF # Concentration in fat (µg/L)
    CVF = CF/PF # Concentration leaving fat (µg/L) 
    CK = AK/VK # Concentration in kidneys (µg/L) 
    CVK = CK/PK # Concentration leaving kidneys (µg/L) 
    Cfil = Afil/Vfil # Concentration in filtrate compartment (µg/L) 
    CSk = ASk/VSk # Concentration in skin compartment (µg/L) 
    CVSk = CSk/PSk # Concentration leaving skin compartment (µg/L)  
    CR = AR/VR # Concentration in rest of the body (µg/L) 
    CVR = CR/PR # Concentration leaving rest of the body (µg/L) 
    Clun <- Alun/Vlun
    CVlun <- Clun/PLun
    Csc <- Asc/Vsc
    Ctrans <- Atrans/0.0000001
    Cve <- Ave/Vve
    CvenFree <- Aven/Vven_plas
    Cven <- CvenFree/Free
    CartFree <- Aart/Vart_plas
    Cart <- CartFree/Free
    
    # Lung compartment
    dAlun <- QCP*Cven*Free - QCP*CVlun*Free
    
    # Stratum corneum
    dAsc <- 2*CLsc*Csurf - 4*CLsc*Csc + 2*CLsc*Ctrans
    
    # Transfer compartment
    dAtrans <- 2*CLsc*Csc - 2*CLsc*Ctrans - 2*CLve*(Ctrans/Kscve) + 2*CLve*Cve
    
    # Viable epidermis
    dAve <- 2*CLve*(Ctrans/Kscve) - 2*CLve*Cve - CLcell*Cve/Kver
    
    # Skin compartment
    dASk <- CLcell*Cve/Kver + QSk*(Cart*Free-CSk*FreeSk)    # Rate of change in skin(µg/h) 
    
    # Venous blood (plasma) compartment      
    dAven <- QF*CF*FreeF + (QL+QG)*CL*FreeL + QR*CR*FreeR + QSk*CSk*FreeSk + QK*CK*FreeK - QCP*Cven*Free
    
    # Arterial blood (plasma) compartment      
    dAart <- QCP*CVlun*Free - QCP*Cart*Free - Qfil*Cart*Free 
    
    # Gut compartment 
    dAG <- QG*(Cart*Free-CG*FreeG) + Oraldose + Drinkdose 
    
    # Liver compartment 
    dAL <- (QL*(Cart*Free)) + (QG*CG*FreeG) - ((QL+QG)*CL*FreeL) # Rate of change in liver (ug/h) 
    
    # Fat compartment 
    dAF <- QF*(Cart*Free-CF*FreeF)   # Rate of change in fat (µg/h) 
    
    # Kidney compartment 
    dAK <- QK*(Cart*Free-CK*FreeK) + (Tm*Cfil)/(Kt+Cfil) # Rate of change in kidneys (µg/h) 
    
    # Filtrate compartment 
    dAfil <- Qfil*(Cart*Free-Cfil) - (Tm*Cfil)/(Kt+Cfil) # Rate of change in filtrate compartment (ug/h) 
    
    # Storage compartment for urine 
    dAdelay <- Qfil*Cfil - kurine*Adelay   
    
    # Urine 
    dAurine <- kurine*Adelay 
    
    # Rest of the body 
    dAR <- QR*(Cart*Free-CR*FreeR)   # Rate of change in rest of the body (µg/h) 
    
    list(c(dAlun, dAsc, dAtrans, dAve, dASk, dAven, dAart,
           dAG, dAL, dAF, dAK, dAfil, dAdelay, dAurine, dAR),
         c(Cart=Cart, Cven=Cven, CVlun=CVlun,
           Csurf=Csurf,Csc=Csc,Ctrans=Ctrans,Cve=Cve, CSk=CSk,
           CG=CG, CVG=CVG, CL=CL, CVL=CVL, CF=CF, CVF=CVF,
           CK=CK, CVK=CVK, Cfil=Cfil, CR=CR, CVR=CVR))
  })
}

#### parms ----
# Parameters used in the model
parms_PFAS_generic <- c(Kt, Free, FreeL, FreeF, FreeK, FreeSk, FreeR, FreeG,
                        PL, PF, PK, PSk, PR, PG, PLun,
                        Kscve,Kver)

#### A_init ----
# Initial values
A_init_PFAS_generic <- c(Alun = ALunbirth,
                         Asc = 0,
                         Atrans = 0,
                         Ave = 0,
                         ASk = ASkbirth,
                         Aven = Avenbirth,
                         Aart = Aartbirth,
                         AG = AGbirth,
                         AL = ALbirth, 
                         AF = AFbirth,
                         AK = AKbirth,
                         Afil = Afilbirth,
                         Adelay = Adelaybirth,
                         Aurine = Aurinebirth,
                         AR = ARbirth)

#### OUTPUT ----
# PFOA
output_PFAS <- lsoda(A_init_PFAS_generic, TIME, PFAS_generic, parms_PFAS_generic)
output.df <- as.data.frame(output_PFAS)
output.df <- Variables_df %>%
  rename(time = TIME) %>%
  left_join(output.df)

#### PLOTS ----
# BW over time
Figure_BW <- ggplot() +
  geom_line(data=output.df, aes(x=age, y=BW), color='black') +
  
  theme_bw() +
  labs(title="BW over time",x="\nAge (y)", y="Bodyweight (kg)\n") +
  theme(plot.title = element_text(hjust = 0.5))
Figure_BW

# Oraldose
Figure_Oraldose <- ggplot() +
  geom_line(data=output.df, aes(x=age, y=Oraldose), color='black') +
  
  theme_bw() +
  labs(title="Oral dose",x="\nAge (y)", y="Dose (\u03BCg)\n") +
  theme(plot.title = element_text(hjust = 0.5))
Figure_Oraldose

# Dermaldose
Figure_Dermaldose <- ggplot() +
  geom_line(data=output.df, aes(x=age, y=Dermaldose), color='black') +
  
  theme_bw() +
  labs(title="Dermal dose",x="\nAge (y)", y="Dose (\u03BCg)\n") +
  theme(plot.title = element_text(hjust = 0.5))
Figure_Dermaldose

# Csurf
Figure_Csurf <- ggplot() +
  geom_line(data=output.df, aes(x=age, y=Csurf), color='black') +
  
  theme_bw() +
  labs(title="Csurf",x="\nAge (y)", y="Amount (\u03BCg)\n") +
  theme(plot.title = element_text(hjust = 0.5))
Figure_Csurf

# Amount in stratum corneum
Figure_Asc <- ggplot() +
  geom_line(data=output.df, aes(x=age, y=Asc), color='black') +
  
  theme_bw() +
  labs(title="Amount in stratum corneum",x="\nAge (y)", y="Amount (\u03BCg)\n") +
  theme(plot.title = element_text(hjust = 0.5))
Figure_Asc

# Amount in gut
Figure_AG <- ggplot() +
  geom_line(data=output.df, aes(x=age, y=AG), color='black') +
  
  theme_bw() +
  labs(title="Amount in gut",x="\nAge (y)", y="Amount (\u03BCg)\n") +
  theme(plot.title = element_text(hjust = 0.5))
Figure_AG

# Amount in liver
Figure_AL <- ggplot() +
  geom_line(data=output.df, aes(x=age, y=AL), color='black') +
  
  theme_bw() +
  labs(title="Amount in liver",x="\nAge (y)", y="Amount (\u03BCg)\n") +
  theme(plot.title = element_text(hjust = 0.5))
Figure_AL

# Amount in skin
Figure_ASk <- ggplot() +
  geom_line(data=output.df, aes(x=age, y=ASk), color='black') +
  
  theme_bw() +
  labs(title="Amount in skin",x="\nAge (y)", y="Amount (\u03BCg)\n") +
  theme(plot.title = element_text(hjust = 0.5))
Figure_ASk

# Concentration in plasma (ug/L)
Figure_PFAS <- ggplot() +
  geom_line(data=output.df, aes(x=age, y=Cart, color="Cart", lty="Cart")) +
  geom_line(data=output.df, aes(x=age, y=Cven, color="Cven", lty="Cven")) +
  
  scale_colour_manual(name='Output',
                      values=c('Cart'='red',
                               'Cven'='blue'),
                      labels=c('Cart'='Concentration in arterial plasma',
                               'Cven'='Concentration in venous plasma')) +
  
  scale_linetype_manual(name='Output',
                        values=c('Cart'='solid',
                                 'Cven'='solid'),
                        labels=c('Cart'='Concentration in arterial plasma',
                                 'Cven'='Concentration in venous plasma')) +
  
  theme_bw() +
  labs(title="PFAS concentration-time profiles",x="\nAge (y)", y="Concentration (\u03BCg/L)\n") +
  theme(plot.title = element_text(hjust = 0.5))

Figure_PFAS

# tiff("Westerhout 2024 - Figure PFAS kinetics.tif",
#      res=600, compression = "lzw",
#      height=6,
#      width=8,
#      units="in")
# Figure_PFAS
# 
# dev.off()
