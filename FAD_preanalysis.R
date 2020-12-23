#### FAD_dataAnalysis ####
#
###Purpose-------------------------------------------------
#
#code for preprocessing the questionnaire data from the 
#revision of a Chinese version of Free Will and Determinism 
#Plus scale 
# 
#
# for prerigistration
# 
# Qing-Lan Liu
#
###-------------------------------------------------

rm(list = ls())     # remove all variables

# get the directory of the current file and set working directory to the folder
curDir = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(curDir)

pkgTest <- function(x)
{
  if(!require(x,character.only = TRUE))
  {
    install.packages(x,dep = TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}
#packages
pkgNeeded <- (c("tidyverse","psych","xlsx","magrittr","stats","lavaan",
                "semTools","semPlot"))

lapply(pkgNeeded,pkgTest)
rm('pkgNeeded') # remove the variable 'pkgNeeded';



####----------------------------Section1------------------------------
#### this section aims to compare EFA factor loadings of 3 Chinese 
#### translations, English version, Japanese version and French version.
####------------------------------------------------------------------


####cleaning data-----------------------------------------------------
#get the data
df1 <- read.table("FAD_1.csv",header = TRUE, sep = ",")    #n=8248, version 1 of Chinese FAD+                  
df2 <- read.table("FAD_2.csv",header = TRUE, sep = ",")    #n=1333,version 2 of Chinese FAD+,excluded 141,check item, age < 18,not native,not complete,
df3 <- read.table("FAD_3.csv",header = TRUE, sep = ",")    #n=715, version 3 of Chinese FAD+, high school students
df4 <- read.table("FAD_A.csv",header = TRUE, sep = ",")    #n=419 (collected from OSF,America data),exluded 8 cases, time < 60 sec
df5 <- read.table ("FAD_F.csv",header = TRUE, sep = ",")   #n=711, French FAD+
df6 <- read.table ("FAD_J.csv",header = TRUE, sep = ",")   #n=3000,Japanese FAD+
#recode df1
df1 <- df1 %>%                                                                # recode gender_gene,gender_self_report
  mutate_at(c("gender_gene","gender_self_report"),
            funs(dplyr::recode(., `male`=1,`female`=2)))
df1 <- df1 %>%                                                                # recode education level,father's eudcation level,mother's education level
  mutate_at(c("education_level","father_education_level","mother_education_level"),
            funs(dplyr::recode(., `小学及以下`=1,`初中`=2,`中专或职高`=3,`高中`=4,
                               `大专`=5,`本科`=6,`硕士`=7,`博士`=8)))
#note:"小学及以下"=primary schoool or less，"初中"=middle school，"中专或职高"=secondary school or specialized school，
#"高中"=high school，"大专"=college，"本科"= bachelor，"硕士"=master，"博士"=doctor
df1 <- df1 %>%                                                                # recode education level,father's eudcation level,mother's education level
  mutate_at(c("birth_level"),
            funs(dplyr::recode(., `大城市`=1,`中等城市`=2,`小城市`=3,`县城`=4,
                               `镇`=5,`村`=6)))


####descriptive----------------------------------------------------------
psych::describe(df1$age)  #mean = 26.51, sd = 5.9
psych::describe(df2$age)  #mean = 23.22, sd = 5.41
gendertable1 <- with(df1, table(gender_gene,education_level))
print(gendertable1) #male = 2449, female = 2192; edul1=14,edul2=38,edul3=51,edul4=95,edul5=421,edul6=2245,edul7=705,edul8=109
gendertable2 <- with(df2, table(gender,edu))             
print(gendertable2) #male = 548,female = 691;eduL1=105,eduL2=80,eduL3=27,eduL4=45,eduL5=354,eduL6=101,eduL7=27



####EFA-----------------------------------------------------------------
#EFA for df1
df1.fad <- subset(df1, select = c(FD1:UP27))
df1.fadcor <- cor(df1.fad,use = "pairwise.complete.obs")   #pairewise deletion for missing value
psych::fa(df1.fadcor, nfactors = 4, n.obs = 8248, rotate = "oblimin", fm = "ml")   #method: oblimin rotation and maximun likelihood extraction

#EFA for df2
df2.fad <- subset(df2, select = c(FD1:FW21,SD22:UP27))
df2.fadcor <- cor(df2.fad,use = "complete.obs")   #listwise for missing value
psych::fa(df2.fadcor, nfactors = 4, n.obs = 1333, rotate = "oblimin", fm = "ml")   #method: oblimin rotation and maximun likelihood extraction

#EFA for df3
df3.fad <- subset(df3, select = c(FD1:UP27))
df3.fadcor <- cor(df3.fad,use = "complete.obs")   #listwise for missing value
psych::fa(df3.fadcor, nfactors = 4, n.obs = 714, rotate = "oblimin", fm = "ml")   #method: oblimin rotation and maximun likelihood extraction

#EFA for df4
df4.fad <- subset(df4,select = c(FD1:UP27))
df4.factor <- cor(df4.fad, use = "complete.obs")  #listwise for missing value
psych::fa(df4.factor, nfactors = 4, n.obs = 711, rotate = "oblimin", fm = "ml")  #method: oblimin rotation and maximun likelihood extraction

#EFA for df5
df5.fad <- subset(df5,select = c(FD1:UP27))
df5.factor <- cor(df5.fad, use = "complete.obs")  #listwise for missing value
psych::fa(df5.factor, nfactors = 4, n.obs = 3000, rotate = "oblimin", fm = "ml")  #method: oblimin rotation and maximun likelihood extraction



####----------------------------Section2------------------------------
#### this section aims to evaluating psychological properties of the final 
#### Chinese adoption, including reliability (alpha, omega), validity (CFA,
#### correlation with BFI, locus of control） and measurement invariance
####------------------------------------------------------------------


####cleaning data-----------------------------------------------------
#load data
df2.1 <- read.table("FAD_2.csv",header = TRUE, sep = ",")    #Chinese version 1 FAD+                  
df2.2 <- read.table("FAD_A.csv",header = TRUE, sep = ",")    #English FAD+
df2.3 <- read.table("FAD_J.csv",header = TRUE, sep = ",")    #Japanese FAD+
df2.4 <- read.table ("FAD_F.csv",header = TRUE, sep = ",")   #French FAD+
df2.5 <- read.table("FAD_1.csv",header = TRUE, sep = ",")    # 
df2.6 <- read.table("FAD_3.csv",header = TRUE, sep = ",")    #                 

#rename columns
demoNames <- c('age','gender','edu','famIncome','faEdu','moEdu','faOccu','moOccu')
BFI_names <- c(paste('BFI_A',1:9,sep = ''),
               paste('BFI_C',1:9,sep = ''),
               paste('BFI_N',1:8,sep = ''),
               paste('BFI_O',1:10,sep = ''),
               paste('BFI_E',1:8,sep = ''))
IPC_Names <- c(paste("IPC",1:24,sep = '_')) 
FAD_Name  <- c("FD1", "SD2", "UP3", "FW4", "FD5", "SD6", "UP7",
               "FW8", "FD9", "SD10", "UP11", "FW12", "FD13", 
               "SD14", "UP15", "FW16", "FD17", "SD18", "UP19", 
               "UP20","FW21", "SD22", "FW23", "SD24", "UP25", 
               "FW26", "UP27")
reFADNames <- c("reFD1","reSD2","reUP3","reFW4","reFD5","reSD6",               
                "reUP7","reFW8","reFD9","reSD10","reUP11","reFW12",  
                "reFD13","reSD14", "reUP15", "reFW16","reFD17",   
                "reSD18", "reUP19","reUP20","reFW21", "reSD22", 
                "reFW23","reSD24", "reUP25","reFW26", "reUP27")


colnames(df2.1)
df2.1 <- df2.1 %>%
  dplyr::filter(complete.cases(subjID)) %>%
  dplyr::select(-which(names(df2.1) %in% c('ip','email','weixin_nickname','weixin_sex',
                                         'weixin_addr','start','finish','status')))
names(df2.1)[2:35] <- FAD_Name
names(df2.1)[36:79] <- BFI_names    #精简题目
names(df2.1)[80]   <- checkItem
names(df2.1)[81:104] <- IPC_Names
names(df2.1)[105:112]<- demoNames

colnames(df2.5)
df2.5 <- df2.5 %>%
  dplyr::filter(complete.cases(id)) %>%
  dplyr::select(-which(names(df2.5) %in% c('ip','email','weixin_nickname','weixin_sex',
                                           'weixin_addr','start','finish','status')))
names(df2.5)[2:35] <- reFADNames
names(df2.1)[36:43]<- demoNames

#check item
df2.1 <- subset(df2.1, check == 3 & age > 17)             # check the minimal attention and age 
df2.5 <- subset(df2.5, check == 3 & age > 17)   

#merge the test and retest
colnames(df2.1) 
first.col <- c("subjID","gender","age","edu","faEdu",
               "moEdu","faOccu","moOccu",'famIncome')
df2.1FAD <- df2.1[,c(first.col, dplyr::setdiff(colnames(df2.1), first.col))]

colnames(df2.5)
df2.5FAD <- df2.5[,c(first.col, dplyr::setdiff(colnames(df2.5), first.col))]

length(intersect(df2.1FAD$subjID,df2.5FAD$subjID))    
intersect(colnames(df2.1FAD),colnames(df2.5FAD))        

dfFAD.total <- dplyr::full_join(x = df2.1FAD, y = df2.5FAD, by = intersect(colnames(df2.1FAD),colnames(df2.5FAD)))    


####descriptive------------------------------------------------------------
#scores
FDNames <- c("FD1","FD5","FD9","FD13", "FD17")
SDNames <- c("SD2","SD6","SD10","SD14","SD18","SD22","SD24")
UPNames <- c("UP3","UP7","UP11","UP15","UP19","UP20","UP25","UP27")
FWNames <- c("FW4","FW8","FW12","FW16","FW21","FW23","FW26")

df1$FD <- rowSums(df1[,FDNames],na.rm = F)/length(FDNames)
df1$SD <- rowSums(df1[,SDNames],na.rm = F)/length(SDNames)
df1$UP <- rowSums(df1[,UPNames],na.rm = F)/length(UPNames)
df1$FW <- rowSums(df1[,FWNames],na.rm = F)/length(FWNames)

FADNames <- c("FD","SD","UP","FW")
describe(df2.1[,FADNames])  
#difference between gender
t.test(df1$FD ~ df1$gender_gene)        
t.test(df1$SD ~ df1$gender_gene)       
t.test(df1$UP ~ df1$gender_gene)        
t.test(df1$FW ~ df1$gender_gene)        


####reliability-------------------------------------------------------------------
#alpha, omega, test-retest
FDNames 
reFDNames <- c("reFD1","reFD5","reFD9","reFD13", "reFD17")
FDKeys <- c(1,2,3,4,5)

df2.1$FD <- rowSums(df2.1[,FDNames],na.rm = F)/length(FDNames) # average score for time 1
df2.5$reFD <- rowSums(df2.5[,reFDNames],na.rm = F)/length(reFDNames) # average score for time 2
FDAlpha <-  psych::alpha(df2.1[,FDNames], keys=FDKeys)  # calculate the alpha coefficient of Fatalistic Determinism
print(FDAlpha$total)  # print the alpha for Fatalistic Determinism

FDOmega <- psych::omega(df2.1[,FDNames])  
print(c(FDOmega$omega_h,FDOmega$omega.tot)) 

stats::cor.test(dfFAD.total$FD,dfFAD.total$reFD,alternative = "greater") # caculate test-retest reliability

SDNames 
reSDNames <- c("reSD2","reSD6","reSD10","reSD14","reSD18","reSD22","reSD24")
SDKeys <- c(1,2,3,4,5,6,7)

df2.1$SD <- rowSums(df2.1[,SDNames],na.rm = F)/length(SDNames) # average score for time 1
df2.5$reSD <- rowSums(df2.5[,reSDNames],na.rm = F)/length(reSDNames) # average score for time 2
SDAlpha <-  psych::alpha(df2.1[,SDNames], keys=SDKeys)  # calculate the alpha coefficient of Scientific Determinism
print(SDAlpha$total)  # print the alpha for Scientific Determinism

SDOmega <- psych::omega(df2.1[,SDNames])  
print(c(SDOmega$omega_h,SDOmega$omega.tot)) 

stats::cor.test(dfFAD.total$SD,dfFAD.total$reSD,alternative = "greater") # caculate test-retest reliability


UPNames 
reUPNames <- c("reUP3","reUP7","reUP11","reUP15","reUP19","reUP20","reUP25","reUP27")
UPKeys <- c(1,2,3,4,5,6,7,8)

df2.1$UP <- rowSums(df2.1[,UPNames],na.rm = F)/length(UPNames) # average score for time 1
df2.5$reUP <- rowSums(df2.5[,reUPNames],na.rm = F)/length(reUPNames) # average score for time 2
UPAlpha <-  psych::alpha(df2.1[,UPNames], keys=UPKeys)  # calculate the alpha coefficient of Unpredictability
print(UPAlpha$total)  # print the alpha for Unpredictability

UPOmega <- psych::omega(df2.1[,UPNames])  
print(c(UPOmega$omega_h,UPOmega$omega.tot)) 

FWNames
reFWNames <- c("reFW4","reFW8","reFW12","reFW16","reFW21","reFW23","reFW26")
FWKeys <- c(1,2,3,4,5,6,7)

df2.1$FW <- rowSums(df2.1[,FWNames],na.rm = F)/length(FWNames) # average score for time 1
df2.1$reFW <- rowSums(df2.1[,reFWNames],na.rm = F)/length(reFWNames) # average score for time 2
FWAlpha <-  psych::alpha(df2.1[,FWNames], keys=FWKeys)  # calculate the alpha coefficient of Free Will
print(FWAlpha$total)  # print the alpha for Free Will

FWOmega <- psych::omega(df2.1[,FWNames])  
print(c(FWOmega$omega_h,FWOmega$omega.tot)) 

stats::cor.test(dfFAD.total$FW,dfFAD.total$reFW,alternative = "greater") # caculate test-retest reliability

#reliability of the Big Five Inventory
BFI_ANames <- c("BFI_A1","BFI_A2","BFI_A3","BFI_A4","BFI_A5","BFI_A6","BFI_A7","BFI_A8","BFI_A9")
BFI_AKeys <- c(1,1,1,1,1,-1,-1,-1,-1)         # reverse coded as negative

BFI_AAlpha <-  psych::alpha(df2.1[,BFI_ANames], keys=BFI_AKeys)  # calculate the alpha coefficient of Agreeableness
print(BFI_AAlpha$total)  # print the alpha for Agreeableness

BFI_AOmega <- psych::omega(df2.1[,BFI_ANames])  
print(c(BFI_AOmega$omega_h,BFI_AOmega$omega.tot)) 


BFI_CNames <- c("BFI_C1","BFI_C2","BFI_C3","BFI_C4","BFI_C5","BFI_C6","BFI_C7","BFI_C8","BFI_C9")
BFI_CKeys <- c(1,1,1,1,1,-1,-1,-1,-1)         # reverse coded as negative

BFI_CAlpha <-  psych::alpha(df2.1[,BFI_CNames], keys=BFI_CKeys)  # calculate the alpha coefficient of Conscientiousness
print(BFI_CAlpha$total)  # print the alpha for Conscientiousness

BFI_COmega <- psych::omega(df2.1[,BFI_CNames])  
print(c(BFI_COmega$omega_h,BFI_COmega$omega.tot)) 


BFI_NNames <- c("BFI_N1","BFI_N2","BFI_N3","BFI_N4","BFI_N5","BFI_N6","BFI_N7","BFI_N8")
BFI_NKeys <- c(1,1,1,1,1,-1,-1,-1)         # reverse coded as negative

BFI_NAlpha <-  psych::alpha(df2.1[,BFI_NNames], keys=BFI_NKeys)  # calculate the alpha coefficient of Neuroticism
print(BFI_NAlpha$total)  # print the alpha for Neuroticism

BFI_NOmega <- psych::omega(df2.1[,BFI_NNames])  
print(c(BFI_NOmega$omega_h,BFI_NOmega$omega.tot)) 


BFI_ONames <- c("BFI_O1","BFI_O2","BFI_O3","BFI_O4","BFI_O5","BFI_O6","BFI_O7","BFI_O8","BFI_O9","BFI_O10")
BFI_OKeys <- c(1,1,1,1,1,1,1,1,1,1)         # reverse coded as negative

BFI_OAlpha <-  psych::alpha(df2.1[,BFI_ONames], keys=BFI_OKeys)  # calculate the alpha coefficient of Openness
print(BFI_OAlpha$total)  # print the alpha for Openness


BFI_OOmega <- psych::omega(df2.1[,BFI_ONames])  
print(c(BFI_OOmega$omega_h,BFI_OOmega$omega.tot)) 


BFI_ENames <- c("BFI_E1","BFI_E2","BFI_E3","BFI_E4","BFI_E5","BFI_E6","BFI_E7","BFI_E8")
BFI_EKeys <- c(1,1,1,1,1,1,1,1)      # reverse coded as negative   

BFI_EAlpha <-  psych::alpha(df2.1[,BFI_ENames], keys=BFI_EKeys)  # calculate the alpha coefficient of Extraversion
print(BFI_EAlpha$total)  # print the alpha for Extraversion

BFI_EOmega <- psych::omega(df2.1[,BFI_ENames])  
print(c(BFI_EOmega$omega_h,BFI_EOmega$omega.tot)) 


#from the MLOC
MLOC_INames <- c("MLOC1","MLOC4","MLOC5","MLOC9","MLOC18","MLOC19","MLOC21","MLOC23")
MLOC_IKeys <-  c(1,1,1,1,1,1,1,1)        

MLOC_IAlpha <-  psych::alpha(df2.1[,MLOC_INames], keys=MLOC_IKeys)  # calculate the alpha coefficient of Internal
print(MLOC_IAlpha$total)  # print the alpha for Internal

MLOC_IOmega <- psych::omega(df2.1[,MLOC_INames])  
print(c(MLOC_IOmega$omega_h,MLOC_IOmega$omega.tot)) 


MLOC_PNames <- c("MLOC3","MLOC8","MLOC11","MLOC13","MLOC15","MLOC17","MLOC20","MLOC22")
MLOC_PKeys <-  c(1,1,1,1,1,1,1,1)        

MLOC_PAlpha <-  psych::alpha(df2.1[,MLOC_PNames], keys=MLOC_PKeys)     # calculate the alpha coefficient of Powerful others
print(MLOC_PAlpha$total)  # print the alpha for Powerful others

MLOC_POmega <- psych::omega(df2.1[,MLOC_PNames])  
print(c(MLOC_POmega$omega_h,MLOC_POmega$omega.tot)) 


MLOC_CNames <- c("MLOC2","MLOC6","MLOC7","MLOC10","MLOC12","MLOC14","MLOC16","MLOC24")
MLOC_CKeys <-  c(1,1,1,1,1,1,1,1)        

MLOC_CAlpha <-  psych::alpha(df2.1[,MLOC_CNames], keys=MLOC_CKeys)     # calculate the alpha coefficient of Chance
print(MLOC_CAlpha$total)  # print the alpha for Chance

MLOC_COmega <- psych::omega(df2.1[,MLOC_CNames])  
print(c(MLOC_COmega$omega_h,MLOC_COmega$omega.tot)) 

####validity------------------------------------------------------------------------
#CFA
model <- 'FD =~ FD1 + FD5 + FD9 + FD13 + FD17;
          SD =~ SD2 + SD6 + SD10 + SD14 + SD18 + SD22 + SD24;
          UP =~ UP3 + UP7 + UP11 + UP15 + UP19 + UP20 + UP25 + UP27;
          FW =~ FW4 + FW8 + FW12 + FW16 + FW21 + FW23 + FW26'

fit <- cfa(model, data = df2.1)
summary(fit, standardized = TRUE, fit.measures = TRUE)
moreFitIndices(fit)
#correlation
df2.1$BFI_A <- rowSums(df2.1[,BFI_ANames])/length(BFI_ANames)    # scoring BFI
df2.1$BFI_C <- rowSums(df2.1[,BFI_CNames])/length(BFI_CNames) 
df2.1$BFI_N <- rowSums(df2.1[,BFI_NNames])/length(BFI_NNames) 
df2.1$BFI_O <- rowSums(df2.1[,BFI_ONames])/length(BFI_ONames) 
df2.1$BFI_E <- rowSums(df2.1[,BFI_ENames])/length(BFI_ENames) 

df2.1$MLOC_I <- (rowSums(df2.1[,MLOC_INames])+24)/length(MLOC_INames)  #scoring IPC
df2.1$MLOC_P <- (rowSums(df2.1[,MLOC_PNames])+24)/length(MLOC_PNames)  
df2.1$MLOC_C <- (rowSums(df2.1[,MLOC_CNames])+24)/length(MLOC_CNames)  

corScale <- subset(df2.1[,c('FD','SD','UP','FW','BFI_A',
                            'BFI_C','BFI_N','BFI_O',
                            'BFI_E','MLOC_I','MLOC_P'
                            ,'MLOC_C')])

model2 <- 'FAD =~ FD + SD + UP + FW;
          BFI =~ BFI_C + BFI_N + BFI_O + BFI_E + BFI_A;
          IPC =~ IPC_I + IPC_P + IPC_C'

fit <- cfa(model, data = df2.1)
summary(fit, standardized = TRUE, fit.measures = TRUE)
moreFitIndices(fit)

#measurement invariance
df2.11 <- df2.1[FAD_Name] %>%
  dplyr::mutate(ver = 1)
df2.51 <- df2.5 [FAD_Name]%>%
  dplyr::mutate(ver = 2)
df2.61 <- df2.6[FAD_Name] %>%
  dplyr::mutate(ver = 3)


intersect(df2.11$ver,df2.51$ver)     
length(intersect(colnames(df2.11),colnames(df2.51)))  
fad1 <- dplyr::full_join(x = df2.11, y = df2.51, by = intersect(colnames(df2.11),colnames(df2.51))) 


intersect(fad1$ver,df2.61$ver)     
length(intersect(colnames(fad1),colnames(df2.61)))  
fad2 <- dplyr::full_join(x = fad1, y = df2.61, by = intersect(colnames(fad1),colnames(df2.61))) 


str(fad2)

Mgcfa <- cfa(model, data = fad2, group = "ver")
weak <- cfa(model, data = fad2, group = "ver", group.equal = "loadings")
strong <- cfa(model, data = fad2, group = "ver", group.equal = c("loadings","intercepts"))
strict <- cfa(model, data = fad2, group = "ver", group.equal = c("loadings","intercepts","residuals"))

anova(Mgcfa, weak, strong, strict)
compareFit(Mgcfa, weak, nested = FALSE)
compareFit(weak, strong, nested = FALSE)
compareFit(strong, strict, nested = FALSE)

modWeak <- modificationIndices(weak)
modWeak[modWeak$op == "~1",]

models <- list(fit.configural = Mgcfa, fit.loadings = weak)
partialInvariance(models, "metric", free = "UP3")
