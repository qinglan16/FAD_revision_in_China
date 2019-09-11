#### FAD_EFA ####
#
###Purpose-------------------------------------------------
#
# code for preprocessing the questionnaire data from the revision of a Chinese version of Free Will and Determinism Plus scale and Experiments on the Perceptual Prioritization of the Good self (FADGS)
# 
#
# 
# 
# 
#
###Preparing-------------------------------------------------

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
pkgNeeded <- (c("tidyverse","psych","xlsx","magrittr","stats"))

lapply(pkgNeeded,pkgTest)
rm('pkgNeeded') # remove the variable 'pkgNeeded';
#get the data
df1 <- read.table("FAD_1.csv",header = TRUE, sep = ",")    #n=8248                  
df2 <- read.table("FAD_2.csv",header = TRUE, sep = ",")    #n=1333,excluded 141,check item, age < 18,not native,not complete,
df3 <- read.table("FAD_3.csv",header = TRUE, sep = ",")    #n=419 (collected from OSF),exluded 8 cases, time < 60 sec
df4 <- read.table ("FAD_F.csv",header = TRUE, sep = ",")   #n=711, data of French FAD
df5 <- read.table ("FAD_J.csv",header = TRUE, sep = ",")   #n=3000, data of Japanese FAD
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
#descriptive
psych::describe(df1$age)  #mean = 26.51, sd = 5.9
psych::describe(df2$age)  #mean = 23.22, sd = 5.41
gendertable1 <- with(df1, table(gender_gene,education_level))
print(gendertable1) #male = 2449, female = 2192; edul1=14,edul2=38,edul3=51,edul4=95,edul5=421,edul6=2245,edul7=705,edul8=109
gendertable2 <- with(df2, table(gender,edu))             
print(gendertable2) #male = 548,female = 691;eduL1=105,eduL2=80,eduL3=27,eduL4=45,eduL5=354,eduL6=101,eduL7=27
#scores
FDNames <- c("FD1","FD5","FD9","FD13", "FD17")
SDNames <- c("SD2","SD6","SD10","SD14","SD18","SD22","SD24")
UPNames <- c("UP3","UP7","UP11","UP15","UP19","UP20","UP25","UP27")
FWNames <- c("FW4","FW8","FW12","FW16","FW21","FW23","FW26")
df1$FD <- rowSums(df1[,FDNames],na.rm = F)/length(FDNames)
df2$FD <- rowSums(df2[,FDNames],na.rm = F)/length(FDNames)
df1$SD <- rowSums(df1[,SDNames],na.rm = F)/length(SDNames)
df2$SD <- rowSums(df2[,SDNames],na.rm = F)/length(SDNames)
df1$UP <- rowSums(df1[,UPNames],na.rm = F)/length(UPNames)
df2$UP <- rowSums(df2[,UPNames],na.rm = F)/length(UPNames)
df1$FW <- rowSums(df1[,FWNames],na.rm = F)/length(FWNames)
df2$FW <- rowSums(df2[,FWNames],na.rm = F)/length(FWNames)
#difference between gender
FADNames <- c("FD","SD","UP","FW")
t.test(df1$FD ~ df1$gender_gene)        #m=1.34,f=1.56,[-0.27, -0.18],p < 2.2e-16
t.test(df1$SD ~ df1$gender_gene)        #m=2.16,f=2.17,[-0.04,0.02],p = 0.62
t.test(df1$UP ~ df1$gender_gene)        #m=2.18,f=2.19,[-0.05,0.02], p = 0.4812
t.test(df1$FW ~ df1$gender_gene)        #m=2.21,f=2.24,[-0.06,0.00],p = 0.0511
t.test(df2$FD ~ df2$gender)             #m=2.35,f=2.38,[-0.11,0.04], p = 0.3652
t.test(df2$SD ~ df2$gender)             #m=3.06,f=3.09,[-0.09,0.03], p = 0.2987
t.test(df2$UP ~ df2$gender)             #m=3.06,f=3.07,[-0.07,0.06], p = 0.8724
t.test(df2$FW ~ df2$gender)             #m=3.27,f=3.38,[-0.18,-0.05], p = 0.0006301
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
psych::fa(df3.fadcor, nfactors = 4, n.obs = 419, rotate = "oblimin", fm = "ml")   #method: oblimin rotation and maximun likelihood extraction


#EFA for df4
df4.fad <- subset(df4,select = c(FD1:UP27))
df4.factor <- cor(df4.fad, use = "complete.obs")  #listwise for missing value
psych::fa(df4.factor, nfactors = 4, n.obs = 711, rotate = "oblimin", fm = "ml")  #method: oblimin rotation and maximun likelihood extraction

#EFA for df5
df5.fad <- subset(df5,select = c(FD1:UP27))
df5.factor <- cor(df5.fad, use = "complete.obs")  #listwise for missing value
psych::fa(df5.factor, nfactors = 4, n.obs = 3000, rotate = "oblimin", fm = "ml")  #method: oblimin rotation and maximun likelihood extraction



