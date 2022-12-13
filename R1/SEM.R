# SEM

library(lavaan)
library(blavaan)
library(fastDummies)

# 
load('../Stress_aligned.RData')


model.res <-'

GID =~ identity_1_0neutral + identity_2_0neutral+
  identity_3_0neutral+identity_4_0neutral

SS =~ secondary_stressors__1 + secondary_stressors__2+secondary_stressors__3+
  secondary_stressors__4
  
  SPS =~ perceived_support_1_0neutral + perceived_support_2_0neutral +
  perceived_support_3_0neutral

RES =~ resilience_1 + resilience_2+resilience_3+
  resilience_4+resilience_5+resilience_6 + GID + SS+
  SPS 
  
  
  
SS ~~ GID
SS ~~ SPS
GID ~~ SPS
'
fit.res <- sem(model.res,data=data.filtered, estimator='DWLS')
fitmeasures(fit.res)

model.pss <- '

GID =~ identity_1_0neutral + identity_2_0neutral+
  identity_3_0neutral+identity_4_0neutral

SS =~ secondary_stressors__1 + secondary_stressors__2+secondary_stressors__3+
  secondary_stressors__4

SPS =~ perceived_support_1_0neutral + perceived_support_2_0neutral +
  perceived_support_3_0neutral

PSS =~ perceived_stress_sca_1 + perceived_stress_sca_2+
  perceived_stress_sca_3 + perceived_stress_sca_4 + perceived_stress_sca_5+
  perceived_stress_sca_6 + perceived_stress_sca_7 + perceived_stress_sca_8+
  perceived_stress_sca_9 + perceived_stress_sca_10 + GID + SS+
  SPS 
                    
  
SS ~~ GID
SPS ~~ GID
SPS ~~ SS
'

fit.pss <- sem(model.pss,data=data.filtered,estimator='DWLS'
               )
fitmeasures(fit.pss)



### correlation (Reviewer 3)
library(psych)
# variable list
vars <- c('primary_stressor_avg','secondary','sps','identity',
          'pss','resilience')
cors<-corr.test(data.filtered[,vars])

# draw sempath
library(semPlot)
library(semTools)
library(semptools)

fig.pss <- semPaths(fit.pss, whatLabels = 'est',
                    sizeMan=5,edge.label.cex = 1.5,label.prop=1.5)
fig.res <- semPaths(fit.res, whatLabels = 'est',
                    sizeMan=5,edge.label.cex = 1.5,
                    style='ram',label.prop=1.5)

# demographics
countries <- length(n.include.c)

# age M/SD
ages <- matrix(nrow=countries, ncol = 2)
# sex
# male female other %
gender <- matrix(nrow = countries, ncol=3)
# education
# PhD/BA/College/12/9/6/None
education <- matrix(nrow=countries,ncol=7)
for (i in 1:countries){
  current <- data.filtered[data.filtered$residing_country == country.include[i],]
  ages[i,1] <- mean(current$age, na.rm = T)
  ages[i,2] <- sd(current$age, na.rm=T)
  current.gender <- table(current$gender)
  gender.tot <- sum(current.gender)
  # male female other
  gender[i,1] <- sum(current$gender=='Male',na.rm = T)/gender.tot*1
  gender[i,2] <- sum(current$gender=='Female',na.rm = T)/gender.tot*1
  gender[i,3] <- sum(current$gender=='Other/Would rather not say',na.rm = T)/
    gender.tot*1
  # education
  current.educ <- table(current$education)
  educ.tot <- sum(current.educ)
  # phd to non
  education[i,1 ] <- sum(current$education=='PhD / Doctorate',na.rm = T)/
    educ.tot*1
  education[i,2 ] <- sum(current$education=='University degree (e.g., MA, MSc, BA, BSc)',na.rm = T)/
    educ.tot*1
  education[i,3 ] <- sum(current$education=='Some University or equivalent \n(still ongoing, or completed a module or more, but did not graduate)',na.rm = T)/
    educ.tot*1
  education[i,4 ] <- sum(current$education=='Up to 12 years of school',na.rm = T)/educ.tot
  education[i,5 ] <- sum(current$education=='Up to 9 years of school',na.rm = T)/educ.tot*1
  education[i,6 ] <- sum(current$education=='Up to 6 years of school',na.rm = T)/educ.tot
  education[i,7 ] <- sum(current$education=='None',na.rm = T)/educ.tot
}
gender[is.na(gender)] <- 0
education[is.na(education)] <-0
# summary
TABLE <- cbind(ages,gender,education)
TABLE <- as.data.frame(TABLE)
rownames(TABLE)<- country.include

# whole dataset
mean(data.filtered$age,na.rm=T)
sd(data.filtered$age,na.rm=T)
table(data.filtered$gender)/sum(table(data.filtered$gender))
table(data.filtered$education)/sum(table(data.filtered$education))