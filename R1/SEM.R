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

