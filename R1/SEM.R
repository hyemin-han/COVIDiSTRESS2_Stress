# SEM

library(lavaan)
library(blavaan)

# 
load('../Stress_aligned.RData')


model.res <-'

GID =~ identity_1_0neutral + identity_2_0neutral+
  identity_3_0neutral+identity_4_0neutral

SS =~ secondary_stressors__1 + secondary_stressors__2+secondary_stressors__3+
  secondary_stressors__4

RES =~ resilience_1 + resilience_2+resilience_3+
  resilience_4+resilience_5+resilience_6 + GID + SS
  
SS ~~ GID
'
fit.res <- sem(model.res,data=data.filtered, estimator='DWLS')
fitmeasures(fit.res)

model.pss <- '

GID =~ identity_1_0neutral + identity_2_0neutral+
  identity_3_0neutral+identity_4_0neutral

SS =~ secondary_stressors__1 + secondary_stressors__2+secondary_stressors__3+
  secondary_stressors__4

PSS =~ perceived_stress_sca_1 + perceived_stress_sca_2+
  perceived_stress_sca_3 + perceived_stress_sca_4 + perceived_stress_sca_5+
  perceived_stress_sca_6 + perceived_stress_sca_7 + perceived_stress_sca_8+
  perceived_stress_sca_9 + perceived_stress_sca_10 + GID + SS
  
SS ~~ GID
'

fit.pss <- sem(model.pss,data=data.filtered,estimator='DWLS')
fitmeasures(fit.pss)

