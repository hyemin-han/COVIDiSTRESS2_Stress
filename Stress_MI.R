library(psych)
library(lavaan)
library(sirt)
library(MASS)

# function for factor score adjustment
aligned.factor.scores <- function(lambda,nu,y){
  # calculate inverse matrix
  lambda1 <- ginv((lambda))
  # create matrix for nu
  ns <- nrow(y)
  nus <- matrix(nu,nrow=ns,ncol=length(nu),byrow=T)
  # y - nu
  y_nu <- y - nus
  # f = inv(lambda)*(y-nu)
  F <- lambda1 %*% t(as.matrix(y_nu))
}

# Load the cleaned csv file

# load data
data<-read.csv('Final_COVIDiSTRESS_Vol2_cleaned.csv')

# item names to be aligned
items.pss <- colnames(data)[45:54]
items.sps <- colnames(data)[76:78]
items.identity <- colnames(data)[173:176]
items.resilience <- colnames(data)[120:125]

# reverse coded items
data[,items.resilience[2]] <- 8-data[,items.resilience[2]]
data[,items.resilience[4]] <- 8-data[,items.resilience[4]]
data[,items.resilience[6]] <- 8-data[,items.resilience[6]]

# extract languages with n >= 100
n.langs <- table(data$UserLanguage)
list.langs <- labels(n.langs)[[1]]
langs.include <- list.langs[n.langs>=100]
n.include <- n.langs[n.langs>=100]

# extract data
for (i in 1:length(langs.include)){
  if (i == 1){
    data.mi <- data[data$UserLanguage == langs.include[i],]
  }else{
    current <- data[data$UserLanguage == langs.include[i],]
    data.mi <- rbind(data.mi,current)
  }
}

# set and examine fitmeasures
fits <- c('rmsea.scaled','srmr','cfi.scaled','tli.scaled')

#####
# 1. PSS

# general CFA: PSS
cfa.model.pss <- 'PSS =~ perceived_stress_sca_1 + perceived_stress_sca_2+
  perceived_stress_sca_3 + perceived_stress_sca_4 + perceived_stress_sca_5+
  perceived_stress_sca_6 + perceived_stress_sca_7 + perceived_stress_sca_8+
  perceived_stress_sca_9 + perceived_stress_sca_10'
cfa.whole.pss <- cfa(model=cfa.model.pss,data=data.mi,estimator='WLSMV', group = 
                       'UserLanguage')
fitMeasures(cfa.whole.pss)[fits]
# msea.scaled         srmr   cfi.scaled   tli.scaled 
#  0.09464355   0.06379443   0.87413322   0.83817129 
# not good -> alignment

# measurement alignment test
# extract parameters
par.pss <- invariance_alignment_cfa_config(dat = data.mi[,items.pss], 
                                       group = data.mi$UserLanguage)
# do alignment
mod1.pss <- invariance.alignment(lambda = par.pss$lambda, nu =
                                   par.pss$nu, align.scale = c(0.2, 0.4), align.pow = c(0.25, 0.25))
# test performance
mod1.pss$es.invariance['R2',]
#loadings intercepts 
#0.9850629  0.9930655  Very good


#####
# 2. sps

# general CFA: sps
cfa.model.sps <- 'SPS =~ perceived_support_1_midneutral + 
  perceived_support_2_midneutral+perceived_support_3_midneutral'
cfa.whole.sps <- cfa(model=cfa.model.sps,data=data.mi,estimator='WLSMV', group = 
                       'UserLanguage')
fitMeasures(cfa.whole.sps)[fits]
# msea.scaled         srmr   cfi.scaled   tli.scaled 
#  0.00000e+00  8.35017e-08  1.00000e+00  1.00000e+00 
# good -> metric

cfa.metric.sps <- cfa(model=cfa.model.sps,data=data.mi,estimator='WLSMV', group = 
                        'UserLanguage', group.equal='loading')
fitMeasures(cfa.metric.sps)[fits]
# msea.scaled         srmr   cfi.scaled   tli.scaled 
#    0.05463580   0.02031403   0.98934654   0.98342795 
# not acceptable due to huge fit indicator change -> alignment

# measurement alignment test
# extract parameters
par.sps <- invariance_alignment_cfa_config(dat = data.mi[,items.sps], 
                                           group = data.mi$UserLanguage)
# do alignment
mod1.sps <- invariance.alignment(lambda = par.sps$lambda, nu =
                                   par.sps$nu, align.scale = c(0.2, 0.4), align.pow = c(0.25, 0.25))
# test performance
mod1.sps$es.invariance['R2',]
#loadings intercepts 
# 0.9924918  0.9987810   Very good


#####
# 3. identity


# general CFA: identity
cfa.model.identity <- 'Identity =~ identity_1_0neutral + identity_2_0neutral+
  identity_3_0neutral+identity_4_0neutral'
cfa.whole.identity<- cfa(model=cfa.model.identity,data=data.mi,estimator='WLSMV', group = 
                       'UserLanguage')
fitMeasures(cfa.whole.identity)[fits]
# msea.scaled         srmr   cfi.scaled   tli.scaled 
#  0.10233495   0.03053738   0.96461385   0.89384156

# measurement alignment test
# extract parameters
par.identity <- invariance_alignment_cfa_config(dat = data.mi[,items.identity], 
                                           group = data.mi$UserLanguage)
# do alignment
mod1.identity <- invariance.alignment(lambda = par.identity$lambda, nu =
                                        par.identity$nu, align.scale = c(0.2, 0.4), align.pow = c(0.25, 0.25))
# test performance
mod1.identity$es.invariance['R2',]
#loadings intercepts 
# 0.9516856  0.9923646     good


#####
# 4. resilience


# general CFA: resilience
cfa.model.resilience <- 'Identity =~ resilience_1 + resilience_2+resilience_3+
  resilience_4+resilience_5+resilience_6'
cfa.whole.resilience<- cfa(model=cfa.model.resilience,data=data.mi,estimator='WLSMV', group = 
                           'UserLanguage')
fitMeasures(cfa.whole.resilience)[fits]
# msea.scaled         srmr   cfi.scaled   tli.scaled 
#  0.09295767   0.03737383   0.95253039   0.92088398 not good

# measurement alignment test
# extract parameters
par.resilience <- invariance_alignment_cfa_config(dat = data.mi[,items.resilience], 
                                                group = data.mi$UserLanguage)
# do alignment
mod1.resilience <- invariance.alignment(lambda = par.resilience$lambda, nu =
                                        par.resilience$nu, align.scale = c(0.2, 0.4), align.pow = c(0.25, 0.25))
# test performance
mod1.resilience$es.invariance['R2',]
#loadings intercepts 
# 0.9799195  0.9971573     good



#####
# factor score calculation

for (i in 1:length(langs.include)){
  if (i == 1){
    # create new matrix
    data.aligned <- data.mi[data.mi$UserLanguage==langs.include[i],]
    # aligned factor score
    F.pss <- aligned.factor.scores(mod1.pss$lambda.aligned[i,],
                               mod1.pss$nu.aligned[i,],
                               data.mi[data.mi$UserLanguage==langs.include[i],items.pss])
    F.sps <- aligned.factor.scores(mod1.sps$lambda.aligned[i,],
                                   mod1.sps$nu.aligned[i,],
                                   data.mi[data.mi$UserLanguage==langs.include[i],items.sps])
    F.id <- aligned.factor.scores(mod1.identity$lambda.aligned[i,],
                                  mod1.identity$nu.aligned[i,],
                                   data.mi[data.mi$UserLanguage==langs.include[i],items.identity])
    F.rs <- aligned.factor.scores(mod1.resilience$lambda.aligned[i,],
                                  mod1.resilience$nu.aligned[i,],
                                  data.mi[data.mi$UserLanguage==langs.include[i],items.resilience])
    data.aligned$pss <- t(F.pss)
    data.aligned$sps <- t(F.sps)
    data.aligned$identity <- t(F.id)
    data.aligned$resilience <- t(F.id)
  }else
  {
    # bind
    current <- data.mi[data.mi$UserLanguage==langs.include[i],]
    F.pss <- aligned.factor.scores(mod1.pss$lambda.aligned[i,],
                                   mod1.pss$nu.aligned[i,],
                               current[,items.pss])
    F.sps <- aligned.factor.scores(mod1.sps$lambda.aligned[i,],
                                   mod1.sps$nu.aligned[i,],
                                   current[,items.sps])
    F.id <- aligned.factor.scores(mod1.identity$lambda.aligned[i,],
                                  mod1.identity$nu.aligned[i,],
                                   current[,items.identity])
    F.rs <- aligned.factor.scores(mod1.resilience$lambda.aligned[i,],
                                  mod1.resilience$nu.aligned[i,],
                                  current[,items.resilience])
    current$pss <- t(F.pss)
    current$sps <- t(F.sps)
    current$identity <- t(F.id)
    current$resilience <- t(F.rs)
    data.aligned <- rbind(data.aligned,current)
  }
}

# save aligned datafile
save.image(file='Stress_aligned.RData')


# alphas
psych::alpha(data[,items.pss],check.keys=TRUE)
psych::alpha(data[,items.sps],check.keys=TRUE)
psych::alpha(data[,items.identity],check.keys=TRUE)
psych::alpha(data[,items.resilience],check.keys=TRUE)
