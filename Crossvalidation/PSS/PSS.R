library(psych)
library(lavaan)
library(sirt)
library(MASS)



aligned.factor.scores <- function(lambda,nu,y){
  #calculate inverse matrix
  lambda1 <- ginv((lambda))
  #create matrix for nu
  ns <- nrow(y)
  nus <- matrix(nu,nrow=ns, ncol=length(nu), byrow=T)
  # y - nu
  y_nu <- y - nu
  # f = inv(lambda)*(y-nu)
  F <- lambda1 %*% t(as.matrix(y_nu))
}

#extract compliance_related variable names
# new model: without compliance_3 and 7

items.PSUP <- c('perceived_stress_sca_1','perceived_stress_sca_2',
                'perceived_stress_sca_3','perceived_stress_sca_4',
                'perceived_stress_sca_5','perceived_stress_sca_6',
                'perceived_stress_sca_7','perceived_stress_sca_8',
                'perceived_stress_sca_9','perceived_stress_sca_10')
cfa.model.PSUP <-'RES =~ perceived_stress_sca_1 + perceived_stress_sca_2+
  perceived_stress_sca_3 + perceived_stress_sca_4 + perceived_stress_sca_5+
  perceived_stress_sca_6 + perceived_stress_sca_7 + perceived_stress_sca_8+
  perceived_stress_sca_9 + perceived_stress_sca_10'


#Load the cleaned cvs file

#load data
data<-read.csv('../../Final_COVIDiSTRESS_Vol2_cleaned.csv')

# reverse
data$perceived_stress_sca_4 <- 4- data$perceived_stress_sca_4
data$perceived_stress_sca_5 <- 4- data$perceived_stress_sca_5
data$perceived_stress_sca_7 <- 4- data$perceived_stress_sca_7
data$perceived_stress_sca_8 <- 4- data$perceived_stress_sca_8

#extract languages with n>=100
n.langs <- table(data$UserLanguage)
list.langs <-labels(n.langs)[[1]]
langs.include <-list.langs[n.langs>=100]
n.include <- n.langs[n.langs>=100]

#extract data
for (i in 1:length(langs.include)){
  if (i == 1){
    data.mi <-data[data$UserLanguage == langs.include[i],]
 }else{
   current <- data[data$UserLanguage == langs.include[i],]
   data.mi <- rbind(data.mi,current)
 }
}

#set and examine fitmeasures
fits <- c('rmsea.scaled','srmr','cfi.scaled','tli.scaled')


#####
#resilience

#general CFA: resilience
cfa.whole.sps <-cfa(model=cfa.model.PSUP, data=data.mi,estimator='WLSMV', group=
                      'UserLanguage')
fitMeasures(cfa.whole.sps)[fits]
# 0.09464355   0.06379443   0.87413322   0.83817129 mediocre

# then metric
cfa.metric.sps <-cfa(model=cfa.model.PSUP, data=data.mi,estimator='WLSMV', group=
                      'UserLanguage',group.equal ='loadings')
fitMeasures(cfa.metric.sps)[fits]
#   0.07243391   0.08021644   0.90799462   0.90521114  acceptable

# scalar
cfa.scalar.sps <-cfa(model=cfa.model.PSUP, data=data.mi,estimator='WLSMV', group=
                       'UserLanguage',group.equal =c('loadings','intercepts'))
fitMeasures(cfa.scalar.sps)[fits]
# 0.1210491    0.1174040    0.6919933    0.7352739  not good

# extract parameters
par.sps <- invariance_alignment_cfa_config(dat = data.mi[,items.PSUP], 
                                            group = data.mi$UserLanguage)
# do alignment
mod1.sps <- invariance.alignment(lambda = par.sps$lambda, nu =
                                   par.sps$nu, align.scale = c(0.2, 0.4), align.pow = c(0.25, 0.25))
# test performance
mod1.sps$es.invariance['R2',]
#  loadings intercepts 
#  0.9850629  0.9930655   
# all ≥ 75%. well addressed.

# monte carlo. IRT test (optional)
#--- find parameter constraints for prespecified tolerance
cmod1.sps <- sirt::invariance_alignment_constraints(model=mod1.sps, nu_parm_tol=.4,
                                                lambda_parm_tol=.4 )

par.cle <- par.sps

# simulation function
simulation_CLE <- function(times,n,data,n.include,seed=1){
  # get data
  data.mi <- data
  items.PSUP <- c('perceived_stress_sca_1','perceived_stress_sca_2',
                  'perceived_stress_sca_3','perceived_stress_sca_4',
                  'perceived_stress_sca_5','perceived_stress_sca_6',
                  'perceived_stress_sca_7','perceived_stress_sca_8',
                  'perceived_stress_sca_9','perceived_stress_sca_10')
  compliance <- items.PSUP
  
  # list for return
  cor.mean <- rep(0,times)
  cor.var <- rep(0,times)
  R2.loading <- rep(0,times)
  R2.intercept <- rep(0,times)
  
  # simulation replication
  # repeat
  for (j in 1:times){
    set.seed(seed)
    G <- n.include # groups
    I <- 10 # items
    
    # lambda, nu, and error_var
    err_var.cle <- matrix(1, nrow=G,ncol=I)
    # simulate data
    # enter group mu and sigma
    data.mi$COMP <- rowMeans(data.mi[,compliance])
    mu<-scale(aggregate(x=data.mi$COMP,
              by = list(data.mi$UserLanguage),
              FUN=mean, na.rm=T)[,2])[,1]
    sigma <- (aggregate(x=data.mi$COMP,
                       by = list(data.mi$UserLanguage),
                       FUN=sd, na.rm=T)[,2])
    N <- rep(n,G)
    dat <- invariance_alignment_simulate(
      par.cle$nu,par.cle$lambda,err_var.cle,mu,sigma,N
    )
    par.simul <- invariance_alignment_cfa_config(dat = dat[,compliance], 
                                               group = dat$group,
                                               estimator = 'WLSMV')
    
    # paramester estimation
    #cfa.test <-cfa(cfa.model.simul,dat,estimator='WLSMV',group='group')
    i#pars <- parameterEstimates(cfa.test)
    
    mod1.simul <- invariance.alignment(lambda = par.simul$lambda, nu =
                                         par.simul$nu, align.scale = c(0.2, 0.4), align.pow = c(0.25, 0.25),
                                       optimizer='nlminb')
    
    # true vs aligned scores
    cfa.model.PSUP <-'PSUP =~ perceived_stress_sca_1 + perceived_stress_sca_2+
  perceived_stress_sca_3 + perceived_stress_sca_4 + perceived_stress_sca_5+
  perceived_stress_sca_6 + perceived_stress_sca_7 + perceived_stress_sca_8+
  perceived_stress_sca_9 + perceived_stress_sca_10'
    cfa.model.simul <- cfa.model.PSUP
    cfa.simul <- cfa(cfa.model.simul,dat,estimator='WLSMV',group='group',
                     group.equal=c('loadings','intercepts'),meanstructure=T)
    
    # get group mean
    params.simul <- parameterEstimates(cfa.simul)
    alpha.simul <- params.simul[(params.simul$op=='~1')&(params.simul$lhs=='PSUP'),'est']
    
    # group mean correlation (Muthen 2018)
    correlation <- corr.test(alpha.simul,mod1.simul$pars$alpha0,method='spearman')$r
    
    
    
    psi.simul <- params.simul[(params.simul$op=='~~')&(params.simul$lhs=='PSUP'),'est']
    correlation.psi <- corr.test(psi.simul,mod1.simul$pars$psi0,method='spearman')$r
    
    cor.mean[j] <- correlation
    cor.var[j] <- correlation.psi
    
    R2.loading[j] <- mod1.simul$es.invariance['R2',1]
    R2.intercept[j] <- mod1.simul$es.invariance['R2',2]
  }
  # make matrix
  to.return <-cbind(cor.mean,cor.var,R2.loading,R2.intercept)
  to.return <- data.matrix(to.return)
  
  message  (sprintf('%d/%d Done',j,times))
  
  return(to.return)
  
}



# use five cores

library(foreach)
library(parallel)
library(doParallel)

cores <- 5
times <- 500

cl <- parallel::makeCluster(cores,type='FORK')
doParallel::registerDoParallel(cl)

# n = 100

start_100 <-Sys.time()
now <- foreach (i = seq(1,times)) %dopar%{
  # get result
  simulation_CLE(1,150,data.mi,28,i)
  #  message(sprintf('%d',i))
}
end_100<-Sys.time()
elapsed_100 <- end_100 - start_100
# merge result
for (i in 1:times){
  if (i == 1){
    simulate_100 <- now[[1]]
  }else{
    simulate_100 <- rbind(simulate_100,now[[i]])
  }
}
# save n = 100
write.csv(data.frame(simulate_100),file='simulate_100.csv',row.names = F)


# n=200
# multicore processing for n=200
start_200 <-Sys.time()
now <- foreach (i = seq(1,times)) %dopar%{
  # get result
  simulation_CLE(1,200,data.mi,28,i)
#  message(sprintf('%d',i))
}
end_200<-Sys.time()
elapsed_200 <- end_200 - start_200
# merge result
for (i in 1:times){
  if (i == 1){
    simulate_200 <- now[[1]]
  }else{
    simulate_200 <- rbind(simulate_200,now[[i]])
  }
}
# save n = 200
write.csv(data.frame(simulate_200),file='simulate_200.csv',row.names = F)

# n= 500
start_500 <-Sys.time()
now <- foreach (i = seq(1,times)) %dopar%{
  # get result
  simulation_CLE(1,500,data.mi,28,i)
  #  message(sprintf('%d',i))
}
end_500<-Sys.time()
elapsed_500 <- end_500 - start_500
# merge result
for (i in 1:times){
  if (i == 1){
    simulate_500 <- now[[1]]
  }else{
    simulate_500 <- rbind(simulate_500,now[[i]])
  }
}
# save n = 200
write.csv(data.frame(simulate_500),file='simulate_500.csv',row.names = F)
# stop cluster
parallel::stopCluster(cl)