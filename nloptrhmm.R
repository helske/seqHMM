gradfn<-function(pars,model,estimate){
  
  if(any(!is.finite(exp(pars))))
    return(.Machine$double.xmax^075)
  if(npTM>0){
    model$transitionMatrix[maxTM]<-maxTMvalue     
    model$transitionMatrix[paramTM]<-exp(pars[1:npTM])
    rowSumsA<-rowSums(model$transitionMatrix)
    model$transitionMatrix<-model$transitionMatrix/rowSumsA         
  }
  if(sum(npEM)>0){            
    for(i in 1:model$numberOfChannels){
      emissionArray[,1:model$numberOfSymbols[i],i][maxEM[[i]]]<-maxEMvalue[[i]]    
      emissionArray[,1:model$numberOfSymbols[i],i][paramEM[[i]]]<-
        exp(pars[(npTM+1+c(0,cumsum(npEM))[i]):(npTM+cumsum(npEM)[i])])      
      rowSumsB[,i]<-rowSums(emissionArray[,1:model$numberOfSymbols[i],i, drop = FALSE])
      emissionArray[,1:model$numberOfSymbols[i],i]<-
        emissionArray[,1:model$numberOfSymbols[i],i]/rowSumsB[,i]
    }
  }
  
  if(npIP>0){
    model$initialProbs[maxIP]<-maxIPvalue
    model$initialProbs[paramIP]<-exp(pars[(npTM+sum(npEM)+1):(npTM+sum(npEM)+npIP)])
    sumInit<-sum(model$initialProbs)
    model$initialProbs[]<-model$initialProbs/sumInit
  } 
  
  - gradientMC(model$transitionMatrix, emissionArray, model$initialProbs, obsArray,
    rowSumsA,rowSumsB,sumInit,transNZ,emissNZ,initNZ,exp(pars))
  
}

library(nloptr)
fit1 <- nloptr(initialvalues, likfn, gradfn,lb = rep(-10,14), ub=rep(5,14),
  model=model,estimate=TRUE,opts=list(ranseed = 1,print_level=1, algorithm="NLOPT_GD_MLSL",maxeval=10000,
    local_opts=list(algorithm="NLOPT_LD_LBFGS",xtol_rel = 1e-4)))
fit1#5407.97110558737 


HMM2 <- fitHMM(bHMM, use.em=FALSE)
HMM2$logLik #-5508.148
HMM2 <- fitHMM(bHMM, use.optimx = FALSE)
HMM2$logLik #-5507.003
HMM2 <- fitHMM(bHMM)
HMM2$logLik #-5507.003
HMM2 <- fitHMM(bHMM, use.em=FALSE,method="Nelder-Mead")
HMM2$logLik #-5512.487
HMM2 <- fitHMM(bHMM, use.em=FALSE,optimx.control = list(fnscale=1))
HMM2$logLik #-5476.575 !!!!!!!

HMM2 <- fitHMM(bHMM, use.em=FALSE,soft = FALSE)

likfn<-function(pars,model,estimate=TRUE){
  if((sum(pars[1:2])>1) || any(pars<0) || any(pars>1)) return(Inf)
    model$transitionMatrix[upper.tri(diag(3))]<-pars[1:3]
    diag(model$transitionMatrix) <- 1 - (rowSums(model$transitionMatrix) - diag(model$transitionMatrix))
    
    for(i in 1:model$numberOfChannels){
      emissionArray[,2,i]<-pars[3+(i-1)*3+1:3]
      emissionArray[,1,i]<-1-pars[3+(i-1)*3+1:3]
    }
    model$initialProbs[1:2]<-pars[3+3*3+1:2]
    model$initialProbs[3] <- 1 - sum(model$initialProbs[1:2])
  
  
  if(estimate){
    - sum(logLikMCHMM(model$transitionMatrix, emissionArray, model$initialProbs, obsArray))   
  } else {
      for(i in 1:model$numberOfChannels){
        model$emissionMatrix[[i]][]<-emissionArray[,1:model$numberOfSymbols[i],i]
      }
    model
  }        
}  

eval_g_ineq <- function(pars,model,estimate=TRUE){
  c(sum(pars[1:2])-1)
}
initialvalues <- c(model$transitionMatrix[upper.tri(diag(3),FALSE)],
  model$emissionMatrix[[1]][,2],
  model$emissionMatrix[[2]][,2],
  model$emissionMatrix[[3]][,2],model$initialProbs[1:2])

fit <- nloptr(initialvalues, likfn, lb = rep(0,14), ub=rep(1,14),# eval_g_ineq = eval_g_ineq, 
  model=model,estimate=TRUE,opts=list(print_level=1, algorithm="NLOPT_GN_MLSL",maxeval=10000,
    local_opts=list(algorithm="NLOPT_LN_BOBYQA",xtol_rel = 1e-4)))

fit



fit <- nloptr(initialvalues, likfn, gradfn,lb = rep(-10,14), ub=rep(10,14),
  model=model,estimate=TRUE,opts=list(print_level=3, algorithm="NLOPT_GD_STOGO",maxeval=10000))


initialvalues2 <- initialvalues
initialvalues2[initialvalues2< -5] <- -5
initialvalues2[initialvalues2> 5] <- 5
fit1 <- nloptr(initialvalues, likfn, gradfn,lb = rep(-500,14), ub=rep(5,14),
  model=model,estimate=TRUE,opts=list(print_level=1, algorithm="NLOPT_GD_MLSL",maxeval=10000,
    local_opts=list(algorithm="NLOPT_LD_LBFGS",xtol_rel = 1e-4)))
fit1#5407.97110558737 
best<-likfn(fit1$sol,model,F)

fit <- nloptr(initialvalues, eval_f=likfn, lb = rep(-500,14), ub=rep(5,14),
  model=model,estimate=TRUE,opts=list(print_level=1, algorithm="NLOPT_GN_MLSL",maxeval=10000,
    local_opts=list(algorithm="NLOPT_LN_NELDERMEAD",xtol_rel = 1e-4)))
fit#5407.97110558737 

fit <- nloptr(initialvalues, likfn, gradfn,
  model=model,estimate=TRUE,opts=list(print_level=1, algorithm="NLOPT_LD_LBFGS",maxeval=1000,
    xtol_rel = 1e-8)) #5493.270657


fit <- nloptr(initialvalues, likfn, lb = rep(-500,14), ub=rep(10,14),
  model=model,estimate=TRUE,opts=list(print_level=1, algorithm="NLOPT_GN_ISRES",maxeval=1000,
    xtol_rel = 1e-8)) #5493.270657

fit
