
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

library(nloptr)

likfn<-function(pars,model,estimate=TRUE){
  
    model$transitionMatrix[upper.tri(diag(3),TRUE)]<-pars[1:6]
    for(i in 1:model$numberOfChannels){
      emissionArray[,1:model$numberOfSymbols[i],i]<-pars[6+(i-1)*6+1:6]
    }
    model$initialProbs[]<-pars[6+3*6+1:3]
  
  
  if(estimate){
    - sum(logLikMCHMM(model$transitionMatrix, emissionArray, model$initialProbs, obsArray))   
  } else {
      for(i in 1:model$numberOfChannels){
        model$emissionMatrix[[i]][]<-emissionArray[,1:model$numberOfSymbols[i],i]
      }
    model
  }        
}  

eval_eq <- function(pars,model,estimate=TRUE){
  d<-diag(3)
  d[upper.tri(d,TRUE)]<-pars[1:6]
  c(rowSums(d)-1, 
    rowSums(matrix(pars[6+1:6],3,2))-1,
    rowSums(matrix(pars[6+6+1:6],3,2))-1,
    rowSums(matrix(pars[6+12+1:6],3,2))-1,
    sum(pars[6+3*6+1:3])-1)
}
initialvalues <- c(model$transitionMatrix[upper.tri(diag(3),TRUE)],unlist(model$emissionMatrix),model$initialProbs)

fit <- nloptr(initialvalues, likfn, lb = rep(0,27), ub=rep(1,27), eval_g_eq = eval_eq, 
  model=model,estimate=TRUE,opts=list(print_level=3, algorithm="NLOPT_GN_ISRES",maxeval=10000000))

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

fit <- nloptr(initialvalues, likfn, gradfn,lb = rep(-500,14), ub=rep(50,14),
  model=model,estimate=TRUE,opts=list(print_level=1, algorithm="NLOPT_GD_MLSL",maxeval=10000,
    local_opts=list(algorithm="NLOPT_LD_LBFGS",xtol_rel = 1e-4)))
fit#5407.97110558737 

fit <- nloptr(initialvalues, likfn, gradfn,
  model=model,estimate=TRUE,opts=list(print_level=1, algorithm="NLOPT_LD_LBFGS",maxeval=1000,
    xtol_rel = 1e-8)) #5493.270657


fit <- nloptr(initialvalues, likfn, lb = rep(-500,14), ub=rep(10,14),
  model=model,estimate=TRUE,opts=list(print_level=1, algorithm="NLOPT_GN_ISRES",maxeval=1000,
    xtol_rel = 1e-8)) #5493.270657

fit
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