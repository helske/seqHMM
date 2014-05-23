#' Transform Multichannel Hidden Markov Model to Single Channel Representation
#' 
#' @export
#' @param model Object of class HMModel.
MCtoSC<-function(model){
  
  if(model$numberOfChannels==1)
    return(model)
    
  B<-matrix(0,model$numberOfStates,prod(model$numberOfSymbols))
  
  colnames(B)<-apply(
    expand.grid(colnames(model$emissionMatrix[[1]]),
                colnames(model$emissionMatrix[[2]]),
                colnames(model$emissionMatrix[[3]])),
    1,paste0,collapse="/")
  rownames(B)<-rownames(model$emissionMatrix[[1]])
  for(i in 1:model$numberOfStates){
    B[i,]<-apply(expand.grid(model$emissionMatrix[[1]][i,],
                             model$emissionMatrix[[2]][i,],
                             model$emissionMatrix[[3]][i,]),1,prod)
    
    
  }
  
  
  
  modelx<-model
  modelx$emissionMatrix<-B
  modelx$numberOfSymbols<-ncol(B)
  modelx$numberOfChannels<-as.integer(1)
  modelx$symbolNames<-colnames(B)
  modelx$channelNames<-NULL
  modelx$observations<-
    mapply(paste,model$obs[[1]],model$obs[[2]],
           model$obs[[3]],MoreArgs=list(sep="/"))
  modelx$observations[is.na(model$obs[[1]]) | is.na(model$obs[[2]]) | 
                        is.na(model$obs[[3]])]<-NA
  modelx
}

# Jos logLik(modelx)!=logLik(model), pit?isi johtua puuttuvista
# (toinen sallii osittain puuttuvan tiedon):
# modelz<-model
# modelz$observations[[1]][is.na(fit$model$obs[[1]]) | 
#                           is.na(fit$model$obs[[2]]) | 
#                           is.na(fit$model$obs[[3]])]<-NA
# modelz$observations[[2]][is.na(fit$model$obs[[1]]) | 
#                           is.na(fit$model$obs[[2]]) | 
#                           is.na(fit$model$obs[[3]])]<-NA
# modelz$observations[[3]][is.na(fit$model$obs[[1]]) | 
#                           is.na(fit$model$obs[[2]]) | 
#                           is.na(fit$model$obs[[3]])]<-NA
# logLik(modelx)-logLik(modelz)