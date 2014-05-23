L1<-function(models){
  
  K<-length(models) # klusterien lkm
  
  # Laskee L1-"uskottavuuden" klusteriratkaisulle kaikista klusterin malleista 
  # käyttäen WITMSE-paperin kaavaa (5)
  
  # Klusterin todennäköisyys (empiirinen, sekvenssien lkm / kaikkien sekvenssien lukumäärällä
  clusterSizes<-sapply(models,"[[","numberOfSequences")
  pLambda<-clusterSizes/sum(clusterSizes)
  
  
  L<-numeric(K)
  
  for(i in 1:K){
    
    n<-models[[i]]$lengthOfSequences
    k<-models[[i]]$numberOfSequences
    
    #polkujen todennäköisyydet (yksi jokaiselle sekvenssille)
    viterbi<-mostProbablePath(models[[i]])
    
    #pStates<-viterbi$logP
    
    #sekvenssien ja piilotilojen todennäköisyydet ehdolla polut
    pSeq<-pStates<-numeric(k)    
    names(models[[i]]$initialProbs)<-rownames(models[[i]]$transitionMatrix)
    for(s in 1:models[[i]]$numberOfChannels)
      models[[i]]$observations[[s]]<-sapply(models[[i]]$observations[[s]],as.character)
    for(j in 1:k){
      pStates[j]<-log(models[[i]]$initialProbs[viterbi$mpp[j,1]])+
        sum(log(models[[i]]$trans[cbind(viterbi$mpp[j,-n],viterbi$mpp[j,-1])]))
                                                              
      for(s in 1:models[[i]]$numberOfChannels)
        pSeq[j]<-pSeq[j]+sum(log(
          models[[i]]$emissionMatrix[[s]][cbind(viterbi$mpp[j,],
                                                models[[i]]$observations[[s]][j,])]),na.rm=TRUE)
  
    }
   
    L[i]<-k*log(pLambda[i])+sum(pStates)+sum(pSeq)
  }
  L  
}
