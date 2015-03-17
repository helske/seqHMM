
library(seqHMM)
set.seed(1)

onlyfair <- sample(0:1, size = 100, replace = TRUE)

obs <- matrix(0,100,100)

obs[onlyfair==1,] <- sample(c("heads","tails"),size=length(obs[onlyfair==1,]),replace=TRUE)
obs[sample(which(obs==0), size=sum(obs==0)/1.5)] <- 
  sample(c("heads","tails"), prob = c(0.9,0.1), size=sum(obs==0)/1.5, replace=TRUE)
obs[which(obs==0)] <- 
  sample(c("heads","tails"), size=sum(obs==0)/2, replace=TRUE)

A1 <- matrix(1,1,1)
B1 <- matrix(c(0.5,0.5),1,2)

A2 <- matrix(c(0.5,0.5,0.5,0.5),2,2)
B2 <- matrix(c(0.5,0.9,0.5,0.1),2,2)

library(TraMineR)
obs <- seqdef(obs)
hmm<-buildMixHMM(obs=obs,init=list(1,c(0.5,0.5)),
                 trans=list(A1,A2),emission=list(B1,B2),X=matrix(onlyfair,ncol=1), beta=matrix(0:1,1,2))

fit<-fitMixHMM(hmm)
