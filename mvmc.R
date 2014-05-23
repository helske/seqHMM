#symbolNames<-lapply(observations,function(x) sort(unique(as.character(unlist(x)))))


load(file="V:/NEPS/Data/bioseq_m.Rda")
load(file="V:/NEPS/Data/parentseq_m.Rda")
load(file="V:/NEPS/Data/partnerseq_m.Rda")

library(TraMineR)

w.labels <- c("School", "VocPrep", "VocTrain", "Emp", "Unemp", "Service", "ParLeave", "Gap")
w.alphabet <- c("SC", "VP", "VT", "EM", "UN", "ML", "PL", "GP")
w.seq <- seqdef(bioseq, var=20:453, alphabet=w.alphabet, labels = w.labels, 
                start=15, right=NA)

p.labels <- c("Single", "Cohabiting", "Married", "Divorced", "Widowed")
p.alphabet <- c("S", "C", "M", "D", "W")
p.seq <- seqdef(partnerseq, var=20:453, alphabet=p.alphabet, labels = p.labels,
                start=15, right=NA)

c.labels <- c("No children", "Children")
c.alphabet <- c("NC", "CH")
c.seq <- seqdef(parentseq, var=22:455, alphabet=c.alphabet, 
                labels = c.labels, start=15, right=NA)






library(RColorBrewer)

attr(w.seq, "cpal") <- brewer.pal(8, "Set1")[8:1]
attr(p.seq, "cpal") <- brewer.pal(12, "Set3")[c(5,7,12,4,10)]
attr(c.seq, "cpal") <- brewer.pal(8, "Set3")[c(1,6)]



# seqdplot(w.seq, border=NA)
# seqdplot(p.seq, border=NA)
# seqdplot(c.seq, border=NA)




p.distHAM <- seqsubm(p.seq, method="CONSTANT", with.missing=TRUE, miss.cost=0)
p.distHAM[1:5,1:5] <- matrix(c(0,2,2,4,4,
                               2,0,1,2,2,
                               2,1,0,2,2,
                               4,2,2,0,1,
                               4,2,2,1,0), byrow=TRUE, nrow=5, ncol=5)

c.distHAM <- seqsubm(c.seq, method="CONSTANT", with.missing=TRUE, miss.cost=0,
                     cval=4)

w.distHAM <- seqsubm(w.seq, method="CONSTANT", with.missing=TRUE, miss.cost=0)
w.distHAM[1:8,1:8] <- matrix(c(0,1,1,4,3,2,2,2,
                               1,0,1,4,3,2,2,2,
                               1,1,0,4,3,2,2,2,
                               4,4,4,0,3,2,2,2,
                               3,3,3,3,0,2,2,1,
                               2,2,2,2,2,0,2,1,
                               2,2,2,2,2,2,0,1,
                               2,2,2,2,1,1,1,0), byrow=TRUE, nrow=8, ncol=8)


dist.u.HAM <- seqdistmc(channels=list(p.seq, c.seq, w.seq), method="HAM",
                        norm=FALSE, sm =list(p.distHAM, c.distHAM, w.distHAM),
                        full.matrix=TRUE, link="mean", cweight=c(1,1,1.5),
                        with.missing=TRUE)

library(cluster)

clusterward.u.HAM <- agnes(dist.u.HAM, diss = TRUE, method = "ward")

plot(clusterward.u.HAM, ask=FALSE, which.plot=2)

clHAM.8 <- cutree(clusterward.u.HAM, k = 8)


load("V:/NEPS/KKsekvenssit/HMMmalleja/malli_cl1_5.rda")

model<-malli_cl1_5[[4]]$model
save(model,file="U:/valgrind/model.rda")


# lap <- lapply(malli_cl1_5, `[[`, "opt") 
# lapply(lap, "[[", "minimum")
# lapply(lap, "[[", "iterations")

fit <- malli_cl1_5[[4]]

HMM.cl1 <- fit$model

initprobs <- round(HMM.cl1$initialProbs,9)

A <- HMM.cl1$transitionMatrix
A <- A/rowSums(A)

B <- HMM.cl1$emission

Blist <- rep(B, each=491)

wsplit <- split(bioseq[clHAM.8==1,20:453], f=1:491)
psplit <- split(partnerseq[clHAM.8==1,20:453], f=1:491)
csplit <- split(parentseq[clHAM.8==1,22:455], f=1:491)

obs <- c(wsplit, psplit, csplit)

symbolNames<-rep(lapply(B,colnames),each=491)
HMMrep <- HMModel(observations=obs,
                  transitionMatrix=A, 
                  emissionMatrix=Blist, 
                  initialProbs=initprobs,symbolNames=symbolNames)

repseq <- mostProbablePath(HMMrep)
mpp.seq <- seqdef(mpp)
seqIplot(repseq)
seqs<-seqdef(mostProbablePath(HMM.cl1))
seqIplot(seqs)

x<-posteriorProbs(HMM.cl1)

