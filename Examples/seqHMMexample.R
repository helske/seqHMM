library(seqHMM)
library(TraMineR)

data(biofam)

# Only the first 500 individuals
biofam <- biofam[1:500,]

## Sequence data for the first six individuals
head(biofam[10:25])

## Building one channel per type of event (left, children or married)
bf <- as.matrix(biofam[, 10:25]) 
children <-  bf==4 | bf==5 | bf==6 
married <- bf == 2 | bf== 3 | bf==6 
left <- bf==1 | bf==3 | bf==5 | bf==6 | bf==7

## Labels
children[children==TRUE] <- "Children" 
children[children==FALSE] <- "Childless"

married[married==TRUE] <- "Married" 
married[married==FALSE] <- "Single"

left[left==TRUE] <- "Left home" 
left[left==FALSE] <- "With parents"

## Building sequence objects (starting at age 15)
child.seq <- seqdef(children, start=15) 
marr.seq <- seqdef(married, start=15) 
left.seq <- seqdef(left, start=15)

## Choosing colours for states
attr(child.seq, "cpal") <- c("#66C2A5", "#FC8D62")
attr(marr.seq, "cpal") <- c("#E7298A", "#E6AB02")
attr(left.seq, "cpal") <- c("#A6CEE3", "#E31A1C")

## Plotting state distribution plots of observations
ssplot(list(child.seq, marr.seq, left.seq), type="d", 
            plots="obs", title="State distribution plots")

## The same with first defining the plot with function ssp 
## and then plotting with function plot
ssp1 <- ssp(list(child.seq, marr.seq, left.seq), type="d", 
                    plots="obs", title="State distribution plots")
plot(ssp1)

## Preparing plots for state distributios and index plots of observations for women
#  Sorting by scores from multidimensional scaling
ssp_f2 <- ssp(list(child.seq[biofam$sex=="woman",],
                           marr.seq[biofam$sex=="woman",], 
                           left.seq[biofam$sex=="woman",]),
                      type="d", plots="obs", border=NA,
                      title="State distributions for women", title.n=FALSE,
                      ylab=c("Children", "Married", "Left home"), 
                      withlegend=FALSE, ylab.pos=c(1,2,1))
ssp_f3 <- ssp(list(child.seq[biofam$sex=="woman",],
                           marr.seq[biofam$sex=="woman",], 
                           left.seq[biofam$sex=="woman",]),
                      type="I", sortv="mds.obs", plots="obs", 
                      title="Sequences for women",
                      ylab=c("Children", "Married", "Left home"), withlegend=FALSE,
                      ylab.pos=c(1.5,2.5,1.5))

## Preparing plots for state distributios and index plots of observations for men
ssp_m2 <- ssp(list(child.seq[biofam$sex=="man",], 
                           marr.seq[biofam$sex=="man",], 
                           left.seq[biofam$sex=="man",]), 
                      type="d", plots="obs", border=NA,
                      title="State distributions for men", title.n=FALSE,
                      ylab=c("Children", "Married", "Left home"), 
                      withlegend=FALSE, ylab.pos=c(1,2,1))
ssp_m3 <- ssp(list(child.seq[biofam$sex=="man",],
                           marr.seq[biofam$sex=="man",], 
                           left.seq[biofam$sex=="man",]),
                      type="I", sortv="mds.obs", plots="obs", 
                      title="Sequences for men",
                      ylab=c("Children", "Married", "Left home"), withlegend=FALSE,
                      ylab.pos=c(1.5,2.5,1.5))

## Plotting state distributions and index plots of observations for women and men in two columns 
gridplot(list(ssp_f2, ssp_f3, ssp_m2, ssp_m3), cols=2, byrow=TRUE, 
           row.prop=c(0.42,0.42,0.16))


# Initial values for emission matrices 
B_child <- matrix(NA, nrow=4, ncol=2) 
B_child[1,] <- seqstatf(child.seq[,1:4])[,2]+0.1 
B_child[2,] <- seqstatf(child.seq[,5:8])[,2]+0.1 
B_child[3,] <- seqstatf(child.seq[,9:12])[,2]+0.1 
B_child[4,] <- seqstatf(child.seq[,13:16])[,2]+0.1 
B_child <- B_child/rowSums(B_child)

B_marr <- matrix(NA, nrow=4, ncol=2) 
B_marr[1,] <- seqstatf(marr.seq[,1:4])[,2]+0.1 
B_marr[2,] <- seqstatf(marr.seq[,5:8])[,2]+0.1 
B_marr[3,] <- seqstatf(marr.seq[,9:12])[,2]+0.1 
B_marr[4,] <- seqstatf(marr.seq[,13:16])[,2]+0.1 
B_marr <- B_marr/rowSums(B_marr)

B_left <- matrix(NA, nrow=4, ncol=2) 
B_left[1,] <- seqstatf(left.seq[,1:4])[,2]+0.1 
B_left[2,] <- seqstatf(left.seq[,5:8])[,2]+0.1 
B_left[3,] <- seqstatf(left.seq[,9:12])[,2]+0.1 
B_left[4,] <- seqstatf(left.seq[,13:16])[,2]+0.1 
B_left <- B_left/rowSums(B_left)

# Initial values for transition matrix 
A <- matrix(c(0.9,   0.06, 0.03, 0.01,
              0,    0.9, 0.07, 0.03, 
              0,      0,  0.9,  0.1, 
              0,      0,    0,    1), 
            nrow=4, ncol=4, byrow=TRUE)

# Initial values for initial state probabilities 
initialProbs <- c(0.9, 0.07, 0.02, 0.01)

## Building hidden Markov model with initial parameter values 
bHMM <- buildHMM(observations=list(child.seq, marr.seq, left.seq), 
                 transitionMatrix=A, 
                 emissionMatrix=list(B_child, B_marr, B_left),
                 initialProbs=initialProbs)

## Fitting hidden Markov model 
HMM <- fitHMM(bHMM, em.control=list(maxit=100,reltol=1e-8), 
              itnmax=10000, method="BFGS")

## Plot HMM
plot(HMM$model)

## Prettier version
plot(HMM$model, 
     # larger vertices 
     vertex.size=50, 
     # thicker edges with varying curvature 
     cex.edge.width=3, edge.curved=c(0,-0.7,0.6,0,-0.7,0),
     # Show only states with emission prob. > 0.1
     combine.slices=0.1, 
     # Label for combined states
     combined.slice.label="States with probability < 0.1",
     # Less space for legend
     legend.prop=0.3)


## Plotting observations and hidden states
ssplot(HMM$model)

## Prettier version
ssplot(HMM$model, type="I",
                plots="both",
                # Sorting subjects according to multidimensional
                # scaling scores of the most probable hidden state paths
                sortv="mds.mpp", 
                # Naming the channels
                ylab=c("Children", "Married", "Left home"), 
                # Title for the plot
                title="Observed sequences and the 
                most probable paths of hidden states",
                # Labels for hidden states (most common states)
                mpp.labels=c("1: Childless single, with parents", 
                             "2: Childless single, left home",
                             "3: Married without children",
                             "4: Married parent, left home"),
                # Colours for hidden states
                mpp.col=c("olivedrab", "bisque", "plum", "indianred"),
                # Labels for x axis
                xtlab=15:30,
                # Proportion for legends
                legend.prop=0.45)

## Likelihood
logLik(HMM$model)

# -4103.938

## BIC
BIC(HMM$model)

# 8591.137

## Trimming HMM
trimmedHMM <- trimHMM(HMM$model, maxit=100, zerotol=1e-02)

## Emission probabilities of the original HMM
HMM$model$emiss

# $`1`
# symbolNames
# stateNames  Childless      Children
# 1 1.00000000 1.454112e-177
# 2 1.00000000  3.134555e-21
# 3 1.00000000  1.597482e-14
# 4 0.01953035  9.804697e-01
# 
# $`2`
# symbolNames
# stateNames      Married     Single
# 1 1.220309e-16 1.00000000
# 2 1.602734e-02 0.98397266
# 3 9.833996e-01 0.01660043
# 4 9.464670e-01 0.05353301
# 
# $`3`
# symbolNames
# stateNames    Left home With parents
# 1 1.037512e-18 1.000000e+00
# 2 1.000000e+00 3.735228e-13
# 3 7.127756e-01 2.872244e-01
# 4 1.000000e+00 3.631100e-43

## Emission probabilities of the trimmed HMM
trimmedHMM$emiss

# $`1`
# symbolNames
# stateNames  Childless  Children
# 1 1.00000000 0.0000000
# 2 1.00000000 0.0000000
# 3 1.00000000 0.0000000
# 4 0.01953053 0.9804695
# 
# $`2`
# symbolNames
# stateNames    Married     Single
# 1 0.00000000 1.00000000
# 2 0.01603142 0.98396858
# 3 0.98340199 0.01659801
# 4 0.94646698 0.05353302
# 
# $`3`
# symbolNames
# stateNames Left home With parents
# 1 0.0000000    1.0000000
# 2 1.0000000    0.0000000
# 3 0.7127736    0.2872264
# 4 1.0000000    0.0000000

## Converting multichannel model to single channel model
scHMM <- MCtoSC(HMM$model)

ssplot(scHMM, sortv="from.end", sort.channel=0, legend.prop=0.45)
