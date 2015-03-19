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

## Defining the plot for state distribution plots of observations
mcsp1 <- defineMCSP(list(child.seq, marr.seq, left.seq), type="d", 
                    plots="obs", title="State distribution plots")

## Plotting mcsp1
plot(mcsp1)

## Preparing plots for state distributios and index plots of observations for women
#  Sorting by scores from multidimensional scaling
mcsp_f2 <- defineMCSP(list(child.seq[biofam$sex=="woman",],
                           marr.seq[biofam$sex=="woman",], 
                           left.seq[biofam$sex=="woman",]),
                      type="d", plots="obs", border=NA,
                      title="State distributions for women", title.n=FALSE,
                      ylab=c("Children", "Married", "Left home"), 
                      withlegend=FALSE, ylab.pos=c(1,2,1))
mcsp_f3 <- defineMCSP(list(child.seq[biofam$sex=="woman",],
                           marr.seq[biofam$sex=="woman",], 
                           left.seq[biofam$sex=="woman",]),
                      type="I", sortv="mds.obs", plots="obs", 
                      title="Sequences for women",
                      ylab=c("Children", "Married", "Left home"), withlegend=FALSE,
                      ylab.pos=c(1.5,2.5,1.5))

## Preparing plots for state distributios and index plots of observations for men
mcsp_m2 <- defineMCSP(list(child.seq[biofam$sex=="man",], 
                           marr.seq[biofam$sex=="man",], 
                           left.seq[biofam$sex=="man",]), 
                      type="d", plots="obs", border=NA,
                      title="State distributions for men", title.n=FALSE,
                      ylab=c("Children", "Married", "Left home"), 
                      withlegend=FALSE, ylab.pos=c(1,2,1))
mcsp_m3 <- defineMCSP(list(child.seq[biofam$sex=="man",],
                           marr.seq[biofam$sex=="man",], 
                           left.seq[biofam$sex=="man",]),
                      type="I", sortv="mds.obs", plots="obs", 
                      title="Sequences for men",
                      ylab=c("Children", "Married", "Left home"), withlegend=FALSE,
                      ylab.pos=c(1.5,2.5,1.5))

## Plotting state distributions and index plots of observations for women and men in two columns 
gridplot(list(mcsp_f2, mcsp_f3, mcsp_m2, mcsp_m3), cols=2, byrow=TRUE, 
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
     # Legend with two columns and less space
     combined.slice.label="States with probability < 0.05")


# Plotting observations and hidden states
plot(defineMCSP(HMM$model))

# Prettier version
plot(defineMCSP(HMM$model, type="I", 
                    plots="both", 
                    # Sorting subjects according to multidimensional
                    # scaling scores of the most probable hidden state paths
                    sortv="mds.mpp", 
                    # Naming the channels
                    ylab=c("Children", "Married", "Left home"), 
                    # Title for the plot, number of sequences removed
                    title="Observed sequences and the 
most probable paths of hidden states",
                    xtlab=15:30))

