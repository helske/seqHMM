################## Testing ssp ######################

library(seqHMM)

### Default settings
# Single-channel stslist
ssplot(hmm_biofam$observations[[1]])
# Multiple stslists
ssplot(hmm_biofam$observations)
# HMM
ssplot(hmm_biofam)

### type = "I"
# Single-channel stslist
ssplot(hmm_biofam$observations[[1]], type = "I")
# Multiple stslists
ssplot(hmm_biofam$observations, type = "I")
# HMM
ssplot(hmm_biofam, type = "I")


######################### Legends ###############################

### bottom
# Single-channel stslist
ssplot(hmm_biofam$observations[[1]], with.legend = "bottom")
# Multiple stslists
ssplot(hmm_biofam$observations, with.legend = "bottom")
# HMM
ssplot(hmm_biofam, with.legend = "bottom")

### right.combined
# Single-channel stslist
ssplot(hmm_biofam$observations[[1]], with.legend = "right.combined")
# Multiple stslists
ssplot(hmm_biofam$observations, with.legend = "right.combined")
# HMM
ssplot(hmm_biofam, with.legend = "right.combined")

### bottom.combined
# Single-channel stslist
ssplot(hmm_biofam$observations[[1]], with.legend = "right.combined")
# Multiple stslists
ssplot(hmm_biofam$observations, with.legend = "bottom.combined")
# HMM
ssplot(hmm_biofam, with.legend = "bottom.combined")


### ncol.legend

# Single-channel stslist
ssplot(hmm_biofam$observations[[1]], ncol.legend = 2)
# Multiple stslists
ssplot(hmm_biofam$observations, ncol.legend = 2)
# HMM
ssplot(hmm_biofam, ncol.legend = 2)

## bottom
# Single-channel stslist
ssplot(hmm_biofam$observations[[1]], with.legend = "bottom", ncol.legend = 2)
# Multiple stslists
ssplot(hmm_biofam$observations, with.legend = "bottom", ncol.legend = 2)
# HMM
ssplot(hmm_biofam, with.legend = "bottom", ncol.legend = 2)

## right.combined
# Multiple stslists
ssplot(hmm_biofam$observations, with.legend = "right.combined", ncol.legend = 2)
# HMM
ssplot(hmm_biofam, with.legend = "right.combined", ncol.legend = 2)
# HMM + hidden paths
ssplot(hmm_biofam, plots = "both", with.legend = "right.combined", ncol.legend = 2)

## bottom.combined
# Multiple stslists
ssplot(hmm_biofam$observations, with.legend = "bottom.combined", ncol.legend = 2)
# HMM
ssplot(hmm_biofam, with.legend = "bottom.combined", ncol.legend = 2)
# HMM + hidden paths
ssplot(hmm_biofam, plots = "both", with.legend = "bottom.combined", ncol.legend = 2)

######################### Plots ############################

### HMM
# plots = "both"
ssplot(hmm_biofam, plots = "both")
# plots = "hidden.paths"
ssplot(hmm_biofam, plots = "hidden.paths")


######################### Sorting ############################

### from.start
# Single-channel stslist
ssplot(hmm_biofam$observations[[1]], type = "I", sortv = "from.start")
# Multiple stslists
ssplot(hmm_biofam$observations, type = "I", sortv = "from.start")
# HMM
ssplot(hmm_biofam, type = "I", sortv = "from.start")

### from.start, hidden paths
# Single-channel stslist (expect error)
ssplot(hmm_biofam$observations[[1]], type = "I", sortv = "from.start", sort.channel = 0)
# Multiple stslists (expect error)
ssplot(hmm_biofam$observations, type = "I", sortv = "from.start", sort.channel = 0)
# HMM
ssplot(hmm_biofam, type = "I", sortv = "from.start", sort.channel = 0)

# Hidden paths provided
hp <- hidden_paths(hmm_biofam)
# Single-channel stslist
ssplot(hmm_biofam$observations[[1]], type = "I", sortv = "from.start", sort.channel = 0,
       hidden.paths = hp)
# Multiple stslists
ssplot(hmm_biofam$observations, type = "I", sortv = "from.start", sort.channel = 0,
       hidden.paths = hp)
# HMM
ssplot(hmm_biofam, type = "I", sortv = "from.start", sort.channel = 0,
       hidden.paths = hp)


### from.start, false channel
# Single-channel stslist (expect error)
ssplot(hmm_biofam$observations[[1]], type = "I", sortv = "from.start", sort.channel = 4)
# Multiple stslists (expect error)
ssplot(hmm_biofam$observations, type = "I", sortv = "from.start", sort.channel = 4)
# HMM (expect error)
ssplot(hmm_biofam, type = "I", sortv = "from.start", sort.channel = 4)



### from.end
# Single-channel stslist
ssplot(hmm_biofam$observations[[1]], type = "I", sortv = "from.end")
# Multiple stslists
ssplot(hmm_biofam$observations, type = "I", sortv = "from.end")
# HMM
ssplot(hmm_biofam, type = "I", sortv = "from.end")

### from.end, hidden paths

# Hidden paths not provided
# Single-channel stslist (expect error)
ssplot(hmm_biofam$observations[[1]], type = "I", sortv = "from.end", sort.channel = 0)
# Multiple stslists (expect error)
ssplot(hmm_biofam$observations, type = "I", sortv = "from.end", sort.channel = 0)
# HMM
ssplot(hmm_biofam, type = "I", sortv = "from.end", sort.channel = 0)

# Hidden paths provided
hp <- hidden_paths(hmm_biofam)
# Single-channel stslist
ssplot(hmm_biofam$observations[[1]], type = "I", sortv = "from.end", sort.channel = 0,
       hidden.paths = hp)
# Multiple stslists
ssplot(hmm_biofam$observations, type = "I", sortv = "from.end", sort.channel = 0,
       hidden.paths = hp)
# HMM
ssplot(hmm_biofam, type = "I", sortv = "from.end", sort.channel = 0,
       hidden.paths = hp)


### from.end, false channel
# Single-channel stslist (expect error)
ssplot(hmm_biofam$observations[[1]], type = "I", sortv = "from.end", sort.channel = 4)
# Multiple stslists (expect error)
ssplot(hmm_biofam$observations, type = "I", sortv = "from.end", sort.channel = 4)
# HMM (expect error)
ssplot(hmm_biofam, type = "I", sortv = "from.end", sort.channel = 4)


### mds.obs
# Single-channel stslist
ssplot(hmm_biofam$observations[[1]], type = "I", sortv = "mds.obs")
# Multiple stslists
ssplot(hmm_biofam$observations, type = "I", sortv = "mds.obs")
# HMM
ssplot(hmm_biofam, type = "I", sortv = "mds.obs")
ssplot(hmm_biofam, type = "I", sortv = "mds.obs", plots = "both")
ssplot(hmm_biofam, type = "I", sortv = "mds.obs", plots = "hidden.paths")


### mds.hidden
# Hidden paths not provided
# Single-channel stslist (expect error)
ssplot(hmm_biofam$observations[[1]], type = "I", sortv = "mds.hidden")
# Multiple stslists (expect error)
ssplot(hmm_biofam$observations, type = "I", sortv = "mds.hidden")
#####
# Hidden paths provided
hp <- hidden_paths(hmm_biofam)
# Single-channel stslist
ssplot(hmm_biofam$observations[[1]], type = "I", sortv = "mds.hidden", hidden.paths = hp)
# Multiple stslists
ssplot(hmm_biofam$observations, type = "I", sortv = "mds.hidden", hidden.paths = hp)
# HMM
ssplot(hmm_biofam, type = "I", sortv = "mds.hidden")
ssplot(hmm_biofam, type = "I", sortv = "mds.hidden", plots = "both")
ssplot(hmm_biofam, type = "I", sortv = "mds.hidden", plots = "hidden.paths")


### vector
# Single-channel stslist
ssplot(hmm_biofam$observations[[1]], type = "I", sortv = c(2000:1))
# Multiple stslists
ssplot(hmm_biofam$observations, type = "I", sortv = c(2000:1))
# HMM
ssplot(hmm_biofam, type = "I", sortv = c(2000:1))
ssplot(hmm_biofam, type = "I", sortv = c(2000:1), plots = "both")
ssplot(hmm_biofam, type = "I", sortv = c(2000:1), plots = "hidden.paths")



######################### Plot axes ############################

### No x-axis
# Single-channel stslist
ssplot(hmm_biofam$observations[[1]], xaxis = FALSE)
# Multiple stslists
ssplot(hmm_biofam$observations, xaxis = FALSE)
# HMM
ssplot(hmm_biofam, xaxis = FALSE)
ssplot(hmm_biofam, xaxis = FALSE, plots = "both")
ssplot(hmm_biofam, xaxis = FALSE, plots = "hidden.paths")

### With y-axis
# Single-channel stslist
ssplot(hmm_biofam$observations[[1]], yaxis = TRUE)
# Multiple stslists
ssplot(hmm_biofam$observations, yaxis = TRUE)
# HMM
ssplot(hmm_biofam, yaxis = TRUE)
ssplot(hmm_biofam, yaxis = TRUE, plots = "both")
ssplot(hmm_biofam, yaxis = TRUE, plots = "hidden.paths")


######################### Missing states in legends ############################


### Default settings
# Single-channel stslist
ssplot(hmm_biofam$observations[[1]], with.missing.legend = TRUE)
# Multiple stslists
ssplot(hmm_biofam$observations, with.missing.legend = TRUE)
# HMM
ssplot(hmm_biofam, with.missing.legend = TRUE)


### With missing
# Single-channel stslist
ssplot(hmm_biofam$observations[[1]], with.missing.legend = TRUE, with.missing = TRUE)
# Multiple stslists
ssplot(hmm_biofam$observations, with.missing.legend = TRUE, with.missing = TRUE)
# HMM
ssplot(hmm_biofam, with.missing.legend = TRUE, with.missing = TRUE)

