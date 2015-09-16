data(biofam, package = "TraMineR")
biofam3c <- with(biofam, {
## Building one channel per type of event left, children or married
bf <- as.matrix(biofam[, 10:25])
children <-  bf == 4 | bf == 5 | bf == 6
married <- bf == 2 | bf == 3 | bf == 6
left <- bf == 1 | bf == 3 | bf == 5 | bf == 6 | bf == 7

children[children == TRUE] <- "children"
children[children == FALSE] <- "childless"
# Divorced parents
div <- bf[(rowSums(bf == 7)>0 & rowSums(bf == 5) > 0) |
            (rowSums(bf == 7)>0 & rowSums(bf == 6) > 0),]
children[rownames(bf) %in% rownames(div) & bf == 7] <- "children"

married[married == TRUE] <- "married"
married[married == FALSE] <- "single"
married[bf == 7] <- "divorced"

left[left == TRUE] <- "left home"
left[left == FALSE] <- "with parents"
# Divorced living with parents (before divorce)
wp <- bf[(rowSums(bf == 7) > 0 & rowSums(bf == 2) > 0 & rowSums(bf == 3) == 0 &
            rowSums(bf == 5) == 0 & rowSums(bf == 6) == 0) |
           (rowSums(bf == 7) > 0 & rowSums(bf == 4) > 0 & rowSums(bf == 3) == 0 &
              rowSums(bf == 5) == 0 & rowSums(bf == 6) == 0),]
left[rownames(bf) %in% rownames(wp) & bf == 7] <- "with parents"

list("children" = children, "married" = married, 
     "left" = left, "covariates" = biofam[, c(1:9, 26:27)])
})
biofam3c