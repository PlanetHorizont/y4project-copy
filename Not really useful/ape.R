library(ape)
data(bird.orders)
hc <- as.hclust(bird.orders)
tr <- as.phylo(hc)
all.equal(bird.orders, tr) # TRUE
### shows the three plots for tree objects:
dend <- as.dendrogram(hc)
layout(matrix(c(1:3, 3), 2, 2))
plot(bird.orders, font = 1)
plot(hc)
par(mar = c(8, 0, 0, 0)) # leave space for the labels
plot(dend)
as.phylo.formula 
### how to get identical plots with
### plot.phylo and plot.dendrogram:
layout(matrix(1:2, 2, 1))
plot(bird.orders, font = 1, no.margin = TRUE, label.offset = 0.4)
par(mar = c(0, 0, 0, 8))
plot(dend, horiz = TRUE)
layout(1)


data(carnivora)
plot(as.phylo(~SuperFamily/Family/Genus/Species, data=carnivora))


tr <- rtree(30)
ch <- rcoal(30)
plot(ch)
axisPhylo()
plot(tr, "c", FALSE, direction = "u")
axisPhylo(2, las = 1)

data(woodmouse)
base.freq(woodmouse)
base.freq(woodmouse, TRUE)
base.freq(woodmouse, TRUE, TRUE)
GC.content(woodmouse)
Ftab(woodmouse)
Ftab(woodmouse[1, ], woodmouse[2, ]) # same than above
Ftab(woodmouse[14:15, ]) # between the last two

# Generate random phylogeny
n <- 100
phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)
# Generate random data and standardize to have mean 0 and variance 1
X1 <- rTraitCont(phy, model = "BM", sigma = 1)
X1 <- (X1 - mean(X1))/var(X1)
# Simulate binary Y
sim.dat <- data.frame(Y=array(0, dim=n), X1=X1, row.names=phy$tip.label)
36 binaryPGLMM
sim.dat$Y <- binaryPGLMM.sim(Y ~ X1, phy=phy, data=sim.dat, s2=.5,
                             B=matrix(c(0,.25),nrow=2,ncol=1), nrep=1)$Y
sim.dat



tr <- list(edge = matrix(c(3, 1,3,2), 2, 2), tip.label = c("a","b"), Nnode = 1L)
class(tr) <- "phylo"
str(tr)
plot(tr)

tr<-rtree(10)
str(tr)
tr$edge
plot(tr)
Ntip(tr)
Nnode(tr)
Nedge(tr)

which(tr$edge == 5)
tr$edge


start <- which(tr$edge[, 1] == length(tr$tip.label) + 1)
end <- c(start[-1] - 1, dim(tr$edge)[1])
start
end
tr$edge[start, 2]

data(bird.families)
n <- length(bird.families$tip.label)
start <- which(bird.families$edge[, 1] == n + 1)
end <- c(start[-1] - 1, dim(bird.families$edge)[1])
start
end
bird.families$edge[start, 2]
plot(bird.families)
nodelabels()
startB <- which(tr$edge[, 1] == tr$edge[start[1], 2])
endB <- c(startB[-1] - 1, end[1])
startB
endB



tree<-rcoal(10)
plot(tree)
str(tree)
