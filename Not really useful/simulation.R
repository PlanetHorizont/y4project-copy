
tree <- simulate.tree(method = "hudson",sample = 20,current = 100,ancestral = 1000,time = 50)

## This plots (on the screen) the coalescent tree, with labels

plot.tree(tree,labels = TRUE)




 Pca <- function(n, N = 1000) {

      # This function will approximate the probability of no
      # coalescent events in a sample of n chromosomes from
      # a diploid population of size N

      1.0 - (n * (n-1)) / (4 * N)

      }

   PcaExact <- function(n, N = 1000) {

      # This function will calculate the probability of no
      # coalescent events in a sample of n chromosomes from
      # a diploid population of size N

      result = 1.0

      for (i in 1:(n-1)) {
         result <- result * (2*N - i)/ (2*N)
         }
   
      result
      }

# Look at effect of sample size
   Pca(n = 2)
   Pca(n = 4)
   Pca(n = 10)
   Pca(n = 20)

   # Look at effect of population size
   Pca(n = 2)
   Pca(n = 2, 10)
   Pca(n = 2, 100)
   Pca(n = 2, 1000)

   # Or compare exact and approximate results
   Pca(n = 5, 10)
   Pca(n = 5, 100)

   PcaExact(n = 5, 10)
   PcaExact(n = 5, 100)

 # The implementations below are conceptually simple, but crude
   # they will be very slow when N or n is large

   P2 <- function (S, N = 1000, mu = 10^-6) {

      # This function will calculate the probability of
      # exactly S segregating sites in a sample of two
      # chromosomes from a population of size N and
      # with mutation rate mu

      theta <- 4 * N * mu
 
      (theta / (1 + theta)) ^ S * (1/(1 + theta))
      } 

   Qn <- function (S, n, N = 1000, mu = 10^-6) {

      # This function will calculate the probability that
      # exactly S mutations occur before the next coalescent
      # event in a sample of n chromosomes from a population
      # of size N with mutation rate mu

      theta <- 4 * N * mu

      (theta / (theta + n - 1)) ^ S * ((n - 1)/(n - 1 + theta))
      }

   Pn <- function(S, n, N = 1000, mu = 10^-6) {
 
      # This function will calculate the probability that
      # exactly S mutations are present in a sample of n
      # chromosomes from a population of size N with
      # mutation rate mu

      if (n > 2) {

         result <- 0

         for (i in 0:S)
            result <- result + Pn(i, n-1, N, mu) * Qn(S-i, n, N, mu)
         }
      else
         result <- P2(S, N, mu)

      result
      }

SingleLocusTree <- function(n = 10, N = 1000, debug = FALSE) {

   # Length of each branch in the tree
   lengths <- rep(0, n)

   # Number of descendants for each node
   descendants <- rep(1, n)

   # Parent node for each node
   parent <- rep(0, n)

   # Descendant nodes for each node
   activeNodes <- 1:n

   nextNode = n + 1

   for (i in n:2)
      {
      # sample time to next coalescent event
      t <- rexp(1, rate = (i*(i-1)/2)/N)

      # all current nodes increase in length by time t
      lengths[activeNodes] <- lengths[activeNodes] + t

      # next we select two nodes to coalesce at random
      coalescence <- sample(activeNodes, 2)

      if (debug)
         cat("After ", t, "generations: Nodes ", coalescence[1], " and ", coalescence[2], "coalesce\n")

      # assign their parental nodes
      parent[coalescence] = nextNode

      # remove these from the list of active nodes
      activeNodes <- setdiff(activeNodes, coalescence)

      # create a new node
      parent[nextNode] <- 0
      lengths[nextNode] <- 0
      descendants[nextNode] <- descendants[coalescence[1]] + descendants[coalescence[2]]

      # And update list of active nodes
      activeNodes <- c(activeNodes, nextNode)

      nextNode <- nextNode + 1
      }

  list(descendants = descendants, parent = parent, lengths = lengths)
  }

GetTreeSpectrum <- function(tree, n)
  {
  spectrum <- rep(0, n)

  for (i in 1:length(tree$descendants))
      spectrum[tree$descendants[i]] <- spectrum[tree$descendants[i]] + tree$lengths[i]

  spectrum
  }

  n <- 10
   spectrum <- rep(0, n)
 
   for (i in 1:100)
      {
      tree <- SingleLocusTree(n)
      spectrum <- GetTreeSpectrum(tree, n) + spectrum
      }

   spectrum <- spectrum / 100

TwoLocusTree <- function(n, N = 1000, r = 0.0001, debug = FALSE)
      {
      # This function generates trees for two linked loci

      # Then we initialize the tree with n sequences for
      # the right hand locus
      rightDescendants <- rep(1, n)
      rightLengths     <- rep(0, n)
      rightParent      <- rep(0, n)

      # First we initialize the tree with n sequences for
      # the left hand locus
      leftDescendants  <- rep(1, n)
      leftLengths      <- rep(0, n)
      leftParent       <- rep(0, n)

      activePairs       <- 1:n
      activeLeft        <- vector()
      activeRight       <- vector()

      nextNode <- n + 1

      startN <- n

      while (length(activePairs) + max(length(activeLeft), length(activeRight)) > 1) {

          # Calculate probability of recombination
          prec <- length(activePairs) * r

          # Calculate probability of coalescence
          n    <- length(activePairs) + length(activeLeft) + length(activeRight)
          pca  <- n * (n - 1) / (2 * N)

          # Calculate time to next event
          t    <- rexp(1, pca + prec)

          if (debug) cat(length(activePairs), length(activeLeft), length(activeRight),
                         "After ", t, " generations: ")

          # Update times for all nodes
          indexRight <- c(activeRight, activePairs)
          rightLengths[indexRight] <- rightLengths[indexRight] + t

          indexLeft <- c(activeLeft, activePairs)
          leftLengths[indexLeft] <- leftLengths[indexLeft] + t

          # Was this a coalescence or recombination event?
          if (rbinom(1, 1, pca / (pca + prec)) == 1) {

             # Pick two nodes at random ...
             coalescent <- sample(c(activeLeft, activeRight, activePairs), 2)

             if (debug) cat("Nodes ", coalescent[1], " and ", coalescent[2], " coalesce into ", nextNode, "\n")

             # Set parents appropriately
             rightParent[intersect(coalescent, indexRight)] <- nextNode
             leftParent[intersect(coalescent, indexLeft)] <- nextNode

             # Create ancestor node
             rightParent[nextNode] <- leftParent[nextNode] <- 0
             rightLengths[nextNode] <- leftLengths[nextNode] <- 0
             rightDescendants[nextNode] <- sum(rightDescendants[coalescent])
             leftDescendants[nextNode]  <- sum(leftDescendants[coalescent])

             # Delete coalescing nodes from appropriate list
             activePairs <- setdiff(activePairs, coalescent)
             activeLeft  <- setdiff(activeLeft, coalescent)
             activeRight <- setdiff(activeRight, coalescent)

             # Add new ancestral node to appropriate list
             if (rightDescendants[nextNode] == 0) {
                activeLeft <- c(activeLeft, nextNode)
             } else if (leftDescendants[nextNode] == 0) {
                activeRight <- c(activeRight, nextNode)
             } else
                activePairs <- c(activePairs, nextNode)

             # Check for special situations where one locus finishes coalescing
             if (length(activePairs) == 0)  {
                if (length(activeRight) == 1) activeRight <- vector()
                if (length(activeLeft) == 1) activeLeft <- vector()
                }

             if ((length(activePairs) == 1))  {
                if (length(activeRight) == 0) {
                   activeLeft <- c(activeLeft, activePairs)
                   activePairs <- vector()
                   }
                if (length(activeLeft) == 0) {
                   activeRight <- c(activeRight, activePairs)
                   activePairs <- vector()
                   }
                }

             nextNode <- nextNode + 1

          } else {

             # Pick one node at random to recombine
             recombinant <- ifelse(length(activePairs)==1, activePairs, sample(activePairs, 1))

             if (debug) cat("Node ", recombinant, " splits into ", nextNode, " and ", nextNode + 1, "\n")

             # Assign different ancestors to each stretch
             rightParent[recombinant] <- nextNode
             leftParent[recombinant]  <- nextNode+1

             activeRight <- c(activeRight, nextNode)
             activeLeft <- c(activeLeft, nextNode + 1)

             # Create ancestor for right portion
             rightParent[nextNode] <- leftParent[nextNode] <- 0
             rightLengths[nextNode] <- leftLengths[nextNode] <- 0
             rightDescendants[nextNode] <- rightDescendants[recombinant]
             leftDescendants[nextNode]  <- 0
             nextNode <- nextNode + 1

             # Create ancestor for left portion
             rightParent[nextNode] <- leftParent[nextNode] <- 0
             rightLengths[nextNode] <- leftLengths[nextNode] <- 0
             rightDescendants[nextNode] <- 0
             leftDescendants[nextNode]  <- leftDescendants[recombinant]
             nextNode <- nextNode + 1

             # Remove recombinant node from active list
             activePairs <- setdiff(activePairs, recombinant)
             }
          }

      left <- list(descendants = leftDescendants, parent = leftParent, lengths = leftLengths)
      right <- list(descendants = rightDescendants, parent = rightParent, lengths = rightLengths)

      list(left = left, right = right)
      }


 
