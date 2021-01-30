library(dplyr)
library(stringr)
library(argparser)
library(optparse)

hmmparser <- function(filename) {
  infile <- read.csv(filename, header = F, fill = T, row.names = 1, stringsAsFactors = F)
  S <- na.omit(unlist(infile["S",]))
  V <- na.omit(unlist(infile["V",]))
  tmatrix <- matrix(na.omit(unlist(as.numeric(infile["A",]))), nrow = 2, ncol = 2)
  rownames(tmatrix) <- S
  p0 <- na.omit(unlist(as.numeric(infile["m",])))
  ematrix <- matrix(na.omit(unlist(as.numeric(infile["B",]))), nrow = 2, ncol = 5)
  colnames(ematrix) <- V
  rownames(ematrix) <- S
  return(list(tmatrix, ematrix, p0))
}

markov <- function(hmm, chainlength) {
  tmatrix <- hmm[[1]]
  p0 <- hmm[[3]]
  S <- rownames(tmatrix)
  chain <- rep(0, chainlength)
  chain[1] <- sample(S, 1, prob=p0)
  for (i in 2:chainlength) {
    chain[i] <- sample(S, 1, prob=tmatrix[chain[i-1],])
  }
  return(chain)
}

hiddenmarkov <- function(hmm, chainlength) {
  tmatrix <- hmm[[1]]
  pemit <- hmm[[2]]
  p0 <- hmm[[3]]
  hiddenchain <- markov(hmm, chainlength)
  emit.states <- colnames(pemit)
  observed <- lapply(hiddenchain, function(x) {
    sample(emit.states, 1, prob = pemit[x, ])
  })
  return(list(hiddenchain, observed))
}

GCwindowmaker <- function(filename, windowsize, outfile) {
  infile <- readLines(filename)
  sequence <- unlist(sapply(2:length(infile), function(x) {unlist(str_split(infile[x], ""))}, simplify = T))
  sequencesize <- length(sequence)
  GC <- c()
  ind <- 0
  ss <- c("G", "C")
  for (i in seq(1, sequencesize, windowsize)) {
    ind <- ind + 1
    gc <- 0
    for (nt in i:(windowsize*ind)) {
      if (is.element(sequence[nt], ss)) {
        gc <- gc + 1
      }
    }
    gc <- gc / windowsize
    GC <- append(GC, gc)
  }
  
  #Examine GC Distribution
  pdf(outfile)
  par(mfrow = c(1, 2))
  hist(GC)
  plot(1:length(GC), GC, type = "l", xlab = "Window", main = "GC content by window in sequence")
  dev.off()
  length(GC)
  return(GC)
}

forwardalgorithm <- function(V, hmm) {
  tmatrix <- hmm[[1]]
  ematrix <- hmm[[2]]
  p0 <- hmm[[3]]
  #Initialisation
  p = 1
  alpha <- matrix(0, nrow = length(V), ncol = ncol(tmatrix))
  alpha[1, ] <- p0*ematrix[,V[1]]
  #Recursion
  for (n in 2:length(V)) {
    obs_prob <- 0
    tmp <- (alpha[n-1, ] %*% tmatrix) #should this be summed?
    alpha[n,] <- ematrix[ ,V[n]] * tmp
  }
  #Termination
  L <- sum(alpha[length(V), ])
}

scaledforwardalgorithm <- function(V, hmm) {
  tmatrix <- hmm[[1]]
  ematrix <- hmm[[2]]
  p0 <- hmm[[3]]
  # Forward algorithm, with scaling for numerical stability.
  states <- 1:nrow(tmatrix) # number of states
  L <- length(V)
  alpha <- matrix(0, ncol=nrow(tmatrix), nrow=L)
  scaling <- numeric(L)
  last_alpha <- alpha[1,] <- p0
  
  # main recursion component
  for (i in 1:length(V)) {
    for (k in states) {
      alpha[i, k] <- sum(last_alpha * tmatrix[,k]) * ematrix[k, V[i]]
    }
    scaling[i] <- c_i <- sum(alpha[i,])
    last_alpha <- alpha[i,] <- alpha[i,]/c_i
  }
  dimnames(alpha) <- list(seq_along(V), rownames(ematrix))
  return(list(alpha=alpha, scaling=1/scaling, likelihood = sum(log(scaling))))
}

backwardalgorithm <- function(V, hmm) {
  tmatrix <- hmm[[1]]
  ematrix <- hmm[[2]]
  p0 <- hmm[[3]]
  #Initialisation
  beta <- matrix(0, nrow = length(V), ncol = ncol(tmatrix))
  beta[length(V), ] <- rep(1, ncol(beta))
  #Recursion
  for (n in seq((length(V)-1), 1, -1)) {
    beta[n,] <- ematrix[ ,V[n+1]] %*% tmatrix * beta[n+1,]
  }
  #Termination
  return(beta)
}

scaledbackwardalgorithm <- function(V, hmm) {
  tmatrix <- hmm[[1]]
  ematrix <- hmm[[2]]
  p0 <- hmm[[3]]
  # Do forward algorithm to get scaling
  scaling <- scaledforwardalgorithm(V, hmm)[[2]]
  # Backward algorithm, with scaling for numerical stability.
  states <- 1:(nrow(tmatrix))
  L <- length(V)
  beta <- matrix(0, ncol=nrow(tmatrix), nrow=L)
  last_beta <- beta[L,] <- rep(scaling[L], length(states))
  
  # main recursion component
  for (i in seq(L-1, 1)) {
    for (k in states) {
      beta[i, k] <- sum(ematrix[, V[i+1]] * last_beta * tmatrix[k, ])
    }
    beta[i,] <- last_beta <- scaling[i]* beta[i,]
  }
  dimnames(beta) <- list(seq_along(V), rownames(ematrix))
  return(beta)
}

randomtransitionmatrix <- function(num_states){
  tprobs <- c()
  for (i in 1:num_states) {
    splitter <- c(0, sort(runif((num_states - 1))), 1)
    for (j in 2:length(splitter)) {
      tprobs <- append(tprobs, (splitter[j] - splitter[j-1]))
    }
  }
  tmatrix <- matrix(tprobs, byrow = TRUE, nrow = num_states, ncol = num_states)
}

randomemissionmatrix <- function(num_states, num_values) {
  eprobs <- c()
  for (i in 1:num_states) {
    splitter <- c(0, sort(runif((num_values - 1))), 1)
    for (j in 2:length(splitter)) {
      eprobs <- append(eprobs, (splitter[j] - splitter[j-1]))
    }
  } 
  ematrix <- matrix(eprobs, byrow = TRUE, nrow = num_states, ncol = num_values)
}

randomp0 <- function(num_states) {
  p0 <- c()
  splitter <- c(0, sort(runif((num_states - 1))), 1)
  for (j in 2:length(splitter)) {
    p0 <- append(p0, (splitter[j] - splitter[j-1]))
  }
  return(p0)
}

normalise <- function(matrix) {
  for (i in 1:nrow(matrix)){
    sum <- sum(matrix[i,])
    for (j in 1:ncol(matrix)) {
      matrix[i,j] <- matrix[i,j] / sum
    }
  }
  return(matrix)
}

baumwelchalgorithm <- function(V, num_states) {
  num_values <- length(unique(V))
  a <- randomtransitionmatrix(num_states)
  b <- randomemissionmatrix(num_states, num_values)
  initial_distribution <- randomp0(num_states)
  count <- 0
  converged <- F
  hmm <- list(a, b, initial_distribution)
  while (converged == F) {
    N <- length(V)
    M <- nrow(a)
    K <- ncol(b)
    alpha <- scaledforwardalgorithm(V, hmm)[[1]]
    beta <- scaledbackwardalgorithm(V, hmm)
    xi <- array(0, dim=c(M, M, N-1))
    
    likelihood <- scaledforwardalgorithm(V, hmm)[[3]]
    print(likelihood)
    for(n in 1:N-1){
      denominator = ((alpha[n,] %*% a) * b[,V[n+1]]) %*% matrix(beta[n+1,]) 
      for(s in 1:M){
        numerator = alpha[n,s] * a[s,] * b[,V[n+1]] * beta[n+1,]
        xi[s,,n]=numerator/as.vector(denominator)
      }
    }
    
    xi.all.t = rowSums(xi, dims = 2)
    a = xi.all.t/rowSums(xi.all.t)
    
    newp0 <- c(0,0)
    for (j in 1:nrow(a)) {
      newp0 <- newp0 + (alpha[1, ] * a[j,] * b[,V[1]] * beta[2,])
    }
    newp0 <- newp0 / sum(newp0)
    initial_distribution <- newp0
    
    gamma = apply(xi, c(1, 3), sum)  
    gamma = cbind(gamma, colSums(xi[, , N-1]))
    for(l in 1:K){
      b[, l] = rowSums(gamma[, which(V==l)])
    }
    b = b/rowSums(b)
    
    newhmm <- list(a, b, initial_distribution)
    newlikelihood <- scaledforwardalgorithm(V, newhmm)[[3]]
    delta <- -1*(likelihood - newlikelihood)
    hmm <- newhmm
    if (delta < (0.0005)) {
      count <- count + 1
      if (count > 100) {
        converged <- T
      }
    }
    else {
      count <- 0
    }
    
  }
  return(list(a, b, initial_distribution))
}

viterbialgorithm <- function(V,hmm) {
  tmatrix <- hmm[[1]]
  ematrix <- hmm[[2]]
  p0 <- hmm[[3]]
  
  N = length(V)
  M = nrow(tmatrix)
  prev = matrix(0, N-1, M)
  omega = matrix(0, M, N)
  
  omega[, 1] = log(p0 * ematrix[, V[1]])
  for(n in 2:N){
    for(s in 1:M) {
      probs = omega[, n-1] + log(tmatrix[, s]) + log(ematrix[s, V[n]])
      prev[n - 1, s] = which.max(probs)
      omega[s, n] = max(probs)
    }
  }
  
  S = rep(0, N)
  last_state=which.max(omega[,ncol(omega)])
  S[1]=last_state
  
  j=2
  for(i in (N-1):1){
    S[j]=prev[i,last_state] 
    last_state=prev[i,last_state] 
    j=j+1
  }
  
  
  S=rev(S)
  S <- S-1
  return(S)
  
}

#Create parameter file
createparams <- function(states,values,p0, tmatrix, ematrix, outputfile) {
  params <- rbind(S = states, V = values, m = p0, A = tmatrix, B = ematrix)
  write.csv(params, outputfile, row.names = FALSE)
}

#Create emission sequence
emissionseq <- function()
  
  ###################Operations####################
par(mfrow = c(1,1), mar=c(4,4,2,2))
option_list = list(
  make_option(c("-a", "--action"), type="character", default="everything", 
              help="desired action from script", metavar="character"),
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="HMM parameter file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-l", "--chainlength"), type="integer", default="out.txt", 
              help="Length of desired sequence from HMM", metavar="integer"),
  make_option(c("-p", "--plotfile"), type="character", default="out.pdf", 
              help="Name of output file for any plots", metavar="character"),
  make_option(c("-v", "--sequence"), type="character", default="out.pdf", 
              help="File containing mitted sequence for inspection with forward, baumwelch etc algorithms", metavar="character"),
  make_option(c("-s", "--states"), type="integer", default="out.pdf", 
              help="number of states to assume for Baum-Welch algorithm", metavar="integer")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (opt$action == "emissionseq") {
  hmm <- hmmparser(opt$file)
  length <- opt$chainlength
  hiddenmarkovmodel <- hiddenmarkov(hmm, length)
  emittedseq <- hiddenmarkovmodel[[2]]
  foroutput <- (unlist(lapply(emittedseq, paste, collapse="\n")))
  sink(opt$out)
  invisible(writeLines(foroutput))
  sink()
}

if (opt$action == "emissionseqandplots") {
  hmm <- hmmparser(opt$file)
  length <- opt$chainlength
  hiddenmarkovmodel <- hiddenmarkov(hmm, length)
  emittedseq <- unlist(hiddenmarkovmodel[[2]])
  hiddenseq <- unlist(hiddenmarkovmodel[[1]])
  foroutput <- (unlist(lapply(emittedseq, paste, collapse="\n")))
  sink(opt$out)
  invisible(writeLines(foroutput))
  sink()
  combined <- as.numeric(c((hiddenseq), (emittedseq)))
  print(combined)
  pdf(opt$plotfile) 
  plot(1:length, emittedseq, type = "l", ylim = c(min(combined),max(combined)), col="red", ylab = "V (red) and S (blue)", xlab = "Position in Chain")
  lines(1:length, hiddenseq, col = "blue")
  dev.off() 
}

if (opt$action == "forwardalgorithm") {
  emissionseq <- readLines(opt$sequence)
  hmm <- hmmparser(opt$file)
  likelihood <- forwardalgorithm(emissionseq, hmm)
  print(likelihood)
}

if (opt$action == "scaledforwardalgorithm") {
  emissionseq <- readLines(opt$sequence)
  hmm <- hmmparser(opt$file)
  likelihood <- scaledforwardalgorithm(emissionseq, hmm)[[3]]
  print(likelihood)
}

if (opt$action == "baumwelchalgorithm") {
  emissionseq <- as.numeric(readLines(opt$sequence))
  print("starting baum welch")
  baumwelch <- baumwelchalgorithm(emissionseq, opt$states)
  print("finished baum welch")
  states <- nrow(baumwelch[[1]])
  values <- ncol(baumwelch[[2]])
  A <- c()
  for (i in 1:states) {
    A <- append(A, baumwelch[[1]][,i])
  }
  B <- c()
  for (i in 1:values) {
    B <- append(B, baumwelch[[2]][,i])
  }
  m <- baumwelch[[3]]
  V <- sort(unique(emissionseq))
  S <- 1:opt$states
  n <- max(length(A), length(B))
  length(A) <- n
  length(B) <- n
  length(m) <- n
  length(V) <- n
  length(S) <- n
  newparams <- rbind(A = A, B = B, m = m, V = V, S = S)
  write.csv(newparams, opt$out)
}

if (opt$action == "viterbialgorithm") {
  emissionseq <- as.numeric(readLines(opt$sequence))
  hmm <- hmmparser(opt$file)
  print("starting viterbi")
  viterbi <-viterbialgorithm(emissionseq, hmm)
  print("finished viterbi")
  sink(opt$out)
  writeLines(as.character(viterbi))
  sink()
}

if (opt$action == "GCwindows") {
  windowsize <- as.numeric(opt$states)
  infile <- opt$file
  pdffile <- opt$plotfile
  GC <- GCwindowmaker(infile, windowsize, pdffile)
  sink(opt$out)
  writeLines(as.character(GC))
  sink()
}

if (opt$action == "sequencebinner") {
  GC <- as.numeric(readLines(opt$file))
  Bin1 <- 0.24
  Bin2 <- 0.33
  Bin3 <- 0.41
  Bin4 <- 0.48
  Bin5 <- 1
  
  binnedGC <- c()
  for (window in GC) {
    if (window < Bin1) {
      binnedGC <- append(binnedGC, 1)
    }
    else if (window < Bin2) {
      binnedGC <- append(binnedGC, 2)
    }
    else if (window < Bin3) {
      binnedGC <- append(binnedGC, 3)
    }
    else if (window < Bin4) {
      binnedGC <- append(binnedGC, 4)
    }
    else {
      binnedGC <- append(binnedGC, 5)
    }
  }
  sink(opt$out)
  writeLines(as.character(binnedGC))
  sink()
}



