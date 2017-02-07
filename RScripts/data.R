
# FUNCTION; Simulate Data

lvfm <- function(times, state, parms){
  cij <- parms$cij
  x.med <- parms$xmed
  ni <- parms$ni
  
  xi <- matrix(state, nrow = nrow(cij), ncol = length(times))
  for(i in 2:length(times)){
    xi[,i] <- exp(ni[i,] + log(xi[,i-1]) + cij %*% (xi[,i-1] - x.med))
    #xi[xi < 1e-5] <- 0
  }
  
  return(xi)
}

getpars <- function(N, L, times){
  testmat <- matrix(0, N, N)
  xbar <- rlnorm(N, 0, .1)
  diag(testmat) <- runif(N, -1.9, -.1)/xbar
  
  for(i in 1:L){
    r1 <- sample(1:N, 1)
    r2 <- sample(1:N, 1)
    testmat[r1, r2] <- rnorm(1,0,.5)
    testmat[r2, r1] <- rnorm(1,0,.5)
  }
  
  allpar <- list(cij = testmat, xmed = xbar, ni = matrix(rnorm(times*N, 0, .1), ncol = N, nrow = times))
  
  return(list(p = allpar, m = xbar))
}

gp1 <- getpars(15, 15, 2000)

matplot(t(lvfm(1:100, gp1$m, gp1$p)), typ = "l")
dyn <- (lvfm(1:1000, gp1$m, gp1$p))
fmdat <- t(apply(dyn, 2, function(x) x/sum(x)))
matplot(fmdat, typ = "l")
xbar * exp(testmat %*% (xbar - xbar))

# FUNCTION: Get taxonomy

taxa_names <- function(x){
  # x is the raw data
  
  splitnames <- strsplit(rownames(x), ";")
  n1 <- data.frame(kingdom = NA, phylum = NA, class = NA, order = NA, family = NA, genus = NA)
  #n1 <- matrix(nrow = length(splitnames), ncol = length(splitnames[[1]]))
  for(i in 1:length(splitnames)){
    #n1[,i] <- sapply(splitnames, "[[", i)
    n1[i,1:length(splitnames[[i]])] <- splitnames[[i]]
  }
  
  
  if(!is.na(n1$kingdom[1])){n1$kingdom <- sapply(strsplit(n1$kingdom, "k__"), function(x){if(length(x) == 2){x[2]}else{NA}})}
  if(!is.na(n1$phylum[1])){n1$phylum <- sapply(strsplit(n1$phylum, "p__"), function(x){if(length(x) == 2){x[2]}else{NA}})}
  if(!is.na(n1$class[1])){n1$class <- sapply(strsplit(n1$class, "c__"), function(x){if(length(x) == 2){x[2]}else{NA}})}
  if(!is.na(n1$order[1])){n1$order <- sapply(strsplit(n1$order, "o__"), function(x){if(length(x) == 2){x[2]}else{NA}})}
  if(!is.na(n1$family[1])){n1$family <- sapply(strsplit(n1$family, "f__"), function(x){if(length(x) == 2){x[2]}else{NA}})}
  if(!is.na(n1$genus[1])){n1$genus <- sapply(strsplit(n1$genus, "g__"), function(x){if(length(x) == 2){x[2]}else{NA}})}
  
  n1 <- t(apply(n1, 1, function(x){x[is.na(x)] <- paste("na_", tail(unlist(x[!is.na(x)]), 1), sep = "__");return(x)}))
  
  return(n1)
}


# FUNCTION: Get times

timepoints <- function(x){
  # x is the raw data
  
  tp <- as.numeric(sapply(strsplit(colnames(x), "X"), "[[", 2))
  before <- tp[2:length(tp)] - 1
  
  t2 <- tp[-1][before %in% tp]
  t2ind <- which(tp %in% t2)
  t1 <- t2-1
  t1ind <- t2ind-1
  tpairs <- cbind(t1, t2, t1ind, t2ind)
  return(tpairs)
}

# FUNCTION: Rearrange the data

data_process <- function(x, level, o = FALSE){
  # x is the raw data
  # level is taxonomic level for column names (e.g. "phylum")
  x <- apply(x, 2, function(y) y/sum(y))
  rownames(x) <- taxa_names(x)[,level]
  x <- t(x)
  if(o){
    other <- rowSums(x[,-order(apply(x, 2, median), decreasing = T)[1:min(10, ncol(x))]])
    return(data.frame(x[,order(apply(x, 2, median), decreasing = T)[1:10]], other))
  }else{
    return(data.frame(x))
  }
}

# FUNCTION: Setup data for regression

data_setup <- function(data, level, xi, o = FALSE){
  # take in matrix with columns = taxa, rows = timepoints
  # make sure we are dealing with relative abundances
  tp <- timepoints(data)
  data <- data_process(data, level, o)
  
  
  resp <- data[, xi]
  vi <- log(resp[tp[,"t2ind"]]) - log(resp[tp[,"t1ind"]])

  M <- apply(data[tp[,"t1ind"],], 2, function(x){x - median(x)})[!is.nan(vi) & !is.infinite(vi),]
  
  vi <- vi[!is.nan(vi) & !is.infinite(vi)]
  
  return(data.frame(vi = vi, M))
}

data_split <- function(data){
  len1 <- nrow(data)
  len2 <- ceiling(len1/2)
  brk1 <- sample(1:len1, len2)
  vec1 <- rep(FALSE, len1)
  vec1[brk1] <- TRUE
  
  return(list(train = data[vec1,], test = data[!vec1,]))
}

data_split2 <- function(data){
  len1 <- nrow(data)-1
  len2 <- ceiling(len1/2)
  brk1 <- sample(1:len1, len2)
  vec1 <- rep(FALSE, len1)
  vec1[brk1] <- TRUE
  
  return(which(vec1))
}


# FUNCTION LIMITS
limits <- function(fmdat, errT, bag, sim = F){
  if(!sim){
    tp <- as.numeric(rownames(fmdat))
    before <- tp[2:length(tp)] - 1
    t1 <- which(before %in% tp)
    t2 <- which(before %in% tp)+1
  }
  imat <- matrix(nrow = ncol(fmdat), ncol = ncol(fmdat))
  for(xi in 1:ncol(fmdat)){
    resp <- fmdat[, xi]
    
    if(sim){
      vi <- log(resp[2:length(resp)]) - log(resp[1:(length(resp)-1)])
      M <-  apply(fmdat, 2, function(x){x[1:(length(resp) -1)] - median(x)})[!is.nan(vi) & !is.infinite(vi),]
    }else{
      vi <- log(resp[t2]) - log(resp[t1])
      M <- apply(fmdat, 2, function(x){x[t1] - median(x)})[!is.nan(vi) & !is.infinite(vi),]
    }
    
    vi <- vi[!is.nan(vi) & !is.infinite(vi)]
    
    rmat <- matrix(nrow = bag, ncol = ncol(fmdat))
    for(i in 1:bag){
      ds1 <- data_split(data.frame(vi, M))
      rmat[i,] <- step_forward.alt(ds1, xi, errT)
    }
    imat[xi,] <- apply(rmat, 2, median)
  }

  return(imat * mean(apply(fmdat, 2, median)))
}

testdat <- read.csv("~/Desktop/testdat.csv", header = F)
limits(testdat, 5, 100, sim = T)


# Caporaso Male Data Set

m3l2 <- read.csv("Data/m3l2.csv", row.names = 1)
m3l3 <- read.csv("Data/m3l3.csv", row.names = 1)
m3l4 <- read.csv("Data/m3l4.csv", row.names = 1)
m3l5 <- read.csv("Data/m3l5.csv", row.names = 1)
m3l6 <- read.csv("Data/m3l6.csv", row.names = 1)


# Caporaso Female Data Set

f4l2 <- read.csv("Data/f4l2.csv", row.names = 1)
f4l3 <- read.csv("Data/f4l3.csv", row.names = 1)
f4l4 <- read.csv("Data/f4l4.csv", row.names = 1)
f4l5 <- read.csv("Data/f4l5.csv", row.names = 1)
f4l6 <- read.csv("Data/f4l6.csv", row.names = 1)

