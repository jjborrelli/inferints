
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
  
  n1 <- t(apply(n1, 1, function(x){x[is.na(x)] <- paste("f", tail(unlist(x[!is.na(x)]), 1), sep = "__");return(x)}))
  
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
  times1 <- sort(sample(2:len1, len2))
  times2 <- times1 -1
  
  return(list(train = data[times1,], test = data[times2,]))
}

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

