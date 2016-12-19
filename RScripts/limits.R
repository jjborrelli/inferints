# Functions for the LIMITS (Learning Interactions through MIcrobial Time Series) Algorithm
library(magrittr)

spp <- 6
fitmod <- r.m3l2A %>% data_setup(xi = spp) %>% data_split %>% step_forward(x = spp, errorthres = 1)
fitmod
data <-  r.m3l2A %>% data_setup(xi = spp) %>% data_split

data_setup <- function(data, xi){
  # take in matrix with columns = taxa, rows = timepoints
  # make sure we are dealing with relative abundances
  #if(sum(round(rowSums(data), 2) != 1.00) > 0){stop("Needs relative abundance")}
  
  rows <- as.numeric(sapply(strsplit(row.names(data), "X"), "[[", 2))
  chrows <- rows[which((rows-1 )%in% rows)]
  colnames(data) <- sapply(strsplit(colnames(data), ";"), tail, 1)
  
  resp <- data[row.names(data) %in% paste0("X", chrows), xi]
  resp <- resp[resp!=0]
  
  M <- apply(data[row.names(data) %in% paste0("X", chrows-1),], 2, function(x){x - median(x)})[1:nrow(data[row.names(data) %in% paste0("X", chrows-1),]),]
  vi <- log(resp[2:length(resp)]) - log(resp[1:(length(resp)-1)])
  
  return(data.frame(vi = vi, M))
}

data_split <- function(data){
  len1 <- nrow(data)
  len2 <- ceiling(len1/2)
  times1 <- sort(sample(2:len1, len2))
  times2 <- times1 -1
  
  return(list(train = data[times1,], test = data[times2,]))
}



step_forward <- function(data, x, errorthres = 1){
  xnam <- colnames(data$train)[-1]
  
  f1 <- as.formula(paste("vi ~ ", paste(xnam[x], collapse= "+")))
  fit.init <- lm(f1, data = data$train)
  #err <- mean((data$test$vi-predict(fit.init, newdata = data$test))^2)/var(data$train$vi)
  err <- sum((data$test$vi-predict(fit.init, newdata = data$test))^2)
  
  active <- x
  cond <- FALSE
  
  while(!cond){
    form <- sapply((1:length(xnam))[-active], function(y){
      as.formula(paste("vi ~ ", paste0(c(xnam[active], xnam[y]), collapse = "+")))
    })
    
    fit.step <- lapply(form, function(y){
      fit.step <- lm(y, data = data$train)
    })
    
    err.step <- sapply(fit.step, function(y){
      #mean((data$test$vi-predict(y, newdata = data$test))^2)/var(data$test$vi)
      sum((data$test$vi-predict(y, newdata = data$test))^2)
    })
    
    min.err <- which.min(err.step)
    err.diff <- 100 * (err - err.step[min.err])/err
    
    if(err.diff > errorthres){
      err <- err.step[min.err]
      active <- c(active,(1:length(xnam))[-x][min.err])
    }else{
      cond <- TRUE
    }
    
  }
  
  form.fin <- as.formula(paste("vi ~ ", paste0(c(xnam[active]), collapse = "+")))
  fit.fin <- lm(form.fin, data = data$train)
  return(fit.fin)
}

step_forward.alt <- function(data, x, errorthres = 1){
  trainMAT <- as.matrix(data[[1]][,-1])
  trainRESP <- data[[1]]$vi
  testMAT <- as.matrix(data[[2]][,-1])
  testRESP <- data[[2]]$vi
  
  #
  allsp <- 1:ncol(trainMAT)
  active <- x
  
  # Fit initial regression
  ci <- ginv(trainMAT[,c(x)]) %*% trainRESP
  e1 <- mean((testMAT[,c(x)] %*% ci-testRESP)^2)/var(testRESP)
  
  cond <- FALSE
  while(!cond){
    errs <- sapply(allsp[-active], function(y){
      ci2 <- ginv(trainMAT[,c(active,y)]) %*% trainRESP
      e2 <- mean(((testMAT[,c(active,y)] %*% ci2)-testRESP)^2)/var(testRESP)
      return(e2)
    })
    
    ediff <- 100 * (e1 - min((errs)))/e1
    
    if(ediff >= errorthres){
      e1 <- min(errs)
      active <- c(active, allsp[-active][which.min(errs)])
    }else{
      cond <- TRUE
    }
    
    if(length(active) == length(allsp)){
      cond <- TRUE
    }
    
  }
  
  cif <- ginv(trainMAT[,active]) %*% trainRESP
  
  res <- matrix(0, nrow = 1, ncol = ncol(trainMAT))
  res[,active] <- cif
  
  return(res)
}


iter <- 200
spp <- 1
comat <- matrix(0, nrow = iter, ncol = ncol(r.m3l2A))
fitmod <- list()
for(i in 1:iter){
  data <- r.m3l2A %>% data_setup(xi = spp) %>% data_split
  comat[i,] <- step_forward.alt(data = data, x = spp, errorthres = 3)
  print(i)
}
colMeans(comat)

c.vec <- function(idata, spp, iter = 200, e.thres = 1){
  comat <- matrix(0, nrow = iter, ncol = ncol(idata))
  for(i in 1:iter){
    dat <- idata %>% data_setup(xi = spp) %>% data_split
    fitmod <- step_forward(dat, x = spp, errorthres = e.thres)
    co1 <- colnames(dat[[1]])[-1]%in%names(fitmod$coefficients[-1])
    comat[i,co1] <- fitmod$coefficients[-1]
  }
  return(apply(comat, 2, median))
}

c.vec.alt <- function(idata, spp, iter = 200, e.thres = 1){
  comat <- matrix(0, nrow = iter, ncol = ncol(idata))
  for(i in 1:iter){
    data <- idata %>% data_setup(xi = spp) %>% data_split
    comat[i,] <- step_forward.alt(data = data, x = spp, errorthres = e.thres)
  }
  return(apply(comat, 2, median))
}


imat <- matrix(0, ncol(r.m3l2A), ncol(r.m3l2A))
imat.alt <- matrix(0, ncol(r.m3l2A), ncol(r.m3l2A))
for(i in 1:ncol(r.m3l2A)){
  imat[i,] <- c.vec(r.m3l2A, spp = i, iter = 200, e.thres = 1)
  imat.alt[i,] <- c.vec.alt(r.m3l2A, spp = i, iter = 1, e.thres = 1)
}
imat 
imat.alt



m3l6 <- as.matrix(read.csv("Data/M3L6gut.csv", row.names = 1))

r.m3l6 <- m3l6[which(apply(m3l6, 1, function(x) sum(x !=0)) >= 332),]
r.m3l6A <- t(apply(r.m3l6, 2, function(x) x/sum(x)))


imat <- matrix(0, ncol(r.m3l6A), ncol(r.m3l6A))
for(i in 1:ncol(r.m3l6A)){
  imat[i,] <- c.vec(r.m3l6A, spp = i, iter = 200, e.thres = 3)
}
imat 

taxa1 <- c("Akkermansia muciniphilia", "AAlistipes putredinis", "Bacteroides acidifaciens", "Bacteroides fragilis", "Bacteroides stercoris", "Bacteroides thetaiotaomicron", "Bacteroides uniformis", "Bacteroides vulgatus", "Escherichia coli", "Eubacterium rectale", "Faecalibacterium prausnitzii", "Parebacteroides distasonis", "Roseburia intestinalis", "Sisymbrium irio")


library(MASS)

ci <- ginv(as.matrix(fitmod[[1]][,c(2,3,5)]))%*%fitmod[[1]][,1]
e1 <- sum((as.matrix(fitmod[[1]][,c(2,3,5)]) %*% (ci)-fitmod[[1]]$vi)^2)

ci2 <- ginv(as.matrix(fitmod[[1]][,c(2,3,5,6)]))%*%fitmod[[1]][,1]
e2 <- sum((rowSums(as.matrix(fitmod[[1]][,c(2,3,5,6)]) %*% (ci2))-fitmod[[1]]$vi)^2)