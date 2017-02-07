# Functions for the LIMITS (Learning Interactions through MIcrobial Time Series) Algorithm
library(magrittr)

spp <- 6
fitmod <- m3l2 %>% data_setup(xi = spp) %>% data_split %>% step_forward(x = spp, errorthres = 1)
fitmod
data <-  r.m3l2A %>% data_setup(xi = spp) %>% data_split




step_forward <- function(data, x, errorthres = 1){
  xnam <- colnames(data$train)[-1]
  
  f1 <- as.formula(paste("vi ~ ", paste(xnam[x], collapse= "+")))
  fit.init <- lm(f1, data = data$train)
  err <- mean((data$test$vi-predict(fit.init, newdata = data$test))^2)/var(data$train$vi)
  #err <- sum((data$test$vi-predict(fit.init, newdata = data$test))^2)
  
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
      mean((data$test$vi-predict(y, newdata = data$test))^2)/var(data$test$vi)
      #sum((data$test$vi-predict(y, newdata = data$test))^2)
    })
    
    min.err <- which.min(err.step)
    err.diff <- 100 * (err - err.step[min.err])/err
    
    if(err.diff >= errorthres){
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
  trainMAT <- (as.matrix(data[[1]][,-1]))
  trainRESP <- data[[1]]$vi
  testMAT <- (as.matrix(data[[2]][,-1]))
  testRESP <- data[[2]]$vi
  
  #
  allsp <- 1:ncol(trainMAT)
  active <- x
  
  # Fit initial regression
  ci <- ginv(trainMAT[,active]) %*% trainRESP
  e1 <- mean(((ci %*% testMAT[,active] )-testRESP)^2)/var(testRESP)
  
  cond <- FALSE
  while(!cond){
    errs <- sapply(allsp[-active], function(y){
      ci2 <- ginv(trainMAT[,c(active,y)]) %*% trainRESP
      e2 <- mean(((testMAT[,c(active,y)] %*% ci2)-testRESP)^2)/var(testRESP)
      return(e2)
    })
    
    ediff <- 100 * ((e1 - min(errs))/e1)
    
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
  
  cif <- ginv(rbind(testMAT, trainMAT)[,active]) %*% c(testRESP, trainRESP)
  
  res <- matrix(0, nrow = 1, ncol = ncol(trainMAT))
  res[,active] <- cif
  
  return(res)
}



taxa1 <- c("Akkermansia muciniphilia", "AAlistipes putredinis", "Bacteroides acidifaciens", "Bacteroides fragilis", "Bacteroides stercoris", "Bacteroides thetaiotaomicron", "Bacteroides uniformis", "Bacteroides vulgatus", "Escherichia coli", "Eubacterium rectale", "Faecalibacterium prausnitzii", "Parebacteroides distasonis", "Roseburia intestinalis", "Sisymbrium irio")


library(MASS)

ci <- ginv(as.matrix(fitmod[[1]][,c(2,3,5)]))%*%fitmod[[1]][,1]
e1 <- sum((as.matrix(fitmod[[1]][,c(2,3,5)]) %*% (ci)-fitmod[[1]]$vi)^2)

ci2 <- ginv(as.matrix(fitmod[[1]][,c(2,3,5,6)]))%*%fitmod[[1]][,1]
e2 <- sum((rowSums(as.matrix(fitmod[[1]][,c(2,3,5,6)]) %*% (ci2))-fitmod[[1]]$vi)^2)








######################################
######################################
######################################
### Using data from Fisher and Mehta
library(MASS)
fmdat <- read.csv("Data/m3unionFM.csv", row.names = 1)
head(fmdat)

zeros <- apply(fmdat, 1, function(x) sum(x == 0) > 0)
fmdat2 <- fmdat[!zeros,]

tp <- as.numeric(rownames(fmdat2))
before <- tp[2:length(tp)] - 1
t1 <- which(before %in% tp)
t2 <- which(before %in% tp)+1


# is this using simulated data? 
sim <- TRUE

strt <- Sys.time()
errT <- 3
imat <- matrix(nrow = ncol(fmdat), ncol = ncol(fmdat))
for(xi in 1:ncol(fmdat)){
  resp <- fmdat[, xi]
  
  if(sim){
    vi <- log(resp[2:length(resp)]) - log(resp[1:(length(resp)-1)])
    M <-  apply(fmdat, 2, function(x){x[1:(length(resp) -1)] - median(x)})[!is.nan(vi) & !is.infinite(vi),]
  }else{
    vi <- log(resp[t2]) - log(resp[t1])
    M <- apply(fmdat2, 2, function(x){x[t1] - median(x)})[!is.nan(vi) & !is.infinite(vi),]
  }
  
  vi <- vi[!is.nan(vi) & !is.infinite(vi)]
  
  rmat <- matrix(nrow = 100, ncol = ncol(fmdat))
  for(i in 1:100){
    ds1 <- data_split(data.frame(vi, M))
    rmat[i,] <- step_forward.alt(ds1, xi, errT)
  }
  imat[xi,] <- apply(rmat, 2, median)
}
ends <- Sys.time()
ends - strt
imat * apply(fmdat, 2, median)

cor.test(imat * apply(fmdat, 2, median), gp1$p$cij)




######################################
######################################
######################################
### Using data from Fisher and Mehta
library(MASS)
fmdat <- read.csv("Data/m3unionFM.csv", row.names = 1)
head(fmdat)

zeros <- apply(fmdat, 1, function(x) sum(x == 0) > 0)
fmdat2 <- fmdat[!zeros,]

tp <- as.numeric(rownames(fmdat2))
before <- tp[2:length(tp)] - 1
t1 <- which(before %in% tp)
t2 <- which(before %in% tp)+1


sdata <- function(data, ds, xi){
  tp <- as.numeric(rownames(data))
  before <- tp[2:length(tp)] - 1
  t1 <- which(before %in% tp)
  t2 <- which(before %in% tp)+1

  resp <- data[, xi]
  vi <- log(resp[t2]) - log(resp[t1])
  M <- apply(data, 2, function(x){x[t1] - median(x)})[!is.nan(vi) & !is.infinite(vi),]
  vi <- vi[!is.nan(vi) & !is.infinite(vi)]
  
  df.f <- data.frame(vi, M)
  df1 <- df.f[ds, ]
  df2 <- df.f[!ds, ]
  
  return(list(df1, df2))
}


data_split2 <- function(data, t1){
  len1 <- nrow(data[t1,])
  len2 <- ceiling(len1/2)
  brk1 <- sample(1:len1, len2)
  vec1 <- rep(FALSE, len1)
  vec1[brk1] <- TRUE
  
  return((vec1))
}


t1 = 1:99;t2 = 2:100
errT <- 5
rmat <- lapply(1:ncol(fmdat), function(x) matrix(nrow = 200, ncol = ncol(fmdat)))
for(i in 1:200){
  ds1 <- data_split2(fmdat, t1)
  for(xi in 1:ncol(fmdat)){
    ds2 <- sdata(fmdat, ds1, t1, t2, xi)
    rmat[[xi]][i,] <- step_forward.alt(ds2, xi, errT)
  }
}

imat <- sapply(rmat, function(x) apply(x, 2, median))
imat * apply(fmdat, 2, median)
sum(sign(imat * apply(fmdat, 2, median)) != sign(gp1$p$cij))

test <- cbind(as.vector(imat * apply(fmdat, 2, median)), as.vector(gp1$p$cij))
cor.test(test[rowSums(test) != 0,1],test[rowSums(test) != 0,2])
cor.test(as.vector(imat * apply(fmdat, 2, median)), as.vector(gp1$p$cij))
