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
  ci <- ginv(trainMAT[,c(x,2)]) %*% trainRESP
  e1 <- mean(((testMAT[,c(x,2)] %*% ci)-testRESP)^2)/var(testRESP)
  
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

fmdat <- read.csv("Data/m3unionFM.csv", row.names = 1)
head(fmdat)

tp <- as.numeric(rownames(fmdat))
before <- tp[2:length(tp)] - 1
t2 <- tp[-1][before %in% tp]
t2ind <- which(tp %in% t2)
t1 <- t2-1
t1ind <- t2ind-1
tpairs <- cbind(t1, t2, t1ind, t2ind)

strt <- Sys.time()
errT.1 <- .1
errT1 <- 1
errT2 <- 2
errT3 <- 3
errT5 <- 5
imat <- matrix(nrow = ncol(fmdat), ncol = ncol(fmdat))
imat1 <- matrix(nrow = ncol(fmdat), ncol = ncol(fmdat))
imat2 <- matrix(nrow = ncol(fmdat), ncol = ncol(fmdat))
imat3 <- matrix(nrow = ncol(fmdat), ncol = ncol(fmdat))
imat5 <- matrix(nrow = ncol(fmdat), ncol = ncol(fmdat))
for(xi in 1:ncol(fmdat)){
  resp <- fmdat[, xi]
  #vi <- log(resp[tpairs[,"t2ind"]]) - log(resp[tpairs[,"t1ind"]])
  vi <- log(resp[2:length(resp)]) - log(resp[1:(length(resp)-1)])
  #M <- apply(fmdat, 2, function(x){x[tpairs[,"t1ind"]] - median(x)})[!is.nan(vi) & !is.infinite(vi),]
  M <-  apply(fmdat, 2, function(x){x[1:(length(resp) -1)] - median(x)})[!is.nan(vi) & !is.infinite(vi),]
  
  vi <- vi[!is.nan(vi) & !is.infinite(vi)]
  
  rmat <- matrix(nrow = 200, ncol = ncol(fmdat))
  rmat1 <- matrix(nrow = 200, ncol = ncol(fmdat))
  rmat2 <- matrix(nrow = 200, ncol = ncol(fmdat))
  rmat3 <- matrix(nrow = 200, ncol = ncol(fmdat))
  rmat5 <- matrix(nrow = 200, ncol = ncol(fmdat))
  for(i in 1:200){
    ds1 <- data_split(data.frame(vi, M))
    rmat[i,] <- step_forward.alt(ds1, xi, errT.1)
    rmat1[i,] <- step_forward.alt(ds1, xi, errT1)
    rmat2[i,] <- step_forward.alt(ds1, xi, errT2)
    rmat3[i,] <- step_forward.alt(ds1, xi, errT3)
    rmat5[i,] <- step_forward.alt(ds1, xi, errT5)
  }
  imat[xi,] <- apply(rmat, 2, median)
  imat1[xi,] <- apply(rmat1, 2, median)
  imat2[xi,] <- apply(rmat2, 2, median)
  imat3[xi,] <- apply(rmat3, 2, median)
  imat5[xi,] <- apply(rmat5, 2, median)
}
ends <- Sys.time()
ends - strt
cor.test(imat * apply(fmdat, 2, median), gp1$p$cij)
cor.test(imat1 * apply(fmdat, 2, median), gp1$p$cij)
cor.test(imat2 * apply(fmdat, 2, median), gp1$p$cij)
cor.test(imat3 * apply(fmdat, 2, median), gp1$p$cij)
cor.test(imat5 * apply(fmdat, 2, median), gp1$p$cij)
