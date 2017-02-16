
library(MASS)
library(magrittr)
library(reshape2)

lv_ricker <- function(times, state, parms){
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

get_pars <- function(N, L, times){
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

data_setup <- function(data, xi){
  if(is.null(rownames(data))){tp <- 1:nrow(data)}else{tp <- as.numeric(rownames(data))}
  before <- tp[2:length(tp)] - 1
  t1 <- which(before %in% tp)
  t2 <- which(before %in% tp)+1
  #t1 <- 1:(nrow(data)-1)
  #t2 <- 2:nrow(data)
  
  resp <- data[, xi]
  vi <- log(resp[t2]) - log(resp[t1])
  M <- apply(data, 2, function(x){x[t1] - median(x)})[!is.nan(vi) & !is.infinite(vi),]
  vi <- vi[!is.nan(vi) & !is.infinite(vi)]
  
  
  df.f <- data.frame(vi, M)
  
  return(df.f)
}

data_split <- function(data){
  len1 <- nrow(data)
  len2 <- ceiling(len1/2)
  brk1 <- sample(1:len1, len2)
  vec1 <- rep(FALSE, len1)
  vec1[brk1] <- TRUE
  
  df1 <- data[vec1,]
  df2 <- data[!vec1,]
  
  return(list(df1, df2))
}

step_forward <- function(data, x, errorthres = 1){
  trainMAT <- (as.matrix(data[[1]][,-1]))
  trainRESP <- data[[1]]$vi
  testMAT <- (as.matrix(data[[2]][,-1]))
  testRESP <- data[[2]]$vi
  
  #
  allsp <- 1:ncol(trainMAT)
  active <- x
  
  # Fit initial regression
  ci <- ginv(trainMAT[,active]) %*% trainRESP
  e1 <- mean(((testMAT[,active] %*% ci)-testRESP)^2)/var(testRESP)
  
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
  
  #cif <- ginv(rbind(testMAT, trainMAT)[,active]) %*% c(testRESP, trainRESP)
  cif <- ginv(trainMAT[,active]) %*% trainRESP
  
  res <- matrix(0, nrow = 1, ncol = ncol(trainMAT))
  res[,active] <- cif
  
  return(res)
}

limits <- function(dat, errT, n){
  rmat <- lapply(1:ncol(dat), function(x) matrix(nrow = n, ncol = ncol(dat)))
  
  for(xi in 1:ncol(dat)){
    ds1 <- data_setup(dat, xi)
    for(i in 1:n){
      ds2 <- data_split(ds1)
      rmat[[xi]][i,] <- step_forward(ds2, xi, errT)
    }
  }
  
  imat <- sapply(rmat, function(x) apply(x, 2, median))
  intmat <- t(imat) * apply(dat, 2, median)
  
  return(intmat)
}

spec_sens <- function(real, pred){
  #diag(real) <- 0
  #diag(pred) <- 0
  t1 <- table(real,pred)
  sns <- (sum(t1[c(1,9)])-nrow(real))/(sum(t1[c(1,9,4,6)])-nrow(real))
  spc <- t1[5]/sum(t1[c(2,5,8)]) 
  return(c(specificity = spc, sensitivity = sns))
}



simfun <- function(N.species, N.links, t.fin, sub){
  cond <- FALSE
  while(!cond){
    gp1 <- get_pars(N.species, N.links, t.fin)
    dyn <- lv_ricker(1:t.fin, gp1$m, gp1$p)
    
    if(sum(is.nan(dyn)) > 0){
      cond <- FALSE
      next
    }else if(sum(dyn[,t.fin] < 10^-5) > 0){
      cond <- FALSE
      next
    }else{
      cond <- TRUE
    }
  }
  
  simdat <- t(apply(dyn, 2, function(x) x/sum(x)))
  simdat2 <- t(apply(dyn[order(apply(dyn, 1, median), decreasing = T)[1:sub],], 2, function(x) x/sum(x)))
  simdat3 <-  t(apply(rbind(dyn[order(apply(dyn, 1, median), decreasing = T)[1:sub],], 
                            colSums(dyn[-order(apply(dyn, 1, median), decreasing = T)[1:sub],])), 2, function(x) x/sum(x)))
  simdat4 <- t(dyn)
  simdat5 <- t(dyn[order(apply(dyn, 1, median), decreasing = T)[1:sub],])   
  simdat6 <- t(rbind(dyn[order(apply(dyn, 1, median), decreasing = T)[1:sub],], colSums(dyn[-order(apply(dyn, 1, median), decreasing = T)[1:sub],])))
  
  simlist <- list(simdat, simdat2, simdat3, simdat4, simdat5, simdat6)
  
  return(list(sl = simlist, par = gp1, ord = order(apply(dyn, 1, median), decreasing = T)[1:sub]))
}



get_ss <- function(lims, dat, ord){
  int.reALL <- sign(dat$p$cij)
  int.re10 <- sign(dat$p$cij[ord,ord])
  int.pred <- lapply(lims, sign)
  
  ss1 <- spec_sens(int.reALL, int.pred[[1]])
  ss2 <- spec_sens(int.re10, int.pred[[2]])
  ss3 <- spec_sens(int.re10, int.pred[[3]][1:10,1:10])
  ss4 <- spec_sens(int.reALL, int.pred[[4]])
  ss5 <- spec_sens(int.re10, int.pred[[5]])
  ss6 <- spec_sens(int.re10, int.pred[[6]][1:10, 1:10])
  
  ssall <- rbind(ss1, ss2, ss3, ss4, ss5, ss6)
  
  rownames(ssall) <- c("Allrel", "rel10", "relOth", "Allabs", "abs10", "absOth")
  
  return(ssall)
}

N <- c(12, 16)
L <- c(.1,.2,.3)
er <- rep(seq(.5, 5, .5), 3)
egdat <- (expand.grid(N, L, er))
head(egdat, 10)
dflist <- list()
for(i in 1:nrow(egdat)){
  links <- floor(egdat[i,1]^2*egdat[i,2])
  test <- simfun(egdat[i,1], links, 1000, 10)
  print("sim done")
  lims1 <- lapply(test$sl, limits, errT = egdat[i,3], 100)
  
  gs1 <- get_ss(lims1, dat = test$par, ord = test$ord)
  dflist[[i]] <- data.frame(melt(gs1), err = egdat[i,3], L = links, N = egdat[i,1])
  print(i)
}