# Functions for the LIMITS (Learning Interactions through MIcrobial Time Series) Algorithm
library(magrittr)

spp <- 6
fitmod <- r.m3l2A %>% data_setup(xi = spp) %>% data_split %>% step_forward(x = spp, errorthres = 1)
fitmod

data_setup <- function(data, xi){
  # take in matrix with columns = taxa, rows = timepoints
  # make sure we are dealing with relative abundances
  if(sum(round(rowSums(data), 2) != 1.00) > 0){stop("Needs relative abundance")}
  
  colnames(data) <- sapply(strsplit(colnames(data), ";"), tail, 1)
  
  resp <- data[data[,xi] != 0, xi]
  
  vi <- log(resp[2:length(resp)]) - log(resp[1:(length(resp)-1)])
  M <- apply(data[data[,xi] != 0,], 2, function(x){x - median(x)})[1:(nrow(data[data[,xi] != 0,])-1),]
  
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


iter <- 200
spp <- 7
comat <- matrix(0, nrow = iter, ncol = ncol(r.m3l2A))
for(i in 1:iter){
  data <- r.m3l2A %>% data_setup(xi = spp) %>% data_split
  fitmod <- step_forward(data = data, x = spp, errorthres = 1)
  comat[i,colnames(data$train)[-1]%in%names(fitmod$coefficients[-1])] <- fitmod$coefficients[-1]
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
  return(colMeans(comat))
}

imat <- matrix(0, ncol(r.m3l2A), ncol(r.m3l2A))
for(i in 1:ncol(r.m3l2A)){
  imat[i,] <- c.vec(r.m3l2A, spp = i, iter = 200, e.thres = 3)
}
imat 



m3l6 <- as.matrix(read.csv("Data/M3L6gut.csv", row.names = 1))

r.m3l6 <- m3l6[which(apply(m3l6, 1, function(x) sum(x !=0)) >= 332),]
r.m3l6A <- t(apply(r.m3l6, 2, function(x) x/sum(x)))


imat <- matrix(0, ncol(r.m3l6A), ncol(r.m3l6A))
for(i in 1:ncol(r.m3l6A)){
  imat[i,] <- c.vec(r.m3l6A, spp = i, iter = 200, e.thres = 3)
}
imat 

taxa1 <- c("Akkermansia muciniphilia", "AAlistipes putredinis", "Bacteroides acidifaciens", "Bacteroides fragilis", "Bacteroides stercoris", "Bacteroides thetaiotaomicron", "Bacteroides uniformis", "Bacteroides vulgatus", "Escherichia coli", "Eubacterium rectale", "Faecalibacterium prausnitzii", "Parebacteroides distasonis", "Roseburia intestinalis", "Sisymbrium irio")
