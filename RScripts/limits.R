# Functions for the LIMITS (Learning Interactions through MIcrobial Time Series) Algorithm
library(magrittr)

spp <- 6
fitmod <- r.m3l2A %>% data_setup(xi = spp) %>% data_split %>% step_forward(x = spp, errorthres = 1)
fitmod

data_setup <- function(data, xi){
  # take in matrix with columns = taxa, rows = timepoints
  # make sure we are dealing with relative abundances
  if(sum(round(rowSums(data), 2) != 1.00) > 0){stop("Needs relative abundance")}
  
  colnames(data) <- sapply(strsplit(colnames(data), ";"), "[[", 2)
  vi <- log(data[2:nrow(data), xi]) - log(data[1:(nrow(data)-1), xi])
  M <- apply(data, 2, function(x) x - median(x))[1:(nrow(data)-1),]
  
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
  err <- mean((data$test$vi-predict(fit.init, newdata = data$test))^2)/var(data$train$vi)
  
  active <- x
  
  while(!cond){
    form <- sapply((1:length(xnam))[-active], function(y){
      as.formula(paste("vi ~ ", paste0(c(xnam[active], xnam[y]), collapse = "+")))
    })
    
    fit.step <- lapply(form, function(y){
      fit.step <- lm(y, data = data$train)
    })
    
    err.step <- sapply(fit.step, function(y){
      mean((data$test$vi-predict(y, newdata = data$test))^2)/var(data$test$vi)
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
  
  form.fin <- as.formula(paste("vi ~ ", paste0(c(xnam[active], xnam[y]), collapse = "+")))
  fit.fin <- lm(form.fin, data = data$train)
  return(fit.fin)
}


iter <- 200
spp <- 2
comat <- matrix(0, nrow = iter, ncol = ncol(data$train)-1)
for(i in 1:iter){
  fitmod <- r.m3l2A %>% data_setup(xi = spp) %>% data_split %>%
    step_forward(x = spp, errorthres = 1)
  comat[i,colnames(data$train)[-1]%in%names(fitmod$coefficients[-1])] <- fitmod$coefficients[-1]
}
colMeans(comat)

c.vec <- function(idata, spp, iter = 200, e.thres = 1){
  comat <- matrix(0, nrow = iter, ncol = ncol(data$train)-1)
  for(i in 1:iter){
    fitmod <- idata %>% data_setup(xi = spp) %>% data_split %>%
      step_forward(x = spp, errorthres = e.thres)
    comat[i,colnames(data$train)[-1]%in%names(fitmod$coefficients[-1])] <- fitmod$coefficients[-1]
  }
  return(colMeans(comat))
}

imat <- matrix(0, 6, 6)
for(i in 1:6){
  imat[i,] <- c.vec(r.m3l2A, spp = i, iter = 200, e.thres = .1)
}
imat 


taxa1 <- c("Akkermansia muciniphilia", "AAlistipes putredinis", "Bacteroides acidifaciens", "Bacteroides fragilis", "Bacteroides stercoris", "Bacteroides thetaiotaomicron", "Bacteroides uniformis", "Bacteroides vulgatus", "Escherichia coli", "Eubacterium rectale", "Faecalibacterium prausnitzii", "Parebacteroides distasonis", "Roseburia intestinalis", "Sisymbrium irio")
