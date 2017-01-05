
d1 <- data_process(m3l2, level = "phylum")
d1.5 <- d1[25:50, ]
ma <- apply(d1.5, 2, mean)

cc1 <- ccrepe::ccrepe(d1.5)

cc2 <- cc1$sim.score 
cc2[cc1$p.values > 0.05] <- 0 
cc2[is.na(cc2)] <- 0

plot(graph.adjacency(cc2, weighted = T), layout = layout.circle, 
     edge.width = abs(E(graph.adjacency(cc2, weighted = T))$weight)*10, 
     edge.color = factor(sign(E(graph.adjacency(cc2, weighted = T))$weight)))


tw1 <- seq(1,(nrow(d1)-20), 10)
g <- list()
for(i in 1:length(tw1)){
  d1.5 <- d1[tw1[i]:(tw1[i]+49),]
  
  cc1 <- ccrepe::ccrepe(d1.5)
  
  cc2 <- cc1$sim.score 
  cc2[cc1$p.values > 0.05] <- 0 
  cc2[is.na(cc2)] <- 0
  
  g[[i]] <- graph.adjacency(cc2, weighted = T)
  
  print(i)
}
lapply(g, V)

par(mfrow = c(4,8), mar = c(.1,.1,.1,.1))
for(i in 1:length(g)){plot(g[[i]], layout = layout.circle, edge.width = abs(E(g[[i]])$weight)*10, edge.color = factor(sign(E(g[[1]])$weight)), edge.curved = T)}