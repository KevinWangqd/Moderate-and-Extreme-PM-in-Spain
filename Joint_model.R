library(INLA)
library(ggplot2)
library(lattice)     
library(lwgeom)
library(dplyr)
library(leaflet)
library(viridis)
library(rgdal)

## Get map data
m <- getData(name = "GADM", country = "Spain", level = 0)
m <- m %>%
  st_as_sf() %>%
  st_cast("POLYGON") %>%
  mutate(area = st_area(.)) %>%
  arrange(desc(area)) %>%
  slice(1)
ggplot(m) + geom_sf() + theme_bw()
m <- m %>% st_transform(25830)
ggplot(m) + geom_sf() + theme_bw() + coord_sf(datum = st_crs(m))

## Read Data.csv
d <- read.csv("Dataset.csv")

ggplot(m) + geom_sf() + coord_sf(datum = st_crs(m)) +
  geom_point(data = d, aes(x = x, y = y)) + theme_bw()


##Mesh construction and SPDE
coo <- cbind(d$x, d$y)
bnd <- inla.nonconvex.hull(st_coordinates(m)[, 1:2])
mesh <- inla.mesh.2d(
  loc = coo, boundary = bnd,
  max.edge = c(100000, 200000), cutoff = 1000
)
plot(mesh)
points(coo, col = "Red")

spde <- inla.spde2.pcmatern(
  mesh = mesh, alpha = 2, constr = TRUE,
  prior.range = c(10000, 0.01), # P(range < 10000) = 0.01
  prior.sigma = c(3, 0.01) # P(sigma > 3) = 0.01
)

## index: train and validation; two indexes for mean and maxima model
index <- d$year
train_index <- which(d$year <= 2020)
test_index <- which(d$year == 2021)


timesn <- length(unique(d$year))
indexs1 <- inla.spde.make.index("s1",
                               n.spde = spde$n.spde,
                               n.group = timesn
)
indexs2 <- inla.spde.make.index("s2",
                                n.spde = spde$n.spde,
                                n.group = timesn
)


## Projection matrix A
group <- d$year - min(d$year) + 1
A.est <- inla.spde.make.A(mesh = mesh, loc = coo[train_index,], group = group[train_index],n.group=timesn)
A.val <- inla.spde.make.A(mesh = mesh, loc = coo[test_index,], group = group[test_index],n.group=timesn)


## Stack for joint model (estimation and validation for mean and maxima model)

stack.est.mean <- inla.stack(data = list(y = cbind(as.vector(d$logmean[train_index]), NA)),
                     A = list(1,A.est), 
                     effects = list(data.frame(Intercept1 = rep(1, length(train_index)),
                                               long1=d$longs[train_index],lat1=d$lats[train_index],
                                               alt1=d$alts[train_index], year1=d$year[train_index], 
                                               meanTem1=d$Ts[train_index], meanPre1 = d$Ps[train_index], 
                                               VP1=d$VPs[train_index], pop1 = d$pds[train_index]),s1 = indexs1),
                     tag = "est.mean") 

stack.est.max <- inla.stack(data = list(y = cbind(NA, d$logmax[train_index])),
                        A = list(1,A.est), 
                        effects = list(data.frame(Intercept2 = rep(1, length(train_index)),
                                                  long2=d$longs[train_index],lat2=d$lats[train_index],
                                                  alt2=d$alts[train_index], year2=d$year[train_index], 
                                                  meanTem2=d$Ts[train_index], meanPre2 = d$Ps[train_index], 
                                                  VP2=d$VPs[train_index], pop2 = d$pds[train_index]),s2 = indexs2),
                        tag = "est.max") 




stack.val.mean <- inla.stack(data = list(y = cbind(rep(NA, length(d$logmean[test_index])), NA)),
                             A = list(1,A.val), 
                             effects = list(data.frame(Intercept1 = rep(1, length(test_index)),
                                                       long1=d$longs[test_index],lat1=d$lats[test_index],
                                                       alt1=d$alts[test_index], year1=d$year[test_index], 
                                                       meanTem1=d$Ts[test_index], meanPre1 = d$Ps[test_index], 
                                                       VP1=d$VPs[test_index], pop1 = d$pds[test_index]),s1 = indexs1),
                             tag = "val.mean") 


stack.val.max <- inla.stack(data = list(y = cbind(rep(NA, length(d$logmax[test_index])), NA)),
                             A = list(1,A.val), 
                             effects = list(data.frame(Intercept2 = rep(1, length(test_index)),
                                                       long2=d$longs[test_index],lat2=d$lats[test_index],
                                                       alt2=d$alts[test_index], year2=d$year[test_index], 
                                                       meanTem2=d$Ts[test_index], meanPre2 = d$Ps[test_index], 
                                                       VP2=d$VPs[test_index],  pop2 = d$pds[test_index]),s2 = indexs2),
                             tag = "val.max") 



stk.full <- inla.stack.join(stack.est.mean, stack.est.max,
                            stack.val.mean,stack.val.max)


## Prior 
rprior <- list(theta = list(prior = "pccor1", param = c(0, 0.9))) ## ar1 prior
hyper <- list(beta = list(prior = 'normal', param = c(0,10))) ## vague Gaussian prior for sharing coefficients
hyper.gev = list(tail = list(initial = 0,fixed = TRUE)) ## Gumbel likelihood



## construct sharing fixed effect index
index1 <- stack.est.mean$effects$nrow
index2 <- stack.est.max$effects$nrow
val_ind <- stack.val.max$effects$nrow + stack.val.max$effects$nrow

idx.alt1 = c(rep(1,index1), rep(NA,index2+val_ind))
idx.temp1 = c(rep(1,index1), rep(NA,index2+val_ind))
idx.Pre1 = c(rep(1,index1), rep(NA,index2+val_ind))
idx.VP1 = c(rep(1,index1), rep(NA,index2+val_ind))
idx.pop1 = c(rep(1,index1), rep(NA,index2+val_ind))

idx.alt2 = c(rep(NA,index1), rep(1,index2), rep(NA,val_ind))
idx.temp2 = c(rep(NA,index1), rep(1,index2), rep(NA,val_ind))
idx.Pre2 = c(rep(NA,index1), rep(1,index2), rep(NA,val_ind))
idx.VP2 = c(rep(NA,index1), rep(1,index2), rep(NA,val_ind))
idx.pop2 = c(rep(NA,index1), rep(1,index2), rep(NA,val_ind))



## formula
formula <- 0 +
  Intercept1 + Intercept2 +
  long1 + lat1 + long2 + lat2 +
  f(idx.alt1, alt1, model = "iid") +
  f(idx.alt2, alt2, copy = "idx.alt1", fixed = FALSE, hyper = hyper) +
  f(idx.VP1, VP1, model = "iid") +
  f(idx.VP2, VP2, copy = "idx.VP1", fixed = FALSE, hyper = hyper) +
  f(idx.pop1, pop1, model = "iid") +
  f(idx.pop2, pop2, copy = "idx.pop1", fixed = FALSE, hyper = hyper) +
  f(idx.temp1, meanTem1, model = "iid") +
  f(idx.temp2, meanTem2, copy = "idx.temp1", fixed = FALSE, hyper = hyper) +
  f(idx.Pre1, meanPre1, model = "iid") +
  f(idx.Pre2, meanPre2, copy = "idx.Pre1", fixed = FALSE, hyper = hyper) +
  f(idx.VP1, VP1, model = "iid") +
  f(idx.VP2, VP2, copy = "idx.VP1", fixed = FALSE, hyper = hyper) +
  f(idx.pop1, pop1, model = "iid") +
  f(idx.pop2, pop2, copy = "idx.pop1", fixed = FALSE, hyper = hyper) +
  f(s1, model = spde, group = s1.group, control.group = list(model = "ar1", hyper = rprior)) +
  f(s2, copy = "s1", fixed = FALSE, hyper = hyper)



## joint model
res_Gaussian_Gumbel <- inla(
  formula,
  family = c("gaussian", "gev"),
  control.family = list(list(), list(hyper = hyper.gev)),
  data = inla.stack.data(stk.full),
  control.predictor = list(link = 1, compute = TRUE, A = inla.stack.A(stk.full)),
  control.inla = list(int.strategy = "eb"),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE, return.marginals.predictor = TRUE),
  verbose = TRUE
)

summary(res_Gaussian_Gumbel)
hist(res_Gaussian_Gumbel$cpo$cpo)
hist(res_Gaussian_Gumbel$cpo$pit)


## posterior distributions 
res_Gaussian_Gumbel$summary.random$idx.temp1
res_Gaussian_Gumbel$summary.random$idx.Pre1
res_Gaussian_Gumbel$summary.random$idx.VP1
res_Gaussian_Gumbel$summary.random$idx.alt1
res_Gaussian_Gumbel$summary.random$idx.pop1

res_Gaussian_Gumbel$summary.random$idx.temp2
res_Gaussian_Gumbel$summary.random$idx.Pre2
res_Gaussian_Gumbel$summary.random$idx.VP2
res_Gaussian_Gumbel$summary.random$idx.alt2
res_Gaussian_Gumbel$summary.random$idx.pop2








## Sharing coefficients

Beta <- list(
 "Temperature" = res_Gaussian_Gumbel$marginals.hyperpar$`Beta for idx.temp2`,
  "Precipitation" = res_Gaussian_Gumbel$marginals.hyperpar$`Beta for idx.Pre2`,
  "Vapour Pressure" = res_Gaussian_Gumbel$marginals.hyperpar$`Beta for idx.VP2`,
  "Populatin Density" = res_Gaussian_Gumbel$marginals.hyperpar$`Beta for idx.pop2`,
  "Altitude" = res_Gaussian_Gumbel$marginals.hyperpar$`Beta for idx.alt2`
)

Beta_r <- data.frame(do.call(rbind, Beta)) 
Beta_r$parameter <- rep(names(Beta), times = sapply(Beta, nrow))

ggplot(Beta_r, aes(x = x, y = y)) +
  geom_line() +
  facet_wrap(~parameter, scales = "free") +
  labs(x = "", y = "Density") +
  theme_bw()


## Diamond plot
res_Gaussian_Gumbel$summary.hyperpar[11,]
res_Gaussian_Gumbel$summary.hyperpar[12,]
res_Gaussian_Gumbel$summary.hyperpar[13,]
res_Gaussian_Gumbel$summary.hyperpar[14,]
res_Gaussian_Gumbel$summary.hyperpar[15,]

par(mfrow=c(1,5))

plot(res_Gaussian_Gumbel$summary.hyperpar$`0.5quant`[11], ylim = c(-1, 1),
     ylab = "", xlab="", xaxt = "n")
title("Altitude")
points(res_Gaussian_Gumbel$summary.hyperpar$`0.025quant`[11], 
       col = "blue", pch = 18, cex = 1.4)    
points(res_Gaussian_Gumbel$summary.hyperpar$`0.975quant`[11], 
       col = "blue", pch = 18, cex = 1.4)    
segments(1, res_Gaussian_Gumbel$summary.hyperpar$`0.025quant`[11],  
         1, res_Gaussian_Gumbel$summary.hyperpar$`0.975quant`[11],
         lwd = 2, lty = 3, col = "blue")
abline(h = 0, lty = 2, lwd = 2, col = "gray")


plot(res_Gaussian_Gumbel$summary.hyperpar$`0.5quant`[13], ylim = c(-1, 1), 
     ylab = "", xlab="", xaxt = "n")
title("Population density")
points(res_Gaussian_Gumbel$summary.hyperpar$`0.025quant`[13], 
       col = "blue", pch = 18, cex = 1.4)    
points(res_Gaussian_Gumbel$summary.hyperpar$`0.975quant`[13], 
       col = "blue", pch = 18, cex = 1.4)    
segments(1, res_Gaussian_Gumbel$summary.hyperpar$`0.025quant`[13],  
         1, res_Gaussian_Gumbel$summary.hyperpar$`0.975quant`[13], 
         lwd = 2, lty = 3, col = "blue")
abline(h = 0, lty = 2, lwd = 2, col = "gray")

plot(res_Gaussian_Gumbel$summary.hyperpar$`0.5quant`[15], ylim = c(-1, 1), 
     ylab = "", xlab="", xaxt = "n")
title("Precipitation")
points(res_Gaussian_Gumbel$summary.hyperpar$`0.025quant`[15], 
       col = "blue", pch = 18, cex = 1.4)    
points(res_Gaussian_Gumbel$summary.hyperpar$`0.975quant`[15], 
       col = "blue", pch = 18, cex = 1.4)    
segments(1, res_Gaussian_Gumbel$summary.hyperpar$`0.025quant`[15],  
         1, res_Gaussian_Gumbel$summary.hyperpar$`0.975quant`[15], 
         lwd = 2, lty = 3, col = "blue")
abline(h = 0, lty = 2, lwd = 2, col = "gray")



plot(res_Gaussian_Gumbel$summary.hyperpar$`0.5quant`[14], ylim = c(-1, 1), 
     ylab = "", xlab="", xaxt = "n")
title("Temperature")
points(res_Gaussian_Gumbel$summary.hyperpar$`0.025quant`[14], 
       col = "blue", pch = 18, cex = 1.4)    
points(res_Gaussian_Gumbel$summary.hyperpar$`0.975quant`[14], 
       col = "blue", pch = 18, cex = 1.4)    
segments(1, res_Gaussian_Gumbel$summary.hyperpar$`0.025quant`[14],  
         1, res_Gaussian_Gumbel$summary.hyperpar$`0.975quant`[14], 
         lwd = 2, lty = 3, col = "blue")
abline(h = 0, lty = 2, lwd = 2, col = "gray")


plot(res_Gaussian_Gumbel$summary.hyperpar$`0.5quant`[12], ylim = c(-1, 1), 
     ylab = "", xlab="", xaxt = "n")
title("Vapour pressure")
points(res_Gaussian_Gumbel$summary.hyperpar$`0.025quant`[12], 
       col = "blue", pch = 18, cex = 1.4)    
points(res_Gaussian_Gumbel$summary.hyperpar$`0.975quant`[12], 
       col = "blue", pch = 18, cex = 1.4)    
segments(1, res_Gaussian_Gumbel$summary.hyperpar$`0.025quant`[12],  
         1, res_Gaussian_Gumbel$summary.hyperpar$`0.975quant`[12], 
         lwd = 2, lty = 3, col = "blue")
abline(h = 0, lty = 2, lwd = 2, col = "gray")

dev.off()









