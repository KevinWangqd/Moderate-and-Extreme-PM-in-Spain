## Library
library(lwgeom)
library(readxl)
library(INLA)
library(raster)
library(ggplot2)
library(reshape2)
library(viridis)
library(sp)
library(lattice)     
library(gridExtra)
library(sf)
library(dplyr)
library(excursions)
library(viridis)
library(RColorBrewer)

## The following part is the same as M1.R, except an additional projection matrix A.pre

m <- getData(name = "GADM", country = "Spain", level = 0)
plot(m)
m <- m %>%
  st_as_sf() %>%
  st_cast("POLYGON") %>%
  mutate(area = st_area(.)) %>%
  arrange(desc(area)) %>%
  slice(1)
ggplot(m) + geom_sf() + theme_bw()
m <- m %>% st_transform(25830)
d <- read.csv("Dataset.csv")

## Mesh construction
coo <- cbind(d$x, d$y)
bnd <- inla.nonconvex.hull(st_coordinates(m)[, 1:2])
mesh <- inla.mesh.2d(
  loc = coo, boundary = bnd,
  max.edge = c(100000, 200000), cutoff = 1000
)
points(coo, col = "red")

## Index: train and validation
index <- d$year
train_index <- which(d$year <= 2020)
test_index <- which(d$year == 2021)

## SPDE (PC prior)
spde <- inla.spde2.pcmatern(
  mesh = mesh, alpha = 2, constr = TRUE,
  prior.range = c(10000, 0.01), # P(range < 10000) = 0.01
  prior.sigma = c(3, 0.01) # P(sigma > 3) = 0.01
)

timesn <- length(unique(d$year))
indexs <- inla.spde.make.index("s",
                               n.spde = spde$n.spde,
                               n.group = timesn
)
lengths(indexs)


## Projection matrix A (training) and Ap (prediction/validation)
group <- d$year - min(d$year) + 1
A <- inla.spde.make.A(mesh = mesh, loc = coo[train_index,], group = group[train_index],n.group=timesn)
Ap <- inla.spde.make.A(mesh = mesh, loc = coo[test_index,], group = group[test_index],n.group=timesn)
dim(A)


## Prediction projection matrix A.pre
dp <- read.csv("dpred.csv")
coop <- dp[, c("x", "y")]
coop <- as.matrix(coop)
A.pre <- inla.spde.make.A(mesh = mesh, loc = coop, group = 1,n.group=timesn)

## Stack
stk.e <- inla.stack(
  tag = "pre",
  data = list(y = d$logmax[train_index]),
  A = list(1, A),
  effects = list(data.frame(b0 = rep(1, length(train_index)),
                            long = d$longs[train_index], lat = d$lats[train_index], alt = d$alts[train_index],
                            year = d$year[train_index], VP = d$VPs[train_index], meanTem = d$Ts[train_index], 
                            meanPre = d$Ps[train_index], pop = d$pds[train_index]), s = indexs)
)

stk.p <- inla.stack(
  tag = "tst",
  data = list(y = NA),
  A = list(1, Ap),
  effects = list(data.frame(b0 = rep(1, length(test_index)),
                            long = d$longs[test_index], lat = d$lats[test_index], alt = d$alts[test_index],
                            year = d$year[test_index], VP = d$VPs[test_index], meanTem = d$Ts[test_index], 
                            meanPre = d$Ps[test_index], pop = d$pds[test_index]), s = indexs)
)



stk.pre <- inla.stack(
  tag = "pred",
  data = list(y = NA),
  A = list(1, A.pre),
  effects = list(data.frame(b0 = rep(1, length(dp$longs)), 
                            long = dp$longs, lat = dp$lats, alt = dp$alts,
                            meanTem= dp$Ts, meanPre=dp$Ps,pop=dp$pds, VP = dp$VPs),s = indexs)
)

stk.full <- inla.stack(stk.e, stk.p, stk.pre)

## Prior 

rprior <- list(
  theta = list(prior = "pccor1", param = c(0, 0.9)),
  theta1  = list(prior = "pcprec", param = c(5, 0.01))
)

hyper.gev = list(
  tail = list(initial = 0, fixed = TRUE),
  prec  = list(prior = "pcprec", param = c(1, 0.01))
)


## Formula
formula1 <- y ~ 0 + b0 + long + lat + alt + VP + meanTem + meanPre + pop +
  f(s,model = spde, group = s.group,control.group = list(model = "ar1", hyper = rprior)
  )


## Gumbel1 model             
res_Gumbel1 <- inla(formula1,
                    data = inla.stack.data(stk.full), 
                    family = "gev",
                    control.family = list(hyper = hyper.gev),
                    control.predictor=list(A=inla.stack.A(stk.full), link = 1, compute=TRUE),
                    control.compute = list(dic = TRUE, waic=TRUE, cpo = TRUE, config=TRUE,return.marginals.predictor=TRUE),
                    verbose = TRUE
)    


## Prediction (excursion functions)
cols <- rev(colorRampPalette(brewer.pal(6,"RdYlGn"))(20))

##positive excursions
excursion_100_pos <- excursions.inla(res_Gumbel1,stack = stk.full,tag="pred",alpha = 0.05,u=log(100),type = ">",method = "QC",F.limit = 1)
excursion_50_pos <- excursions.inla(res_Gumbel1,stack = stk.full,tag="pred",alpha = 0.05,u=log(50),type = ">",method = "QC",F.limit = 1)

##negative excursions
excursion_100_neg <- excursions.inla(res_Gumbel1,stack = stk.full,tag="pred",alpha = 0.05,u=log(100),type = "<",method = "QC",F.limit = 1)
excursion_50_neg <- excursions.inla(res_Gumbel1,stack = stk.full,tag="pred",alpha = 0.05,u=log(50),type = "<",method = "QC",F.limit = 1)


dp$excursion_100_pos <- excursion_100_pos$F
dp$excursion_50_pos <- excursion_50_pos$F
dp$excursion_100_neg <- excursion_100_neg$F
dp$excursion_50_neg <- excursion_50_neg$F

#write.csv(dp,"dp.csv")

### Here we cut the mainland with boundaries for autonomous community map
map <- getData(name = "GADM", country = "Spain", level = 1)
plot(map)
map1 <- map %>%
  st_as_sf() %>%
  st_cast("POLYGON") %>%
  mutate(area = st_area(.)) %>%
  arrange(desc(area)) %>%
  slice(1:15)
ggplot(map1) + geom_sf() + theme_bw()


## plots
ggplot(map1) + geom_sf() + coord_sf(datum = NA) +
  geom_point(
    data = dp, aes(x = long, y = lat, color = excursion_100_pos),
    size = 2
  ) +
  labs(x = "", y = "", colour = "Ex.Function(+)") +
  scale_colour_gradientn(limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1), colors=c("white","yellow","orange","red"))+
  theme_bw()

ggplot(map1) + geom_sf() + coord_sf(datum = NA) +
  geom_point(
    data = dp, aes(x = long, y = lat, color = excursion_50_pos),
    size = 3
  ) +
  labs(x = "", y = "",colour = "Ex.Function(+)") +
  scale_colour_gradientn(limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1), colors=c("white","yellow","orange","red"))+
  theme_bw()


ggplot(map1) + geom_sf() + coord_sf(datum = NA) +
  geom_point(
    data = dp, aes(x = long, y = lat, color = excursion_100_neg),
    size = 3
  ) +
  labs(x = "", y = "",colour = "Ex.Function(-)") +
  scale_colour_gradientn(limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1), colors=c("white","green","cyan"))+
  theme_bw()

ggplot(map1) + geom_sf() + coord_sf(datum = NA) +
  geom_point(
    data = dp, aes(x = long, y = lat, color = excursion_50_neg),
    size = 3
  ) +
  labs(x = "", y = "",colour = "Ex.Function(-)") +
  scale_colour_gradientn(limits=c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1), colors=c("white","green","cyan"))+
  theme_bw()





## The following is about the plot for spatial effect

spatial <- read.csv("spatial.csv")
coord <- spatial[,4:5]
coord <- as.matrix(coord)
coord <- SpatialPoints(coord, proj4string = CRS(as.character(NA)),
                       bbox = NULL)

gproj <- inla.mesh.projector(mesh,  coord)

g.mean <- inla.mesh.project(gproj, res_Gumbel1$summary.random$s$mean[1:mesh$n], projection="longlat")
g.sd <- inla.mesh.project(gproj, res_Gumbel1$summary.random$s$sd[1:mesh$n])


a <- as.data.frame(g.mean)
b <- as.data.frame(g.sd)

spatial$gmean <- a$g.mean
spatial$gsd <- b$g.sd

ggplot(map1) + geom_sf() + coord_sf(datum = NA) +
  geom_point(
    data = spatial, aes(x = long, y = lat, color = gmean),
    size = 4
  )  +
  labs(x = "", y = "",colour = "") +
  scale_color_viridis(option="inferno",begin = 0.5) +  theme_bw()+ggtitle("Mean")

ggplot(map1) + geom_sf() + coord_sf(datum = NA) +
  geom_point(
    data = spatial, aes(x = long, y = lat, color = gsd),
    size = 4
  ) + 
  labs(x = "", y = "",colour = "") +scale_color_viridis()+theme_bw()+ggtitle("Standard Deviation")
