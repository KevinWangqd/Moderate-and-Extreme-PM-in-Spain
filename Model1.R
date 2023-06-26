## Packages
library(INLA)
library(raster)
library(ggplot2)
library(viridis)
library(lattice)
library(gridExtra)
library(viridis)

## Functions (slcpo is the same as Logarithm score)
RMSE <- function(set, outcome, data, fit) {
  res <- data[set, outcome] - fit[set]
  RMSE_val <- sqrt(mean(res^2, na.rm = TRUE))
  return(RMSE_val)
}

slcpo <- function(m, na.rm = TRUE) {
  -sum(log(m$cpo$cpo), na.rm = na.rm)
}

## Get map data
m <- getData(name = "GADM", country = "Spain", level = 0)  # Spain territory
plot(m)

m <- m %>%  # Cut the mainland map
  st_as_sf() %>%
  st_cast("POLYGON") %>%
  mutate(area = st_area(.)) %>%
  arrange(desc(area)) %>%
  slice(1)

ggplot(m) +
  geom_sf() +
  theme_bw()

m <- m %>% st_transform(25830)
ggplot(m) +
  geom_sf() +
  theme_bw() +
  coord_sf(datum = st_crs(m))

## PM and predictors data
d <- read.csv("Dataset.csv")
ggplot(m) +
  geom_sf() +
  coord_sf(datum = st_crs(m)) +
  geom_point(data = d, aes(x = x, y = y)) +
  theme_bw()

## Mesh construction
coo <- cbind(d$x, d$y)
bnd <- inla.nonconvex.hull(st_coordinates(m)[, 1:2])
mesh <- inla.mesh.2d(
  loc = coo, boundary = bnd,
  max.edge = c(100000, 200000), cutoff = 1000
)
points(coo, col = "red")

## Data description and plot
ggplot(d) +
  geom_histogram(mapping = aes(x = value)) +
  facet_wrap(~year, ncol = 1) +
  theme_bw()

ggplot(m) +
  geom_sf() +
  coord_sf(datum = NA) +
  geom_point(
    data = d, aes(x = x, y = y, color = value),
    size = 2
  ) +
  labs(x = "", y = "") +
  scale_color_viridis() +
  facet_wrap(~year) +
  theme_bw()

ggplot(d, aes(x = year, y = value, group = id, color = id)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = c(2017, 2018, 2019, 2020, 2021)) +
  theme_bw() +
  theme(legend.position = "none")

## Index: train and validation
index <- d$year
train_index <- which(d$year <= 2020)
test_index <- which(d$year == 2021)

## SPDE (PC prior)
spde <- inla.spde2.pcmatern(
  mesh = mesh, alpha = 2, constr = TRUE,
  prior.range = c(10000, 0.01),  # P(range < 10000) = 0.01
  prior.sigma = c(3, 0.01)  # P(sigma > 3) = 0.01
)

timesn <- length(unique(d$year))
indexs <- inla.spde.make.index("s",
                               n.spde = spde$n.spde,
                               n.group = timesn
)
lengths(indexs)


## Projection matrix A (training) and Ap (prediction/validation)

group <- d$year - min(d$year) + 1
A <- inla.spde.make.A(mesh = mesh, loc = coo[train_index, ], group = group[train_index], n.group = timesn)
Ap <- inla.spde.make.A(mesh = mesh, loc = coo[test_index, ], group = group[test_index], n.group = timesn)
dim(A)


## Stack (estimation and prediction/validation)
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
stk.full <- inla.stack(stk.e, stk.p)


## Prior 

rprior <- list(theta = list(prior = "pccor1", param = c(0, 0.9)))  # ar1 prior

hyper.gev = list(tail = list(initial = 0, fixed = TRUE))  # fix Gumbel likelihood; \xi==0


## Formula: fixed effect of covariates + spatio-temporal random effect

formula1 <- y ~ 0 + b0 + long + lat + alt + VP + meanTem + meanPre + pop + 
  f(s, model = spde, group = s.group, control.group = list(model = "ar1", hyper = rprior))

## The following formula is for Model 2/4
## formula2 <- y ~ 0 + b0 + long + lat + alt + VP + meanTem + meanPre + pop + f(year, model = "ar1", hyper = rprior) + f(s, model = spde)


# Model 1             
res_Gumbel1 <- inla(formula1,
                    data = inla.stack.data(stk.full), 
                    family = "gev",
                    control.family = list(hyper = hyper.gev),
                    control.predictor = list(A = inla.stack.A(stk.full), link = 1, compute = TRUE),
                    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                    verbose = TRUE
)

# Model 1             
res_Gumbel1 <- inla(formula1,
                    data = inla.stack.data(stk.full), 
                    family = "gev",
                    control.family = list(hyper = hyper.gev),
                    control.predictor = list(A = inla.stack.A(stk.full), link = 1, compute = TRUE),
                    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
                    verbose = TRUE
)


## Model outputs
summary(res_Gumbel1)
res_Gumbel1$summary.hyperpar

## Result processing and data preparation for evaluation
index_train_inla <- inla.stack.index(stk.full, "pre")$data
index_val_inla <- inla.stack.index(stk.full, "tst")$data

result.train1 <- res_Gumbel1$summary.fitted$mean[index_train_inla]
result.val1 <- res_Gumbel1$summary.fitted$mean[index_val_inla]

M_fit_spde1 <- array(NA, nrow(d))
M_fit_spde1[train_index] <- result.train1
M_fit_spde1[test_index] <- result.val1



# Model Comparison (We only include the code for Model 1. For Model 2/3/4, the user may change the likelihood and the formula)

## DIC and WAIC
res_Gumbel1$dic$dic
res_Gumbel1$waic$waic

## RMSE
rmse_train_spde1 <- RMSE(train_index, "logmax", d, M_fit_spde1)
rmse_val_spde1 <- RMSE(test_index, "logmax", d, M_fit_spde1)
rmse_train_spde1
rmse_val_spde1

## CPO and PIT
hist(res_Gumbel1$cpo$pit, xlab = "", main = "Histogram of Model 1 PIT")
cpo_1 <- slcpo(res_Gumbel1, na.rm = TRUE)
cpo_1

## Coverage probability
pre_table_train1 <- as.data.frame(cbind(d$logmax[index_train_inla], 
                                        res_Gumbel1$summary.fitted$`0.025quant`[index_train_inla],
                                        res_Gumbel1$summary.fitted$`0.975quant`[index_train_inla],
                                        res_Gumbel1$summary.fitted$mean[index_train_inla]))
pre_table_val1 <- as.data.frame(cbind(d$logmax[index_val_inla], 
                                      res_Gumbel1$summary.fitted$`0.025quant`[index_val_inla],
                                      res_Gumbel1$summary.fitted$`0.975quant`[index_val_inla],
                                      res_Gumbel1$summary.fitted$mean[index_val_inla]))

cp_train1 <- pre_table_train1[which(pre_table_train1$V1 >= pre_table_train1$V2 & pre_table_train1$V1 <= pre_table_train1$V3), ]
cp_val1 <- pre_table_val1[which(pre_table_val1$V1 >= pre_table_val1$V2 & pre_table_val1$V1 <= pre_table_val1$V3), ]

coverage_prob_train1 <- length(cp_train1[, 1]) / length(pre_table_train1[, 1])
coverage_prob_val1 <- length(cp_val1[, 1]) / length(pre_table_val1[, 1])


## Visualization in training and validation sets
plotting_data1 <- data.frame(Observed = d$logmax,
                             Predicted = M_fit_spde1)
plotting_data1$group[train_index] <- "train"
plotting_data1$group[test_index] <- "validation" 
plotting_data1$group[train_index] <- "train"

p1 <- ggplot(data = plotting_data1, aes(x = Observed, y = Predicted)) +
  geom_point(aes(color = group)) + geom_abline(intercept = 0, slope = 1) + xlim(2, 7) +
  ylim(2, 7) + theme_bw() + ggtitle("")
p1

## Correlation
obs_val <- plotting_data1[which(plotting_data1$group=="validation"),1]
pre_val <- plotting_data1[which(plotting_data1$group=="validation"),2]
scatter.smooth(obs_val,pre_val)
cor(obs_val,pre_val) 

