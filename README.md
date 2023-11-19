# Moderate-and-Extreme-PM-in-Spain
This repository contains a Bayesian spatio-temporal extreme value model, prediction with excursion functions and a joint model with sharing effects.

In the two R scripts of Model1 and Model1 prediction, we only include the formula for M1 (Gumbel model), one may delete the fixed prior of tail==0 to obtain a GEV model. However, note that the tail parameter is sensitive in the model, a proper prior is therefore important.
One may also use the form f(s, model = spde)+f(year, model=ar1)+ f(s, ...spde+ar1) to obatain M2 and M4.
