###BAYESIAN####

setwd('/Users/jobrenn/Documents/Projects/cho-priming/')

library(brms)
library(tidyverse)
library(bayesplot)
library(tidybayes)
library(emmeans)

data2 <- read_csv("English_past_trimmed.csv")

# merge with word frequencies

prime <- read.csv("prime_char.csv")
target <- read.csv("target_char.csv")

data2<-merge(data2, prime, by="prime")
data2<-merge(data2, target, by="target")

# prep condition info

data2 <- data2 %>%
  mutate(type = if_else(condition %in% c("11", "12", "13"), 
                        'Morphological',
                        'Orthographical')) %>%
  mutate(cond = case_when(
    condition %in% c("11", "14") ~ "Identity",
    condition %in% c("12", "15") ~ "Test",
    TRUE ~ "Control" ) ) %>%
  mutate(cond = factor(cond, 
                       levels=c(levels=c("Identity", "Test", "Control")))) 

# sum coding for categorical effects?    
#contrasts(data2$condition2) <- cbind(C1 = c(1, 0, 0), C2 = c(0, 1, 0)) 
#data2$type2<-ifelse((data2$type=="Orthographical"), -0.5, 0.5)
#data2$SOA2<-ifelse((data2$SOA=="masked"), -0.5, 0.5)

# center continuous covariates
data2 <- data2 %>%
  mutate(primelength = str_length(prime)) %>%
  mutate(primelength.c = scale(primelength, center=TRUE, scale=FALSE),
         targetlength.c = scale(wordlength, center=TRUE, scale=FALSE),
         prime_freq.c = scale(prime_freq, center=TRUE, scale=FALSE),
         target_freq.c = scale(target_freq, center=TRUE, scale=FALSE)) 

# to set reasonably weakly-informative priors  
mu=mean(log(data2$rt))
sigma=sd(log(data2$rt))

mu2=exp(mu+sigma^2/2)
sigma2=mu2*sqrt(exp(sigma^2)-1)  

#sigma: exponential distribution -exp(1)

  
  # set priors
  # see: https://vasishth.github.io/bayescogsci/book/sec-trial.html
  # and: https://vasishth.github.io/bayescogsci/book/modeling-a-lexical-decision-task.html
my_priors <- c(
    prior(normal(6, 1), class = 'Intercept'),
    prior(normal(0, 1), class = 'b'),
    prior(normal(0, 1), class = 'sigma'),
    prior(normal(0, 1), class = 'sd')
)
  
  ### SANDBOX
  subjects <- unique(data2$participant)
  
  # test

  data2 %>% 
    filter(participant %in% subjects[1:80]) %>%
    brm(rt ~ 1 + targetlength.c + (1 + targetlength.c | participant),
        data = .,
        prior = my_priors,
        cores = 4,
        iter=8000,
        family= lognormal() ) -> m_test
  
  m_test
#  pairs(m_test)

## NOTES
# <https://discourse.mc-stan.org/t/low-ess-and-high-rhat-for-random-intercept-slope-simulation-rstan-and-rstanarm/9985/2>
# - could the model be too good? no
# - I think the problem is non-independence between intercept and participant sd... and I think it can be addressed by increasing iterations;??
# - IF increased iterations keeps RHat low... it does!?
# - ELSE there is something wrong with the model...
#
# - model: rt ~ 1, subjects[1:20]: looks good!
#
# - model: rt ~ 1 + (1|participant)
# -   subjects[1:5]: ESS ~ 1000
# -   subjects[1:20]: ESS < 300 :( :(
# -   no problem with (1|item)!
# -   population priors too loose?
#  
# - model: rt ~ intercept + targetlength.c + (1|participant)
# - 20 - 40 participants (about 1 min to fit)
# - lognormal() returns more accurate fit than shifted_lognormal(), 
#     intercept = 6.5 rather than 6.0... effect of shift param??
# - regardless Intercept is poorly explored (ESS < 400; rhat >1.01)  
# - ESS *lower* with more participants
  

    
  
summary(m_test)
plot(m_test)
pp_check(m_test)
mcmc_plot(m_test)
mcmc_plot(m_test, "^b_[^I]")

### END SANDBOX

##################
# simple version #
##################
# 
# # fit model
# f2<-brm(rt ~ cond * type * SOA + primeLength.c + targetLength.c + 
#           (1 | participant) + 
#           (1 | item), 
#         data=data2, 
#         prior = my_priors,
#         cores=4,
#         iter = 8000,
#         family=lognormal())

# fit model (frequency included)
f2<-brm(rt ~ cond * type * SOA * prime_freq.c * target_freq.c + primeLength.c + targetLength.c + 
          (1 | participant) + 
          (1 | item), 
        data=data2, 
        prior = my_priors,
        cores=4,
        iter = 8000,
        family=lognormal())

summary(f2)
plot(f2)
pp_check(f2)
mcmc_plot(f2)
mcmc_plot(f2, "^b_[^I]")

save(f2, file='model-simple.Rda')
#load(file='model-simple.Rda')

# plot estimated effects

color_scheme_set('red')
mcmc_areas_ridges(f2, regex_pars="^b_[^I]")



# Test specific hypotheses

# - orthographic, masked : test > identity
H = ' (typeOrthographical + condTest + condTest:typeOrthographical) >
      (typeOrthographical)'
hypothesis(f2, H)

# - masked, test : ortho = morph priming?
H = '(condTest + typeOrthographical + condTest:typeOrthographical) = 
      condTest'
hypothesis(f2, H)

# - masked, control : ortho = morph priming?
H = '(condControl + typeOrthographical + condControl:typeOrthographical) = 
      condControl'
hypothesis(f2, H)

# - masked, test-control : ortho = morph priming?
H = '(condTest + typeOrthographical + condTest:typeOrthographical) - 
     (condControl + typeOrthographical + condControl:typeOrthographical) = 
      (condTest - condControl)'
hypothesis(f2, H)


# - unmasked, test - control : ortho:prime.freq > morph:prime.freq

H = '(SOAunmasked + condTest + typeOrthographical + prime_freq.c + 
    condTest:SOAunmasked + 
    condTest:typeOrthographical + 
    condTest:prime_freq.c + 
    condTest:typeOrthographical:SOAunmasked + 
    condTest:typeOrthographical:prime_freq.c + 
    condTest:SOAunmasked:prime_freq.c +
    typeOrthographical:SOAunmasked:prime_freq.c +
    condTest:typeOrthographical:SOAunmasked:prime_freq.c) -
    (SOAunmasked + condControl + typeOrthographical + prime_freq.c + 
    condControl:SOAunmasked + 
    condControl:typeOrthographical + 
    condControl:prime_freq.c + 
    condControl:typeOrthographical:SOAunmasked + 
    condControl:typeOrthographical:prime_freq.c + 
    condControl:SOAunmasked:prime_freq.c +
    typeOrthographical:SOAunmasked:prime_freq.c +
    condControl:typeOrthographical:SOAunmasked:prime_freq.c) >
    (SOAunmasked + condTest + prime_freq.c + 
    condTest:SOAunmasked + 
    condTest:prime_freq.c + 
    condTest:SOAunmasked:prime_freq.c) -
    (SOAunmasked + condControl + prime_freq.c + 
    condControl:SOAunmasked + 
    condControl:prime_freq.c + 
    condControl:SOAunmasked:prime_freq.c) '
hypothesis(f2, H)


# Many ways to plot expected values across conditions!

# (1)  with conditional_effects
c <- make_conditions(data2, vars=c('SOA'))
conditional_effects(f2, 
                    effects='type:cond',
                    conditions = c)

# (2) - with emmeans... 
# same as 1!
my_emm <- f2 %>% emmeans(~ cond + type + SOA, epred=TRUE) 

my_emm %>%
  summary() %>%
  #  gather_emmeans_draws(ndraws=100)
  ggplot(aes(x=type, y=emmean, ymin=lower.HPD, ymax=upper.HPD, col=cond)) +
  geom_pointinterval(position='dodge') +
  facet_wrap(~SOA) # why do these have less var than conditional_effects?

# (3) - emmeans but now with draws 
# (same again!)
my_emm %>% gather_emmeans_draws(ndraws=500) %>%
  ggplot(aes(x=type, y=.value, col=cond)) +
    stat_pointinterval(position='dodge') +
    facet_wrap(~SOA) # why do these have less var than conditional_effects?

# (4) - fitted values for every data point DON'T DO THIS
# less var! prog not marginalizing correctly!!
my_draws <- data2 %>% add_epred_draws(f2, ndraws=100, re_formula)

my_draws %>%
  group_by(type, cond, SOA, .draw) %>%
  summarize(yhat = mean(.epred)) %>%
  ggplot(aes(x=type, y=yhat, col=cond)) +
    stat_pointinterval(position='dodge') +
    facet_wrap(~SOA) # why do these have less var than conditional_effects?

# (5) - fitted values for target conditions all else NA
# same as conditional_effects!
my_draws_2 <-
  expand_grid(
    type = unique(data2$type),
    cond = unique(data2$cond),
    SOA = unique(data2$SOA),
    targetLength.c = 0,
    primeLength.c = 0,
  ) %>%
  add_epred_draws(f2, ndraws=500, re_formula = NA)

my_draws_2 %>%
  ggplot(aes(x=type, y=.epred, col=cond)) +
  stat_pointinterval(position='dodge') +
  facet_wrap(~SOA) # why do these have MORE var than conditional_effects?


# TAKE-AWAYS: conditional_effects (1), emmeans (2), emmeans+tidybayes (3) and just tidybayes with covariates set to zero/NA (5) are all equivalent and look good; DO NOT average tidybayes draws from every datapoint!!! (4)

#########################
#full version (updated) #
#########################

f3 <- brm(rt ~ cond * type * SOA + 
             primeLength.c + targetLength.c + 
             (1 + cond * type *SOA || participant) + 
             (1 + cond * type *SOA || item),
           data=data2,
           prior = my_priors,
           iter = 8000, 
           family=lognormal())

summary(f3)
plot(f3)
pp_check(f3)

save(f3, file='model-complex.Rda')
