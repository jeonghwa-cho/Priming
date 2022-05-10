## Simulation-based power analysis for Jeonghwa's English past tense priming 

#########
# Setup #
#########

#TODO: speed these up with parallel computation
# <https://cran.r-project.org/web/packages/future.apply/index.html>

setwd('/Users/jobrenn/Documents/Projects/cho-priming/')

library(tidyverse)
library(future.apply)

plan(multisession) # parallel replicate()

prop.se <- function(x) {
  p = mean(x)
  n = length(x)
  sqrt((p * (1-p))/n)
}

########
# Load #
########

data2 <- read_csv("English_past_trimmed.csv")

# merge with word frequencies

prime  <- read.csv("prime_char.csv")
target <- read.csv("target_char.csv")

data2 <- merge(data2, prime, by="prime")
data2 <- merge(data2, target, by="target")

# prep condition info

data2 <- data2 %>%
  as_tibble() %>%
  mutate(type = if_else(condition %in% c("11", "12", "13"), 
                        'Morphological',
                        'Orthographical'),
         cond = case_when(
           condition %in% c("11", "14") ~ "Identity",
           condition %in% c("12", "15") ~ "Test",
           TRUE ~ "Control" ),
         cond = factor(cond, 
                       levels=c(levels=c("Identity", "Test", "Control")))
         )

###############
# Data subset #
###############

# identity vs test, just masked morphological 

data3 <- data2 %>%
  select(number, participant, item, SOA, type, cond, rt) %>%
  filter(type == 'Orthographical',
         cond %in% c('Test', 'Control'),
         SOA == 'masked')

# # focus on data subset: masked, test-identity vs control-identity
# data3 <- data2 %>%
#   pivot_wider(names_from=cond, values_from=rt) %>%
#   mutate(dTest    = Test - Identity,
#          dControl = Control - Identity) %>%
#   select(-Test, -Control, -Identity) %>%
#   pivot_longer(c('dTest', 'dControl'), names_to='dCond', values_to='rt')

###########
# By-hand #
###########
# PRO: model different effect sizes
# CON: simplifies relationships between and within datapoints


# hyperparameters from the data
# sampling plan:
# rt = lognormal(I_a + effect, err_within) - params in log(msec)
#      where I_a = rnorm(I, err_between) - params in msec

hp <- list()

hp$eff <- data3 %>% # raw effect in msec
  group_by(participant, cond) %>%
  summarize(rt = mean(rt)) %>%
  group_by(cond) %>%
  summarize(rt = mean(rt)) %>%
  pull(rt) %>%
  diff()

hp$I <- data3 %>% # in msec
  filter(cond == 'Control') %>%
  pull(rt) %>%
  log() %>%
  mean() %>%
  exp()

hp$err_w <- data3 %>% # in log(msec)
  group_by(participant) %>%
  summarize(err = sd(log(rt))) %>%
  pull(err) %>%
  mean()

hp$err_b <- data3 %>% # in msec
  group_by(participant) %>%
  summarize(err = mean(rt)) %>%
  pull(err) %>%
  sd()

hp

# define sim functions

sim_one <- function(eff=10, Np=40, Nk=100, Nd=25, I = 700, err_b = 120, err_w=0.233) {
  sim_data = list()
  for (p in 1:Np) {
    I_a = rnorm(1, I, err_b)
    sim_data[[p]] <-
      tibble(
        item = rep(1:Nd),
        A    = rlnorm(Nd, log(I_a), err_w),
        B    = rlnorm(Nd, log(I_a + eff), err_w)
      ) %>%
      pivot_longer(cols=c(A,B), names_to = 'cond', values_to = 'rt')
  }
  sim_data <- sim_data %>% bind_rows(.id = 'ID')
  stat <- sim_data %>%
    group_by(ID, cond) %>%
    summarize(rt = mean(log(rt)), .groups='drop') %>%
    t.test(rt ~ cond, paired = TRUE, data = .)
  is_sig <- ifelse(stat$statistic < 0 & stat$p.value <= 0.05, 1, 0)
}

sim_many <- function(eff, Np, Nk, Nd, I, err_b, err_w) {
  cat('Running', Nk, 'simulations with N =', Np, 'and effect =', eff, 'ms\n')
  pwr <- future_replicate(Nk, sim_one(eff=eff, Np=Np, Nk=Nk, Nd=Nd,
                               I=I, err_b=err_b, err_w=err_w) )
  P    <- mean(pwr)
  P.se <- prop.se(pwr)
  tibble(
    Np = Np,
    Nk = Nk,
    Nd = Nd,
    E  = eff,
    P  = P,
    P.upper = P + 2 * P.se,
    P.lower = P - 2 * P.se
  )
}



sim_params <-
  expand_grid(
#    eff = c(10, 20, 30),
    #Np = c(20, 40, 60, 80, 100, 120)
    eff = c(5, 7, 10, 20),
    Np = c(20, 40, 60, 80, 100, 120, 140)
  )

# run sim!
system.time(
all_sims <- map2_dfr(sim_params$eff, 
                     sim_params$Np, 
                     ~ sim_many(eff=.x, Np=.y, Nd=25, Nk=100,
                                I=hp$I, err_b = hp$err_b, err_w=hp$err_w) )
)

save(all_sims,  file='simulated-experiments.Rda')

# plot!

all_sims %>%
  ggplot(aes(y = P, x = Np, col=factor(E), ymin=P.lower, ymax=P.upper)) +
  geom_line() +
#  geom_errorbar(col='darkgrey', width=1) +
  geom_ribbon(aes(fill=factor(E)), alpha=0.2, colour=NA) +
  geom_point() +
  scale_color_brewer('Effect Size (ms)', type='qual', palette=2) +
  scale_fill_brewer('Effect Size (ms)', type='qual', palette=2) +
  scale_x_continuous('# Participants', breaks=sim_params$Np) +
  scale_y_continuous('Estimated Power', limits=c(0, 1)) 

## TESTING PARALLELS

plan(sequential)
system.time(
sim_many(eff=10, Np=100, Nd=25, Nk=100,
         I=hp$I, err_b = hp$err_b, err_w=hp$err_w) 
)

plan(multicore)
system.time(
  sim_many(eff=10, Np=100, Nd=25, Nk=100,
           I=hp$I, err_b = hp$err_b, err_w=hp$err_w) 
)

plan(multisession)
system.time(
  sim_many(eff=10, Np=100, Nd=25, Nk=100,
           I=hp$I, err_b = hp$err_b, err_w=hp$err_w) 
)



##################
# Resample-based #
##################
# PRO: replicates relationships between datapoints
# CON: no control over effect size


sim_resample_one <- function(d, Np) {
  subjs <- sample(unique(d$participant), Np, replace=TRUE)
  sim_data <- list()
  for (p in seq_along(subjs)) {
    sim_data[[p]] <- d %>% filter(participant == subjs[p])
  }
  sim_data <- bind_rows(sim_data, .id='ID')
  stat <- sim_data %>%
    group_by(ID, cond) %>%
    summarize(rt = mean(log(rt))) %>%
    t.test(rt ~ cond, paired=TRUE, data= . )
  is_sig <- ifelse(stat$statistic < 0 & stat$p.value <= 0.05, 1, 0)
}

sim_resample_many <- function(d, Nk=100, Np=40) {
  cat('Running', Nk, 'simulations with N =', Np, '\n')
  pwr <- replicate (Nk, sim_resample_one(d=d, Np=Np) )
  tibble(
    Np = Np,
    Nk = Nk,
    P  = mean(pwr)
  )
}

Np = c(20, 40, 60, 80, 100, 120)

all_resample_sims <- map_dfr(Np, ~ sim_resample_many(data3, 100, .x) )

all_resample_sims %>%
  ggplot(aes(y = P, x = Np)) +
  geom_line() +
  geom_point() +
  ylim(0, 1) +
  scale_x_continuous(breaks=Np) 

save(all_sims, all_resample_sims,  file='simulated-experiments.Rda')


