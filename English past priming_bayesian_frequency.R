library(brms)
library(stringr)
library(tidyverse)

###BAYESIAN####

data2<-read.csv("English_past_trimmed.csv")

prime <- read.csv("prime_char.csv")
target <- read.csv("target_char.csv")

data2<-merge(data2, prime, by="prime")
data2<-merge(data2, target, by="target")

data2$type <- ifelse((data2$condition == "11"|data2$condition == "12"|data2$condition=="13"), "Morphological", "Orthographical")
data2$condition2 <- ifelse((data2$condition == "11"|data2$condition == "14"), "Identity", ifelse((data2$condition == "12"|data2$condition == "15"), "Test", "Control"))
data2$condition2<- factor(data2$condition2, levels=c("Identity", "Test", "Control"))

contrasts(data2$condition2) <- cbind(C1 = c(1, 0, 0), C2 = c(0, 1, 0)) 
data2$type2<-ifelse((data2$type=="Orthographical"), -0.5, 0.5)
data2$SOA2<-ifelse((data2$SOA=="masked"), -0.5, 0.5)

data2$primelength<-str_length(data2$prime)
data2$primeLength.c = scale(data2$primelength, center=TRUE, scale=FALSE)
data2$targetLength.c = scale(data2$wordlength, center=TRUE, scale=FALSE)

data2$prime_freq.c = scale(data2$prime_freq, center=TRUE, scale=FALSE)
data2$target_freq.c = scale(data2$target_freq, center=TRUE, scale=FALSE)
data2$prime_neigh.c = scale(data2$prime_neighbors, center=TRUE, scale=FALSE)
data2$target_neigh.c = scale(data2$target_neighbors, center=TRUE, scale=FALSE)

mu=mean(log(data2$rt))
sigma=sd(log(data2$rt))

mu2=exp(mu+sigma^2/2)
sigma2=mu2*sqrt(exp(sigma^2)-1)  

#sigma: exponential distribution -exp(1)

#masked
masked<-data2%>%filter(SOA=="masked")

m1<-brm(rt~condition2*type2*prime_freq.c*target_freq.c+prime_neigh.c+target_neigh.c+primeLength.c+targetLength.c+(1+condition2*type2||participant)+(1+condition2*type2||item),
           data=masked,
           prior = c(
             prior(normal(6, 1.5), class = Intercept),
             prior(normal(0, 1), class = b)
           ),
           iter = 3000, 
           warmup = 2000,
           control = list(adapt_delta = 0.99),
           family=shifted_lognormal())

#unmasked
unmasked<-data2%>%filter(SOA=="unmasked")

m2<-brm(rt~condition2*type2*SOA2+primeLength.c+targetLength.c+(1+condition2*type2*SOA2||participant)+(1+condition2*type2*SOA2||item),
        data=unmasked,
        prior = c(
          prior(normal(6, 1.5), class = Intercept),
          prior(normal(0, 1), class = b)
        ),
        iter = 3000, 
        warmup = 2000,
        control = list(adapt_delta = 0.99),
        family=shifted_lognormal())

#subset
length(unique(data2$item))

sub <- data2%>%filter(participant %in% unique(participant)[80:120],
                      item %in% unique(item)[30:100])

f2<-brm(rt~condition2*type2*SOA2+primeLength.c+targetLength.c+(1|participant)+(1|item), 
        data=sub, 
        prior = c(
          prior(normal(6, 1.5), class = Intercept),
          prior(normal(0, 1), class = b)
        ),
        family=shifted_lognormal())


summary(f2)
plot(f2)
pp_check(f2)
mcmc_plot(f2)
mcmc_plot(f2, "^b_[^I]")

f3<-brm(rt~condition2*type2*SOA2+primeLength.c+targetLength.c+(1|participant)+(1|item), data=data2)
summary(f3)
plot(f3)
pp_check(f3)
