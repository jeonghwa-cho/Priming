####analysis####
data<-read.csv("Spanish_data_all.csv")

data$condition<-as.factor(data$condition)
data$item<-as.factor(data$item)
data$rt<-data$rt*1000

data<-data%>%filter(condition=="11"|condition=="12"|condition=="13"|condition=="14")

accuracy<-data%>%dplyr::group_by(participant)%>%dplyr::summarise(accuracy=sum(corr)/100*100)
View(accuracy)
data<-data[data$participant!="60d346f1edaf18ae50b94e28_SPA_um_listB",]
data<-data[data$participant!="6075b914e3df3a75e8a5fab8_SPA_um_listC",]
data<-data[data$participant!="5f49d522d8bdca55b61c293c_SPA_um_listB",]
data<-data[data$participant!="5f40382cedd4d214a0dbdd7f_SPA_um_listC",]

length(unique(data$participant))
       
accuracy<-data%>%dplyr::group_by(item)%>%dplyr::summarise(accuracy=sum(corr)/98*100)
View(accuracy)
data<-data[data$item!="32",]

data2<-data[data$corr=="1",]
data2<-data

data3<-data2[data2$rt>200,]
data3<-data3[data3$rt<2000,] #96.9% 
nrow(data3)/nrow(data2)

####Bayesian####

library(brms)
library(tidyverse)

data3$wordlength<-str_length(data3$target)
data3$primelength<-str_length(data3$prime)
data3$primeLength.c = scale(data3$primelength, center=TRUE, scale=FALSE)
data3$targetLength.c = scale(data3$wordlength, center=TRUE, scale=FALSE)

data3$lex<-ifelse((data3$condition=="11"|data3$condition=="12"), -0.5, 0.5)
data3$tense<-ifelse((data3$condition=="11"|data3$condition=="13"), -0.5, 0.5)


f2<-brm(rt~lex*tense+primeLength.c+targetLength.c+(1|participant)+(1|item), 
        data=data3%>%filter(data3$exp=="masked"), 
        prior = c(
          prior(normal(6, 1.5), class = Intercept),
          prior(normal(0, 1), class = b)
        ),
        family=shifted_lognormal())
pp_check(f2)
plot(f2)
summary(f2)
mcmc_plot(f2)
mcmc_plot(f2, "^b_[^I]")

f3<-brm(rt~lex*tense+primeLength.c+targetLength.c+(1+lex*tense|participant)+(1+lex*tense|item), 
        data=data3%>%filter(data3$exp=="unmasked"), 
        prior = c(
          prior(normal(6, 1.5), class = Intercept),
          prior(normal(0, 1), class = b)
        ),
        family=shifted_lognormal())

pp_check(f3)
plot(f3)
summary(f3)
mcmc_plot(f3)
mcmc_plot(f3, "^b_[^I]")

cat( "Dataset is as follows; ",
         "\n", 
         "\n", 
         file = "C:/Users/LG/Desktop/brm_fit_priming.txt") 

capture.output(summary(f2), 
                file = "C:/Users/LG/Desktop/brm_fit_priming.txt", 
                append = TRUE) 

################
