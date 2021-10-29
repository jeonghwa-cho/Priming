####concatenate all data####

library(dplyr)

setwd("C:/Users/LG/Downloads/SpanishRaw")

(temp = list.files(pattern="*.csv"))

O<-list()

for (i in 1:length(temp)){
  file<-read.csv(temp[i])
  
  number<-rep(i, 200)
  participant<-file[32:231, "participant"]
  rt<-file[32:231,"key_resp_7.rt"]
  condition<-file[32:231, "condition"]
  item<-file[32:231, "item"]
  prime<-file[32:231, "prime"]
  target<-file[32:231, "target"]
  list <- file[32:231, "expName"]
  key<-file[32:231, "key_resp_7.keys"]
  corr<-file[32:231, "key_resp_7.corr"]
  
  df<-data.frame(number, participant, rt, item, prime, target, condition, list, key, corr)
  O[[i]]<-df
}

data<-bind_rows(O)

write.csv(data, "Spanish_masked_data.csv")

library(stringr)
data<-read.csv("C:/Users/LG/Downloads/SpanishRaw/Spanish_masked_data.csv")
data$exp<-"masked"

umdata<-read.csv("C:/Users/LG/Downloads/umSpanishRaw/Spanish_unmasked_data.csv")
umdata$exp<-"unmasked"

data<-rbind(data, umdata)

data$participant2<-paste(data$participant, data$list, sep="_")
data$participant<-as.factor(data$participant2)

write.csv(data, "C:/Users/LG/Downloads/SpanishRaw/Spanish_data_all.csv")

####analysis####

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

sum<-data3%>%
  group_by(condition)%>%
  summarise(rt=mean(rt))

library(ggplot2)
data3$LexIsSame <- ifelse((data3$condition == "11"|data3$condition == "12"), "LexId", "LexDiff")
data3$TenseIsSame <- ifelse((data3$condition == "11"|data3$condition == "13"), "TenseId", "TenseDiff")

sum<-data3%>%
  group_by(LexIsSame, TenseIsSame, participant)%>%
  summarise(rt=mean(rt))

sum$Tense <- factor(sum$TenseIsSame, levels=c("TenseId", "TenseDiff"))
levels(sum$Tense)

sum%>%ggplot(aes(x=LexIsSame, y=rt, fill=Tense))+
  geom_bar(aes(fill=Tense), stat="summary",fun=mean, width=0.8,position=position_dodge(0.9))+
  geom_errorbar(stat='summary', fun.data=mean_se,width=.5,position=position_dodge(0.9))+
  scale_x_discrete(limits=c("LexId", "LexDiff"))+
  ylab("RT (ms)")+
  xlab("")+
  theme_bw()

sum%>%ggplot(aes(x=condition, y=rt))+
  geom_bar(stat="identity",position=position_dodge())

####LMER####

library(lmerTest)
library(MuMIn)
data3$wordlength<-str_length(data3$target)
data3$primelength<-str_length(data3$prime)
data3$primeLength.c = scale(data3$primelength, center=TRUE, scale=FALSE)
data3$targetLength.c = scale(data3$wordlength, center=TRUE, scale=FALSE)

data3$lex<-ifelse((data3$condition=="11"|data3$condition=="12"), -0.5, 0.5)
data3$tense<-ifelse((data3$condition=="11"|data3$condition=="13"), -0.5, 0.5)

freq <- read.csv("C:/Users/LG/Desktop/Priming_Psychopy/stimuli info/freq.csv")
data3$condition_item<-paste(data3$condition, data3$item, sep="_")
freq$condition_item<-paste(freq$condition, freq$item, sep="_")

data4<-merge(data3, freq, by="condition_item")
data4$primefreq.c = scale(data4$primefreq, center=TRUE, scale=FALSE)
data4$targetfreq.c = scale(data4$targetfreq, center=TRUE, scale=FALSE)

summary(data4)

f<-lmer(log(rt)~lex*tense*primefreq.c+primeLength.c+targetLength.c+(1|participant)+(1|item.x), data=data4)
summary(f)

f<-lmer(log(rt)~condition+primeLength.c+targetLength.c+(1|participant)+(1|item), data=data3%>%filter(condition=="13"|condition=="14"))
summary(f)

####Bayesian####

library(brms)
library(tidyverse)
f1<-brm(rt~lex*tense*primefreq.c*targetfreq.c+primeLength.c+targetLength.c+(1+lex*tense|participant)+(1+lex*tense|item.x), 
        data=data4%>%filter(data4$exp=="masked"), 
        prior = c(
          prior(normal(6, 1.5), class = Intercept),
          prior(normal(0, 1), class = b)
        ),
        family=shifted_lognormal())



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