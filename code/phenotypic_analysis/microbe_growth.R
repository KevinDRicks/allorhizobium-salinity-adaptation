#code pulls in microbial growth rates and characterizes differences between selective history

#load libraries
library("lme4")
library('lmerTest')
library('broom')
library('clipr')
library('broom.mixed')
microbe_growth <- read.csv("/path/to/microbial_growth_rates.csv")

#currently data has a column for the growth rate in each environment
#we can run models on each of these
summary(lm(Microbegrowth_r_LowN_LowS~HistN*HistS,microbe_growth))
summary(lm(Microbegrowth_r_LowN_HighS~HistN*HistS,microbe_growth))
summary(lm(Microbegrowth_r_HighN_LowS~HistN*HistS,microbe_growth))
summary(lm(Microbegrowth_r_HighN_HighS~HistN*HistS,microbe_growth))

strain.rep <- nrow(microbe_growth)

#we can also run it in the long form, using all the data at the same time in a full model, using invididual strain as a random effect
#this requires restructuring the data, adding in new columns for each contemporary environment
merge.data <- rbind(data.frame(ContSaline=rep('Low',strain.rep),ContNitrogen=rep('Low',strain.rep),
                               r=microbe_growth$Microbegrowth_r_LowN_LowS,Strain=microbe_growth$StrainID,
                               HistSaline=microbe_growth$HistS,HistNitrogen=microbe_growth$HistN),
                    data.frame(ContSaline=rep('High',strain.rep),ContNitrogen=rep('Low',strain.rep),
                               r=microbe_growth$Microbegrowth_r_LowN_HighS,Strain=microbe_growth$StrainID,
                               HistSaline=microbe_growth$HistS,HistNitrogen=microbe_growth$HistN),
                    data.frame(ContSaline=rep('Low',strain.rep),ContNitrogen=rep('High',strain.rep),
                               r=microbe_growth$Microbegrowth_r_HighN_LowS,Strain=microbe_growth$StrainID,
                               HistSaline=microbe_growth$HistS,HistNitrogen=microbe_growth$HistN),
                    data.frame(ContSaline=rep('High',strain.rep),ContNitrogen=rep('High',strain.rep),
                               r=microbe_growth$Microbegrowth_r_HighN_HighS,Strain=microbe_growth$StrainID,
                               HistSaline=microbe_growth$HistS,HistNitrogen=microbe_growth$HistN))
                    
microbial_growth_model <- (lmer(r~ContSaline*ContNitrogen*HistSaline*HistNitrogen + (1|Strain),merge.data))
anova(microbial_growth_model)
summary(microbial_growth_model)

tidy_model <- tidy(microbial_growth_model)
#write_clip(tidy_model)
