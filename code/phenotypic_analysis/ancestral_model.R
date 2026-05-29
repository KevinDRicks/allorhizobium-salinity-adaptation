#code estimates  impact of the ancestral allorhizobium strain on plant growth across two nitrogen treatments
#code takes pixel area output from out image processing, converts it to mm2, and runs a model to characterize growth 

#load libraries
library(lme4)
library(lmerTest)
library(broom)
library(clipr)
library(broom.mixed)

#load in data
plant_growth <- read.csv("/path/to/ancestor_strain_plant_growth.csv")
#convert plate number of factor
plant_growth$PlateNumber <- as.factor(plant_growth$PlateNumber)

#area of plant is in pixels. We convert this to mm2 using our conversion factor
cf <- 211/72722.3
plant_growth$area <- plant_growth$area*cf



#linear model describing impact of ancestral allorhizobium strain and contemporary nitrogen on plant growth
growth_model <- lmer(area~Date.merge*Strain*Nitrogen+ Date.merge*StartingSize + 
                       (1|PlatexWell),
                     plant_growth)
#plate was excluded as a random effect to prevent boundary warning
anova(growth_model)
summary(growth_model)

tidy_model <- tidy(growth_model)
#write_clip(tidy_model)
