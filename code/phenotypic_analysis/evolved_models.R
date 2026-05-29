#code estimates  impact of the evolved allorhizobium strains on plant growth across salinity environments
#code takes pixel area output from out image processing, converts it to mm2, and runs a model to characterize growth 
#data includes sterile and ancestral inoculation, which we exclude for this analysis

#load libraries
library(lme4)
library(lmerTest)
library(broom)
library(clipr)
library(broom.mixed)

#load in data
plant_growth <- read.csv("/path/to/evolved_strain_plant_growth.csv")
plant_growth$Plate <- as.factor(plant_growth$Plate)

#data list area of plant in pixels. We convert this to mm2 
cf <- 211/72722.3
#convert plate number of factor
plant_growth$area <- plant_growth$area*cf

#for the full model, we only considered the evolved strains
plant_growth_evolved <- plant_growth[plant_growth$Status=='Evolved',]

#linear model describing impact of ancestral allorhizobium strain and contemporary nitrogen on plant growth
growth_model <- lmer(area~Date.merge*ContSaline*SalineHistory*NitrogenHistory+Date.merge*StartingSize +  (1|StrainID) + 
               (1|Chamber_Column) + (1|PlatexWell),
               plant_growth_evolved)



anova(growth_model)
summary(growth_model)

tidy_model <- tidy(growth_model)
#write_clip(tidy_model)


#modeling individual strain effects across salinity environments
#coefficients extracted from this model represent the individual strain growth estimates, seen in the strain_partner_quality.csv 
#specifically the strain x time x salinity coefficient represents plant growth rate in each salinity environment with each strain
strain_model <- lmer(area~Date.merge*StrainID*ContSaline+Date.merge*StartingSize +
                       (1|Chamber_Column) + (1|PlatexWell),
                     plant_growth)
anova(strain_model)
summary(strain_model)
