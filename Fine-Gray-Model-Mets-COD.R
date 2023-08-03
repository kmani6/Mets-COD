setwd("~/Desktop/result/FinalResults")
rm(list=ls())

library(MASS)
library(tidyverse)
library(survival)
library(prodlim)
library(riskRegression)
library(ggpubr)
library(QHScrnomo)
library(keras)
library(magrittr)
library(dplyr)
library(xgboost)
library(pROC)
library(caret)

#Code for fitting Fine-Gray competing risks and developing competing risk nomograms, 
#as adapted from the following two papers:

#Zhang, Z., Cortese, G., Combescure, C., Marshall, R., Lee, M., Lim, H. J., Haller, B., 
#& written on behalf of AME Big-Data Clinical Trial Collaborative Group (2018). 
#Overview of model validation for survival regression model with competing risks using melanoma study data. 
#Annals of translational medicine, 6(16), 325. https://doi.org/10.21037/atm.2018.07.38

#Zhang, Z., Geskus, R. B., Kattan, M. W., Zhang, H., & Liu, T. (2017). 
#Nomogram for survival analysis in the presence of competing risks. 
#Annals of translational medicine, 5(20), 403. https://doi.org/10.21037/atm.2017.07.27

#A SEER Case Listing file may be generated from 1992 to 2019 with inclusion criteria
#of all patients that have metastatic cancer. This file cannot be provided as SEER
#has National Cancer Institute privacy agreements that users must adhere to. 

#data = read.csv("MetsCOD1992.csv")

# Assign new column names
colnames(data) <- c("year", "month", "event", "agecat","sex","race","origin","site","bone","brain","liver","lung","TStage3","NStage3","TStage2","NStage2","TStage1","NStage1","siteCOD","otherCOD")

# Combine different AJCC Staging Systems from 2010 to 2019.
data$tstage <- ifelse(data$year >= 2010 & data$year <= 2015, data$TStage1,
                      ifelse(data$year >= 2016 & data$year <= 2017, data$TStage2,
                             ifelse(data$year >= 2018 & data$year <= 2019, data$TStage3, NA)))
data <- data[!(data$tstage %in% c(NA, "Blank(s)", "Not applicable", "88")), ]

data$tstage <- gsub("^.(1.*)", "T1", data$tstage)
data$tstage <- gsub("^.(2.*)", "T2", data$tstage)
data$tstage <- gsub("^.(3.*)", "T3", data$tstage)
data$tstage <- gsub("^.(4.*)", "T4", data$tstage)
data$tstage <- gsub("^.(X.*)", "TX", data$tstage)
data$tstage <- gsub("^.(0.*)", "T0", data$tstage)

data <- data[!(data$tstage %in% c("Tis", "pIS", "Ta", "Tis(LAMN)")), ]

data$nstage <- ifelse(data$year >= 2010 & data$year <= 2015, data$NStage1,
                      ifelse(data$year >= 2016 & data$year <= 2017, data$NStage2,
                             ifelse(data$year >= 2018 & data$year <= 2019, data$NStage3, NA)))

data <- data[!(data$nstage %in% c(NA, "Blank(s)", "Not applicable", "88")), ]

data$nstage <- gsub("^.(1.*)", "N1", data$nstage)
data$nstage <- gsub("^.(2.*)", "N2", data$nstage)
data$nstage <- gsub("^.(3.*)", "N3", data$nstage)
data$nstage <- gsub("^.(4.*)", "N4", data$nstage)
data$nstage <- gsub("^.(X.*)", "NX", data$nstage)
data$nstage <- gsub("^.(0.*)", "N0", data$nstage)

# Change values in the "event" column
data$event[data$event == "Alive"] <- 0
data$event[data$siteCOD == "Dead (attributable to this cancer dx)"] <- 1
data$event[data$otherCOD == "Dead (attributable to causes other than this cancer dx)"] <- 2
data <- data[data$event != "Dead", ]

# Remove "years" from agecat
data$agecat <- gsub(" years", "", data$agecat)

# Create a function to categorize the age ranges
categorize_age <- function(age_range) {
  if (age_range == "85+") {
    return("85+")
  } else if (age_range %in% c("00", "01-04", "05-09", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54")) {
    return("0-54")
  } else if (age_range %in% c("55-59", "60-64")) {
    return("55-64")
  } else if (age_range %in% c("65-69", "70-74")) {
    return("65-74")
  } else if (age_range %in% c("75-79", "80-84")) {
    return("75-84")
  } else {
    return(NA)
  }
}

# Apply the categorize_age function to agecat column
data$agecat <- sapply(data$agecat, categorize_age)

# Remove "Unknown" from months and convert to numeric
data <- subset(data, month != "Unknown")
data$month <- as.numeric(data$month)

# Modify values in bone, brain, liver, and lung columns
data <- data[data$bone != "Unknown" & data$brain != "Unknown" & data$liver != "Unknown" & data$lung != "Unknown", ]
data$bone <- ifelse(data$bone == "Yes", 1, 0)
data$brain <- ifelse(data$brain == "Yes", 1, 0)
data$liver <- ifelse(data$liver == "Yes", 1, 0)
data$lung <- ifelse(data$lung == "Yes", 1, 0)

#Combine different Sites into correct SEER groupings.
data$site <- ifelse(data$site %in% c("NHL - Extranodal", "NHL - Nodal"),
                    "NHL", data$site)

data$site <- ifelse(data$site %in% c("Cecum", "Ascending Colon", "Hepatic Flexure",
                                     "Transverse Colon", "Splenic Flexure", "Descending Colon",
                                     "Sigmoid Colon", "Large Intestine, NOS"),
                    "Colon", data$site)

data$site <- ifelse(data$site %in% c("Rectosigmoid Junction", "Rectum"),
                    "Rectum", data$site)

data$site <- ifelse(data$site %in% c("Liver", "Intrahepatic Bile Duct"),
                    "Liver and Bile Duct", data$site)

data$site <- ifelse(data$site %in% c("Lip", "Tongue", "Salivary Gland", "Floor of Mouth",
                                     "Gum and Other Mouth", "Nasopharyx", "Tonsil",
                                     "Oropharynx", "Hypopharynx", "Other Oral Cavity and Pharynx"),
                    "Oral Cavity and Pharynx", data$site)

# Find the top 15 sites
top_sites <- head(sort(table(data$site), decreasing = TRUE), 15)

# Filter the data to include only the rows with the top 15 sites
data <- data[data$site %in% names(top_sites), ]

# select cols
dt = data %>%
  dplyr::select(year, month, event,
                agecat, sex, race, site,
                bone, brain, liver, lung,
                tstage, nstage)

dt <- dt %>%
  mutate(tstage = plyr::revalue(tstage, 
                                c("T0" = "T0/TX",
                                  "TX" = "T0/TX")),
         nstage = plyr::revalue(nstage, 
                                c("N0" = "N0/NX",
                                  "NX" = "N0/NX"))) %>%
  mutate(agecat = factor(agecat,
                         c("0-54", "55-64", "65-74", 
                           "75-84", "85+")),
         sex = factor(sex, c("Male", "Female")),
         site = factor(site, names(sort(table(dt$site), T))),
         tstage = factor(tstage, c("T0/TX", "T1", "T2",
                                   "T3", "T4")),
         nstage = factor(nstage, levels = c("N0/NX","N1", "N2", "N3")),
         agecat = factor(agecat, levels = c("0-54", "55-64", "65-74", "75-84", "85+")))

# Now "Breast" is the reference level for the "site" variable
dt$site <- factor(dt$site, levels = c("Breast", "Lung and Bronchus", "Ovary", "Urinary Bladder", "Pancreas", "Esophagus", "Prostate", "Kidney and Renal Pelvis", "Colon", "Oral Cavity and Pharynx", "Melanoma of the Skin", "Stomach", "Rectum", "Liver and Bile Duct", "Corpus Uteri"))

# Remove "Unknown" from race
dt <- subset(dt, race != "Unknown")
# Remove rows with NA in the "race" column
dt <- dt[complete.cases(dt$race), ]

dt$race <- factor(dt$race, levels = c("White", "Black", "Asian or Pacific Islander", "American Indian/Alaska Native"))

# Set the seed for reproducibility
set.seed(125)

# Split the data into training and validation sets, using 80/20 split
train_indices <- createDataPartition(dt$month, p = 0.8, list = FALSE)
dt_train <- dt[train_indices, ]
dt_test <- dt[-train_indices, ]

dt_train$race <- recode(dt_train$race, "Asian or Pacific Islander" = "AAPI")
dt_train$race <- recode(dt_train$race, "American Indian/Alaska Native" = "AIAN")
dt_train$site <- recode(dt_train$site, "Lung and Bronchus" = "Lung")
dt_train$site <- recode(dt_train$site, "Kidney and Renal Pelvis" = "Kidney")
dt_train$site <- recode(dt_train$site, "Oral Cavity and Pharynx" = "Oral Cavity")
dt_train$site <- recode(dt_train$site, "Urinary Bladder" = "Bladder")
dt_train$site <- recode(dt_train$site, "Melanoma of the Skin" = "Melanoma")
dt_train$site <- recode(dt_train$site, "Liver and Bile Duct" = "Liver")

#The cause may either be 1 for deadh due to diangosed cancer or 2 for death due to
#non-cancer causes.

FGModel = FGR(Hist(month, event) ~ agecat + sex + race + site +
                bone + brain + liver + lung + tstage + nstage,
              cause = 1, data = dt_train)

FGScore <- Score(list("FGR" = FGModel), 
                 Hist(month, event)~1, data=dt, 
                 cause=1, times=c(1, 3, 5)*12, 
                 plots="ROC", metrics="auc")

#Code for cross-validation method for Fine-Gray competing risks as modeled from
#Zhang et. al "Overview of model validation for survival regression model with 
#competing risks using melanoma study data (Annals of translational medicine, 2018)

FGScore.CV <- Score(list("FGR" = FGModel),
                    Hist(month,event)~1, 
                    data=dt_train, cause=1, times =c(1, 3, 5)*12,
                    split.method = "cv5",
                    B = 1,
                    parallel = c("multicore"),
                    plots = c("ROC", "calibration"))

# Perform cross-validation without explicitly fitting the model
score.cv <- Score(list("FGR" = FGR(Hist(month, event) ~ agecat + sex + race + site +
                                     bone + brain + liver + lung + tstage + nstage,
                                   cause = 1, data = dt)),
                  Hist(month, event) ~ 1,
                  data = dt, cause = 1, times = c(1, 3, 5) * 12,
                  split.method = "cv5",
                  B = 1,
                  parallel = c("multicore"),
                  plots = "calibration")

# Print the results
print(auc)

#Calculate ROC curves for Full Model
FGScore(s2, times = 12)
FGScore(s2, times = 36)
FGScore(s2, times = 60)
title(main = "AUROC of death due to other causes at 1-year",line=1.0)

# Calculate the minimum and maximum values of the 'month' variable
min_month <- min(dt_train$month)
max_month <- max(dt_train$month)

# Generate a sequence of time points covering the range of 'month' variable if needed.
times <- seq(min_month, max_month, by = 5) # Adjust the 'by' argument as needed

# Generate Calibration curves of the full model.
score.cv <- riskRegression::Score(list("Fine-Gray" = FGModel),
                                  formula = Hist(month, event) ~ 1,
                                  cause=1,
                                  data = dt_train,
                                  times = c(1, 3, 5) * 12,
                                  cens.method = "local",
                                  se.fit = 1L,
                                  plots = "cal")

# Plot AUC as a function of time. 
ggplot(data=score.cv$AUC$score,aes(x=times,y=AUC,colour=model))+
  geom_point()+geom_line()

#Change Lines 240-246 to plot results as a function of time. 
times <-  c(1) * 12 # Adjust the time points as needed
par(mar = c(5, 6, 4, 2))  # Adjust the values inside the c() function to increase/decrease the margins
plotCalibration(score.cv, times = times, cens.method = "jacknife", bar = TRUE, col= c("red","blue"), main = "5-year", xlab = "", ylab = "", xlim = c(0,1), ylim = c(0,1))
mtext(side=2, line=4.0, "Estimated actual risk", col="black", font=1,cex=1.2)
mtext(side=1, line=2.5, "Predicted risk", col="black", font=1,cex=1.2)
title(main = "Calibration at 1-year")
1-year-Cause-line

library(pROC)
library(rms)

# Extract the coefficient information from FGModel manually
coefficients <- FGModel$crrFit$coef
exp_coefficients <- exp(coefficients)

# Calculate approximate 95% CI for HRs
lower_CI <- exp(coefficients - 1.96 * coef_summary[, "se(coef)"])
upper_CI <- exp(coefficients + 1.96 * coef_summary[, "se(coef)"])

# Combine the information into a data frame
coef_table <- data.frame(
  Predictor = names(exp_coefficients),
  Hazard_Ratio = round(exp_coefficients, 2)
)

# Generate a Markdown table
md_table <- knitr::kable(coef_table, format = "markdown")

# Write the Markdown table to a text file
writeLines(md_table, "table.md")

#CODE to develop Nomograms for survival analysis in the presence of competing risks
#as adapted from the following paper:
#Zhang, Z., Geskus, R. B., Kattan, M. W., Zhang, H., & Liu, T. (2017). 
#Nomogram for survival analysis in the presence of competing risks. 
#Annals of translational medicine, 5(20), 403. https://doi.org/10.21037/atm.2017.07.27

library(MASS)
library(tidyverse)
library(survival)
library(mstate)
library(rms)
library(QHScrnomo)

# nomogram for cause 1
dt_train = dt_train %>%
  mutate(month = ifelse(month == 0, 0.05, month),
         id = 1:nrow(dt_train))

covariate = colnames(dt_train)[4:13]

df.w = crprep("month", "event",
              data=dt_train, trans=c(1,2),
              cens=0, id="id",
              keep=covariate)

ddist <- datadist(df.w)
options(datadist='ddist')

#Code for COD labeled as diagnosed cancer.
FineGrayMod <- cph(Surv(Tstart,Tstop,status==1)~agecat + sex + 
              race + site + bone +
              brain + liver + lung + 
              tstage + nstage,
            data=df.w,
            weight=weight.cens,
            subset=failcode==1,
            surv=T)

surv <- Survival(FineGrayMod)

nom.sur1<- nomogram(FineGrayMod,
                    fun=list(function(x) 1-surv(12,x),
                             function(x) 1-surv(36,x),
                             function(x) 1-surv(60,x)),
                    funlabel=c("1-year other causes Prob.",
                               "3-year other causes Prob.",
                               "5-year other causes Prob."),
                    lp=F)

plot(nom.sur1, force.label = TRUE)

# Set the desired output file path
output_file <- "nomogram_plot_1.tiff"

# Set the resolution in dots per inch (dpi)
dpi <- 500

# Set the width and height of the plot in inches
width <- 15
height <- 10

# Start the TIFF device and specify the output file, resolution, and dimensions
tiff(output_file, width = width, height = height, units = "in", res = dpi)

# Plot the nomogram
plot(nom.sur1)

# End the TIFF device
dev.off()

# extract coefficients which can be deployed via RSHiny.
p = length(nom.sur1) - 4

coef1 = NULL

for(i in 1:p){
  temp = data.frame(variable = nom.sur1[[i]][[1]],
                    score = nom.sur1[[i]][[3]])
  coef1 = rbind(coef1, temp)
}

risk1 = NULL

for(i in (p+2):(p+4)){
  temp = data.frame(risk = nom.sur1[[i]][[2]],
                    score = nom.sur1[[i]][[1]])
  risk1 = rbind(risk1, temp)
}

write.csv(coef1, "coef1.csv")
write.csv(risk1, "risk1.csv")

#Code for COD labeled as death due to non-cancer causes of death. 
FineGrayMod_OtherCOD <- cph(Surv(Tstart,Tstop,status==2)~agecat + sex + 
              race + site +
              brain + liver + lung + 
              tstage + nstage,
            data=df.w,
            weight=weight.cens,
            subset=failcode==2,
            surv=T)


# nomogram for cause 2
surv <- Survival(FineGrayMod_OtherCOD)
nom.sur2<- nomogram(FineGrayMod_OtherCOD,
                    fun=list(function(x) 1-surv(12,x),
                             function(x) 1-surv(36,x),
                             function(x) 1-surv(60,x)),
                    funlabel=c("1-year other causes Prob.",
                               "3-year other causes Prob.",
                               "5-year other causes Prob."),
                    lp=F)

# extract coefficients
p = length(nom.sur2) - 4

coef2 = NULL

for(i in 1:p){
  temp = data.frame(variable = nom.sur2[[i]][[1]],
                    score = nom.sur2[[i]][[3]])
  coef2 = rbind(coef2, temp)
}

risk2 = NULL

for(i in (p+2):(p+4)){
  temp = data.frame(risk = nom.sur2[[i]][[2]],
                    score = nom.sur2[[i]][[1]])
  risk2 = rbind(risk2, temp)
}

write.csv(coef2, "coef2.csv")
write.csv(risk2, "risk2.csv")

save.image("nomo.RData")
