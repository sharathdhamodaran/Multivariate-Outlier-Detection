#-------------------------------------------------------------------------------
# Script Name: Multivariate Outlier Detection using Hotelling
# Author: Sharath Kumar Dhamodaran
#-------------------------------------------------------------------------------

# setup
# load the libraries & functions
source("scripts/libraries.R")
source("scripts/functions/get_max_pct.R")
source("scripts/functions/mod_p1.R")
source("scripts/functions/mod_p2.R")
source("scripts/functions/boxplot_outliers.R")
source("scripts/functions/contribution_plot.R")
source("scripts/functions/get_params_to_remove.R")
source("scripts/functions/unknown_dec_rule.R")
source("scripts/functions/t2_plot.R")

# Data extraction --------------------------------------------------------------

#  make a file path by specifying the location
filesLocation <- here("data")

# look for all the .csv files and makes a list of them called "files"
files <- dir(filesLocation, pattern = "*.csv")
files

# bind all the .csv's into one dataframe
partData <- files %>%
  # read in all the files, appending the path before the filename
  map(~ read_csv(file.path(filesLocation, .))) %>%
  # reduce with rbind into one dataframe
  reduce(rbind) %>%
  data.frame()

# peak into the windspeed data
head(partData)

# Data Manipulation ------------------------------------------------------------

setDT(partData)

# remove numeric coulmns based on different conditions in the function
removeCols <- get_params_to_remove(subset(partData,
                                          select = -c(CancerPatients)
                                          )
                                   )

# remove the columns returned from the function from the data
partData <- partData[, -which(names(partData) %in% removeCols),
                     with = FALSE
                     ]

# remove rows with missing values
partData <- partData[complete.cases(partData[, sapply(partData, is.numeric),
                                             with = FALSE
                                             ])]

#  Get all the params from the data
params <- colnames(partData)
params <- params[!params == "CancerPatients"]

# convert the data back to data frame for further analysis
partData <- data.frame(partData)

# Modeling & Visualization------------------------------------------------------

# Split the data into train and test for accomplishing phase I & II of MOD
# set a random number to obtain reproducible results
set.seed(123)
# split the data in the ratio mentioned in SplitRatio.
sample <- sample.split(partData, SplitRatio = 0.95)
# create a training dataset with rows which are marked as TRUE
train <- subset(partData, sample == TRUE)
test <- subset(partData, sample == FALSE)

# Phase I - Multivariate Outlier Detection
p1 <- mod_p1(in.data = train, in.vars = params, in.id_var = "CancerPatients")

# Phase II - Multivariate Outlier Detection
p2 <- mod_p2(p1, in.score_data = test)

# Visualizing OOC contributions plots
plots <- contribution_plot(p2)
plots$plot1

# visualizing Hotelling T2 plot
t2_plot(p2)

# clear all workspace objects and retain only what we want
# rm(list = setdiff(ls(), c("p1", "p2")))

