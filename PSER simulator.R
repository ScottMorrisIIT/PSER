# PSER simulator
# Plots the results of the PSER stopping for a range of hypo/hyper parameter values
# Created by (obscured for blind review), June 2018

# The simulaiton requires a sequence of SE estiamtes from responses to a CAT for multiple individuals.
# This can be obtained in several ways:
# 1. Use the results of a CAT that select items using the PSER algorithm (Choi et al, 2010). Read in the sequence of PSER estiamtes.
# 2. Use the results of a CAT that applyies any item selection method. Input the sequence of SE estimates.
# 3. Compute the sequence of PSER values implied by the Fisher Information for set of know theta values.

# Setup -------------------------------------------------------------------

# Required libraries
library(reshape2)

# Path and name of data files
CAT_data_filename <- "CAT_data.csv"  # needed if using actual/simulated CAT input
item_parameter_filename <- "item_pararameters.csv"  # needed if simulating CAT from item parameters
trait_distribution_filename <- "score_data.csv"  # needed if simulaitng CAT from item parameters with empirical trait distribution


# Functions ---------------------------------------------------------------
# Functions required depends on the method used. All methods use the pser() function.
# CAT data input: get_CAT_sequence ()
# Item parameter input: item_informatio_GRM (), simulate_CAT_sequence()

# Function to apply PSER algorighm for a specified range of hypo and hyper levels
pser <- function (dat, hypo, hyper_inc, SE_cutoff, max_cat_length) {
  
  IDList <- unique(dat$ID)
  n_obs <- length(IDList)
  #SE_final <- CAT_length <- array(0,dim=c(length(hypo),length(hyper_inc),n_obs))
  result <- data.frame (ID = integer(), hypo = double(), hyper_inc = double(), SE_final = double(), CAT_length = integer())
  
  for (h1 in 1:length(hypo))
  {
    for (h2 in 1:length(hyper_inc))
    {
      hyper <- hypo[h1] + hyper_inc
      
      for (t in 1:n_obs)
      {
        for (round in 2:max_cat_length)
        {
          CAT_SE <- dat$SE[dat$ID == IDList[t] & dat$Position == round]
          SE_change <- dat$SE_change[dat$ID == IDList[t] & dat$Position == round]
          #SE_final[h1,h2,t] <- CAT_SE 
          #CAT_length[h1,h2,t] <- round
          if(SE_change < hypo[h1] | (CAT_SE < SE_cutoff & SE_change < hyper[h2]))
          {
            break
          } 
        }
        tmp <- data.frame(ID = IDList[t], hypo = hypo[h1], hyper_inc = hyper_inc[h2], SE_final = CAT_SE, CAT_length = round)
        result <- rbind(result, tmp)
        
      }
    }
  }
  return(result)
} # end of pser function

# Function to read expected or actual SE from sequence of CAT responses
get_CAT_sequence <- function (df, vars, ComputeChange = TRUE) {
  
  # This function assigns values from the user's data to the dataframe CATData, and coimputes SE_change.
  newdf <- data.frame(ID = df[variable_labels[1]],Position = NA, SE = NA, SE_change = NA) 
  newdf$Position <-  df[,variable_labels[2]]
  newdf$SE <- df[,variable_labels[3]]
  if (ComputeChange) {
    newdf$SE_change <- NA
  } else {
    newdf$SE_change <- df[variable_labels[4]]
  }
  
  IDList <- unique(newdf$ID)
  first_item <- min(newdf$Position) # index for first round of the CAT
  
  # Compute SE_change
  if (ComputeChange) {
    newdf$SE_change[newdf$Position == first_item] <- 1.0 - newdf$SE[newdf$Position == first_item]  # for the first item SE_change = 1 - SE
    for (t in 1:length(unique(newdf$ID)))
    {
      for (round in newdf$Position[newdf$ID == IDList[t] & newdf$Position != first_item])
      {
        newdf$SE_change[newdf$ID == IDList[t] & newdf$Position == round] <- newdf$SE[newdf$ID == IDList[t] & newdf$Position == (round-1)] - newdf$SE[newdf$ID == IDList[t] & newdf$Position == round]
      }
    }
  }
  return (newdf)
}
# end of get_CAT_sqeuence function

# Function to compute item information function for the GRM
# Alternatively, you can skip this an specify a (#items x #theta ) item information matrix 
item_information_GRM <- function (prm, nTheta, firstTheta, lastTheta) {
  
  max_nCat <- max(prm$nCat) # maximum number of category parameters
  nItm <- nrow(prm)   # number of items
  ccol <- 3  # column of prm containing the first category location (c) parameter
  
  brf <- array(0,dim=c(nItm,nTheta,max_nCat-1 ))
  prf <- array(0,dim=c(nItm,nTheta,max_nCat ))
  info <- array(0,dim=c(nItm,nTheta ))
  
  for (i in 1:nItm)
  {
    for (t in 1:nTheta)
    {
      theta <- firstTheta + (lastTheta - firstTheta)/(nTheta - 1) * (t-1)
      
      for (cat in 1:(prm$nCat[i]-1))
      {
        brf[i,t,cat] <- 1/(1 + exp(-prm$a[i]*theta - prm[i,(ccol+cat-1)]))
      }
      prf[i,t,1] <- 1 - brf[i,t,1]
      prf[i,t,prm$nCat[i]] <- brf[i,t,(prm$nCat[i]-1)]
      info [i,t] <- (brf[i,t,1]*(1-brf[i,t,1]))^2/prf[i,t,1] + (brf[i,t,(prm$nCat[1]-1)]*(1-brf[i,t,(prm$nCat[i]-1)]))^2/prf[i,t,prm$nCat[i]]
      for (cat in 2:(prm$nCat[i]-1))
      {
        prf[i,t,cat] <- brf[i,t,cat-1] - brf[i,t,cat]
        info [i,t] <- info[i,t] + (brf[i,t,cat-1]*(1-brf[i,t,cat-1])-brf[i,t,cat]*(1-brf[i,t,cat]))^2/prf[i,t,cat]
      }
      info[i,t] <- info[i,t]*prm$a[i]^2 
    }
  }
  return(info)
}  # end of item_information_GRM function

#################################################################.
# Function to determine optimal squence of items for each theta and CAT information at each round of the CAT.
# To simulate Bayesian trait estimation, we add the inverse of the prior trait variance to the CAT information (Segall, 2010).
# For standard normal prior, the Bayesian adjustment is 1.0.
simulate_CAT_sequence <- function (info, nTheta, firstTheta, lastTheta) {
  
  nItm <- nrow(info)
  Bayesian_adjustment <- 1.0
  t <- c(1:nTheta)
  theta <- firstTheta + (lastTheta - firstTheta)/(nTheta - 1) * (t-1)
  max_rounds <- nItm
  max_info <- best_item <- CAT_info <- CAT_SE <- SE_change <- array (0,dim=c(nTheta,max_rounds))
  used <- array(0, dim=c(n_theta,nItm))
  
  for (t in 1:nTheta)
  {
    for (round in 1:max_rounds)
    {
      
      for (i in 1:nItm)
      {
        if (info[i,t] > max_info[t,round] & used[t,i] == 0)
        {
          max_info[t,round] <- info[i,t]
          best_item[t,round] <- i
        }
      }
      used[t,best_item[t,round]] <- 1
      if (round == 1)
      {
        #  Bayesian_adjustment: Add 1 to CAT information to reflect influence of standarad normal prior (Brown & Croudace, 2015, p. 313)
        CAT_info[t,round] <- max_info[t,round] + Bayesian_adjustment 
        CAT_SE[t,round] <- 1/sqrt(CAT_info[t,round])
        SE_change[t,round] <- NA
      } else
      {
        CAT_info[t,round] <- CAT_info[t,round-1]+max_info[t,round]
        CAT_SE[t,round] <- 1/sqrt(CAT_info[t,round])
        SE_change[t,round] <- max((CAT_SE[t,(round-1)] - CAT_SE[t,round]),0)
      }
    }
  }
  
  # Create dataframe of results
  tmp1 <- data.frame (ID = seq(1,nTheta), CAT_SE)
  tmp2 <- melt(tmp1, id = "ID")
  tmp2$Position <- rep(seq(1,ncol(CAT_SE)),each=nTheta)
  names(tmp2)[names(tmp2) == 'value'] <- 'SE'
  tmp1 <- data.frame (ID = seq(1,nTheta), SE_change)
  tmp3 <- melt(tmp1, id = "ID")
  tmp3$Position <- rep(seq(1,ncol(SE_change)),each=nTheta)
  names(tmp3)[names(tmp3) == 'value'] <- 'SE_change'
  result <- data.frame(ID = tmp2$ID, Position = tmp2$Position, SE = tmp2$SE, SE_change = tmp3$SE_change)
  result <- result[order(result$ID,result$Position),]
  #View(CATData)
  
  return(result)
  
} # end of simulate_CAT_sequence function


# Main Program ------------------------------------------------------------

# Specify trait range and number of quadrature points 
n_theta <- 41
min_theta <- -4
max_theta <- 4
t <- c(1:n_theta)
theta <- round(min_theta + (max_theta - min_theta)/(n_theta - 1) * (t-1),1)
remove(t)

# Step 1: Create a sequence of CAT results with SE change (expected or actual) after each item.
# Use only one of the following approaches: CAT Input or Item Parameter Input

# Method 1: CAT input ---------------------------------------------------------------
# Generating PSER from a sequence of CAT responses with actural or expected SE. 

# load data file
my_CAT_data <- read.csv(CAT_data_filename, header = TRUE)

# Create a dataframe'CATData' with variables: ID, Position (index for item adminstration order), SE & SE_change. 
# If SE_change does not exist, the function get_CAT_sequence can be used to create the dataframe.

# to manuallyl create CATData
# CATData <- data.frame(ID = my_CAT_data$ID)
# CATData$Position <- my_CAT_data$Position
# CATData$SE <- my_CAT_data$SE
# CATData$SE_change <- my_CAT_data$SE_change

# Use get_CAT_data function to compute SE_change and create CATData
variable_labels <- c("ID","Position","SE")  # specify variable names in my_CAT_data
CATData <- get_CAT_sequence(my_CAT_data, variable_labels, ComputeChange = TRUE)
#View(CATData)

# Compute data characteristics used in simulation
IDList <- unique(CATData$ID)
n_obs <- length(IDList)
max_rounds <- max(CATData$Position)
wt <- array(1, dim = n_obs) # assign unit weights when computing average results

# Method 2: Item Parameter Input ----------------------------------------------------
# Generating PSER from item parameters with no item response data

# Read item parameters and assign to datafram 'item_parameters'
# Note: the following assumes item parameters are from the GRM, and
# in the logit = a * THETA + c form. 

# read item parameters from file
my_parameters <- read.csv(item_parameter_filename, header = TRUE)

# identify relevant columns of my_parameters
nCat <- my_parameters$nCat  # vector indicating number of response cateories for each item
a <- my_parameters$a
cmatrix <- my_parameters[,c("c1","c2","c3","c4")]
colnames(cmatrix) <- paste0("c",1:(max(nCat)-1))
item_parameters <- data.frame(nCat,a,cmatrix)

# Compute data descriptives to use used later
n_items <- nrow(item_parameters)   # number of items
n_obs <- n_theta

# call function to compute item information
info <- item_information_GRM (item_parameters, n_theta, min_theta, max_theta)

# call function to generate the optimal item sequence and CAT information at each level of theta
CATData <- simulate_CAT_sequence (info, n_theta, min_theta, max_theta)
#View(CATData)

# Assign weights for computing RMSE and average length. Two methods available

# 1. Normal theta distribution. Use when computing performance from item parameters and trait distribution is unknown. 
wt <- dnorm(theta,0,1) 

# 2. Empirical theta distribution. Use when normality assumption not appropriate and trait estimates are available.
# Theta distribution will be computed from a data file with trait scores
scores <- read.csv(trait_distribution_filename, header = TRUE)$scores
wt <- hist(scores,breaks = c(theta,(max(theta)+.01)), include.lowest = TRUE, plot = FALSE)$density

# Step 2: Apply PSER algorithm --------------------------------------------

# Set range of hypo & hyper parameters to examine
hypo_levels <- seq(0,.05,.005)   
hyper_increment_levels <- c(0.005,.01, .015, 1.0)  # hyper = hypo + hyper_inc
#hyper_inc <- 1.0  # set hyper_inc = 1 to cause this parameter to be ignored
SE_threshold <- 0.3
max_rounds <- n_items
CAT_result <- pser(CATData, hypo_levels, hyper_increment_levels, SE_threshold, max_rounds)


# Step 3: Computer average results ----------------------------------------
# Compute weighted average SE and CAT length at each hypo/hyper combination

RMSE <- AveLength <- array(NA, dim = c(length(hypo_levels),length(hyper_increment_levels)))
for (h1 in 1:length(hypo_levels))
{
  for (h2 in 1:length(hyper_increment_levels))
  {
    tmp <- subset(CAT_result, (hypo == hypo_levels[h1] & hyper_inc == hyper_increment_levels[h2]))
    RMSE[h1,h2] <- sqrt(sum(tmp$SE_final^2 * wt)/sum(wt))
    AveLength[h1,h2] <- sum(tmp$CAT_length * wt)/sum(wt)
  }
}


# Step 4: Plot results ----------------------------------------------------
# Three types of precision/efficiency plots are available
# a) Average SE by CAT length at varying hypo levels at fixed hyper value
# b) Average SE by CAT length at each combination of hypo/hyper levels
# c) SE by CAT length at varying hypo levles at specific trait levels (available for item parameter input only) 

# Plot average precision/efficiency tradeoff for a single hyper value
h2 <- length(hyper_increment_levels)  # This uses the last value in the hyper list, which should be set at 1.0 
#h1 <- match(c(0,.005,.01,.015),hypo) # examine a subset of hypo levels
h1 <- 1:length(hypo_levels)
xmax <- max(AveLength[,h2])  # upper limit of x axis to plot
yrange <- c(.2,.7)
par(mar = c(5.1,4.1,4.1,4.5))
plot(AveLength[h1,h2],RMSE[h1,h2],type = 'o',ylim = yrange,xlim = c(0,xmax), pch = seq(1:length(h1)), lwd = 1,
     ylab = "Ave SE",xlab = "Ave Number of Items",
     main = "Average Precision/Efficiency of Alternate Hypo Levels")
#text(AveLength[,h2],RMSE[,h2],labels = hypo, pos = 3, cex = 0.7)  # add hypo level labels
legend(x=(xmax * 1.05), y = (yrange[2]),legend = hypo_levels[h1],pch = seq(1,length(h1)), title = "Hypo",
       pt.cex = 1, cex = min(1,(11/length(h1))), xpd = TRUE)
if (hyper_increment_levels[h2] >= 1.0) 
  mtext("Note: Hyper rule not applied", side = 1, outer = TRUE, line = -1, adj = 0, cex = 0.8)


# Plot average precision/efficiency tradeoff at multiple hypo and hyper levels
h1 <- match(c(0,.005,.01,.015),hypo_levels)
#h1 <- 1:length(hypo_levels)
h2 <- 1:length(hyper_increment_levels)
xmax <- max(AveLength)
#xmax <- 20
yrange <- c(min(RMSE),max(RMSE)+.05)
par(mar = c(5.1, 4.1, 4.1, 5.1))
plot(AveLength[h1[1],],RMSE[h1[1],],type='o', ylim = yrange, xlim = c(0,xmax), pch = 1, col = seq(1:length(h2)),
     ylab = "Ave SE", xlab = "Ave Number of Items", 
     main = "Average Precision/Efficiency of Alaternate Hypo/Hyper Levels")
for (h in 2:length(h1))
{
  lines(AveLength[h1[h],],RMSE[h1[h],], lwd = 1)
  points(AveLength[h1[h],],RMSE[h1[h],],pch = h, col = seq(1:length(h2)) )
}
legend("top", horiz = TRUE,legend = hyper_increment_levels, pch = 1, col = seq(1,length(h2)), title = "Hyper Increment", xpd = TRUE, inset = c(0,0))
legend("right",legend = hypo_levels[h1],pch = seq(1,length(h1)), pt.cex = 1, title = "Hypo", cex = min(1,(10/length(h1))), inset = c(-0.17,0), xpd = TRUE)
mtext("Note: Hyper = Hypo + Hyper_Increment", side = 1, outer = TRUE, line = -1, adj = 0, cex = 0.8)

# Plot precision/efficienty tradeoff at several theta levels (at single hyper level)
# NOTE: This plot is not available when useing actual CAT responses
h2 <- length(hyper_increment_levels) # this sets hyper at its last value, which should be set at 1.0
theta_set <- c(-1.2, -1.0, -0.8) # Specify desired theta values. 
tindex <- match(theta_set,theta)
tmp <- subset(CAT_result, hyper_inc == hyper_increment_levels[h2] & ID %in% tindex)
yrange <- c(min(tmp$SE_final),max(tmp$SE_final)+.05)
xmax <- max(tmp$CAT_length) + 1
#xmax = 15   # manually set xmax to zoom in on X axis
par(mar = c(4.8,4.1,4.1,5.1))
plot(tmp$CAT_length[tmp$ID == tindex[1]],tmp$SE_final[tmp$ID == tindex[1]],type = 'o',lty=1, pch = seq(1:length(hypo_levels)), col = 1,
     ylim = yrange, xlim = c(0,xmax),
     ylab = "SE",xlab = "Number of Items",
     main = "Precision/Efficiency of Alternate Hypo Levels by Theta")
for (i in 2:length(tindex))
{
  lines(tmp$CAT_length[tmp$ID == tindex[i]],tmp$SE_final[tmp$ID == tindex[i]],lty = i, lwd = 2, col = i )
  points(tmp$CAT_length[tmp$ID == tindex[i]],tmp$SE_final[tmp$ID == tindex[i]],pch = seq(1:length(hypo_levels)), col = i )
}
legend(x = (xmax + 0.5), y = (yrange[2] + 0.05), legend = theta[tindex],lty = seq(1,length(tindex)), lwd = 2, col = seq(1,length(tindex)), title = "Theta", xpd = TRUE)
legend(x = (xmax + 0.5), y = (yrange[1]+(yrange[2]-yrange[1])/1.5), legend = hypo_levels,pch = seq(1,length(hypo_levels)), pt.cex = 1, title = "Hypo", cex = min(1,(11/length(hypo_levels))), xpd = TRUE)
if (hyper_increment_levels[h2] >= 1.0) 
  mtext("Note: Hyper rule not applied", side = 1, outer = TRUE, line = -1, adj = 0, cex = 0.8)
