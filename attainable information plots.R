# Attainable Information Graphs
# Graphs  bank information and standard error attainable by Computer Adaptive Test,
# based on IRT item parameters from the Graded Response Model. Requires input files with item parameters,
# and optionally trait estimates for a representative sample of the target population.
# Created by (obscured for blind review), May 2019

# Setup -------------------------------------------------------------------

# Path and name of data files
item_parameter_filename <- "item_pararameters.csv"  # needed if simulating CAT from item parameters
trait_distribution_filename <- "score_data.csv"  # needed if simulaitng CAT from item parameters with empirical trait distribution

# Functions ---------------------------------------------------------------
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
  tmp0 <- data.frame (ID = seq(1,nTheta), CAT_info)
  tmp1 <- melt(tmp0, id = "ID")
  tmp1$Position <- rep(seq(1,ncol(CAT_SE)),each=nTheta)
  names(tmp1)[names(tmp1) == 'value'] <- 'Bank_Info'
  tmp0 <- data.frame (ID = seq(1,nTheta), CAT_SE)
  tmp2 <- melt(tmp0, id = "ID")
  tmp2$Position <- rep(seq(1,ncol(CAT_SE)),each=nTheta)
  names(tmp2)[names(tmp2) == 'value'] <- 'SE'
  tmp0 <- data.frame (ID = seq(1,nTheta), SE_change)
  tmp3 <- melt(tmp0, id = "ID")
  tmp3$Position <- rep(seq(1,ncol(SE_change)),each=nTheta)
  names(tmp3)[names(tmp3) == 'value'] <- 'SE_change'
  result <- data.frame(ID = tmp2$ID, Position = tmp2$Position, Bank_Info = tmp1$Bank_Info, SE = tmp2$SE, SE_change = tmp3$SE_change)
  result <- result[order(result$ID,result$Position),]
  #View(CATData)
  
  return(result)
  
} # end of simulate_CAT_sequence function



# Main Program ------------------------------------------------------------

# Specify target precisions - used only to generate reference line
SE_cutoff <- .3   # Specify threshold for SE stopping rule
info_cutoff <- 1/(SE_cutoff^2)  # Test information corresponding to SE cutoff

# Specify trait range and number of quadrature points 
n_theta <- 41
min_theta <- -4
max_theta <- 4
t <- c(1:n_theta)
theta <- round(min_theta + (max_theta - min_theta)/(n_theta - 1) * (t-1),1)
remove(t)

# Read item parameters and assign to datafram 'item_parameters'
# Note: the following assumes item parameters are from the GRM, and
# in the logit = a * THETA + c form. 

# read item parameters from file
my_parameters <- read.csv(item_parameter_filename, header = TRUE)

# extract relevant columns of my_parameters
nCat <- my_parameters$nCat  # vector indicating number of response cateories for each item
a <- my_parameters$a
cmatrix <- my_parameters[,c("c1","c2","c3","c4")]
colnames(cmatrix) <- paste0("c",1:(max(nCat)-1))
item_parameters <- data.frame(nCat,a,cmatrix)

# Compute data descriptives to use used later
n_items <- nrow(item_parameters)   # number of items

# call function to compute item information
info <- item_information_GRM (item_parameters, n_theta, min_theta, max_theta)

# call function to generate the optimal item sequence and CAT information at each level of theta
CATData <- simulate_CAT_sequence (info, n_theta, min_theta, max_theta)
#View(CATData)

# compute full bank information
BankInfo <- CATData[CATData$Position == n_items,"Bank_Info"]
BankSE <- CATData[CATData$Position == n_items,"SE"]

# Load empirical theta distribution. 
# Theta distribution will be computed from a data file with trait scores
scores <- read.csv(trait_distribution_filename, header = TRUE)$scores
wt <- hist(scores,breaks = c(theta,(max(theta)+.01)), include.lowest = TRUE, plot = FALSE)$density


# Graph Attainable Information with Trait Distribution --------------------------
trait_label <- "Attainable Information"
par(lwd=2)
par(mar=c(5,4,4,4) + 0.3)
maxinfo <- max(BankInfo)
t <- c(1:n_theta)
theta <- min_theta + (max_theta - min_theta)/(n_theta - 1) * (t-1)
hist(scores, freq = FALSE, xlim = c(-4,4),density = 10,col='darkgreen',axes = FALSE, xlab = "", ylab = "",main="")
axis(side=4, at = pretty(range(wt)))
par(new = TRUE)
plot(theta,BankInfo,type='l',col='black',lty=1,lwd=3,ylim=c(0,maxinfo),xlab="Trait Level",ylab="Bank Information",main=trait_label)
axis(side=1, at = pretty(range(BankInfo)))
mtext("Trait Frequency",side=4,line=2)

# Graph Attainable Standard Error with Trait Distribution --------------------------
trait_label <- "Attainable Standard Error"
par(lwd=2)
par(mar=c(5,4,4,4) + 0.3)
maxSE <- max(BankSE)
t <- c(1:n_theta)
theta <- min_theta + (max_theta - min_theta)/(n_theta - 1) * (t-1)
hist(scores, freq = FALSE, xlim = c(-4,4),density = 10,col='darkgreen',axes = FALSE, xlab = "", ylab = "",main="")
axis(side=4, at = pretty(range(wt)))
par(new = TRUE)
plot(theta,BankSE,type='l',col='black',lty=1,lwd=3,ylim=c(0,maxSE),xlab="Trait Level",ylab="Full Bank Expected SE",main=trait_label)
axis(side=1, at = pretty(range(BankSE)))
mtext("Trait Frequency",side=4,line=2)

