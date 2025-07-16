#ECON484 Week 10
#16 July 2025
#Kevin Laughren

#####
# Overview of adaptive sampling exercise
#####
#We want to complete several adaptive sampling algorithms over four rounds of four actions
#We want to compare the number of responses (regret) of each algorithm, as well as their posterior distributions

#####
# Libraries
#####
library(stats)

#####
# Set parameters
#####

#create a vector of length (actions) that maps to x-variables
#$(1,0)\rightarrow k=1$, $(1,1)\rightarrow k=2$, $(2,0)\rightarrow k=3$, $(2,1)\rightarrow k=4$
actionad <- c(1,1,2,2)
actionman <- c(0,1,0,1)

#True response rates
actioncount <- 4
theta <- c(0.0050,0.0055,0.0060,0.0065)
actioncolours <- c("darkred","darkgrey","black","darkblue")

#declare a vector of probabilities to plug into the Beta() density function for plotting
#should exceed the upper and lower bounds of theta
prob = seq(0.00, 0.015,length = 151)

#Sample batch parameters
batchcount <- 10
#batchsize is per action. batchsize*4*batchcount should be equal to population
batchsize <- 2500

#Population parameters
population <- 1000000
actionpopulation <- round(population/actioncount,0)

#Name and count of sampling methods
methodcount <- 3
methodname <- c("greedy","Thompson","RCT")

#pre-define seeds to loop through during sampling
seeds <- c(500)
popseed <- 5150


#####
# Data simulation
#####
#(i) we use simulated data for this exercise so that we can 'know' what optimal performance should look like
#(ii) we simulate a large dataset at the beginning and sample from it (rather than simulate data each round)
#(iii) by doing this, we can include a sampling step after each batch as a reminder of when you need to go get more data

#create empty vectors to fill in a loop using "_pop" to represent a population variable
ad_pop <- c()
man_pop <- c()
action_pop <- c()
y_pop <- c()

#loop to generate data
#seed to isolate a single population from which we sample
set.seed(popseed)

for(k in 1:actioncount){
  #append x-variables (known)
  ad_pop <- append(ad_pop,rep(actionad[k],actionpopulation))
  man_pop <- append(man_pop,rep(actionman[k],actionpopulation))
  action_pop <- append(action_pop,rep(k,actionpopulation))
  #generate y-variables (random)
  y_k <- ifelse(runif(actionpopulation)<=theta[k],1,0)
  y_pop <- append(y_pop,y_k)
}

#create population dataframe with an index for sampling
index_pop <- c(1:population)

popdata <- data.frame(ad_pop,
                      man_pop,
                      y_pop,
                      action_pop,
                      index_pop)

#turn into a table to check response rates versus parameters
poptable <- table(popdata$man_pop,
                  popdata$ad_pop,
                  popdata$y_pop,
                  dnn=list("man",
                           "ad",
                           "response"))
poptable
#population data success rates
poptable[5:8]/actionpopulation

#####
#Plotting function
#####
#Plot the beliefs after one round, start with one, then loop the rest after declaring axes
customplot <- function(prob,alphas,betas,method){
  plot(x= prob,
       y = dbeta(prob,
                 shape1 = alphas[1],
                 shape2 = betas[1]),
       col = actioncolours[1], 
       type = "l", 
       lwd = 2.0,
       xlab = "Probability of y=1",
       ylab = "Density",
       ylim = c(0,400),
       family = "serif",
       frame.plot=F,
       yaxt="n",
       xaxt="n",
       main=paste("Sample",
                  s,
                  "probability densities after",
                  "batch",
                  b,
                  toString(method))
       )
  #customize the y-axis
  axis(side = 2,
       at = c(0,200,400)
       )
  #customize the x-axis  
  axis(side = 1,
       at = c(0,0.005,0.01,0.015),
       labels = c(0,0.005,0.01,0.015)
       )
  #plot the remaining actions
  for(k in 2:actioncount){
    lines(x= prob,
          y = dbeta(prob,
                    shape1 = alphas[k],
                    shape2 = betas[k]),
          col = actioncolours[k], 
          type = "l", 
          lwd = 2.0)
  }
  #add a legend
  legend(x="topright",
         legend=c("action 1",
                  "action 2",
                  "action 3",
                  "action 4"),
         lty=c(1,1,1,1),
         lwd = 2.0,
         col=actioncolours)
}


#####
#Nested Loop through (seeds, batches, methods)
#####

#define sample vectors
ad <- c()
man <- c()
y <- c()
action <- c()

#vectors to track when data was sampled
sample <- c()
batch <- c()
method <- c()

#set a single seed while setting up, turn off after
seeds <- c(1)

for(s in 1:length(seeds)){
  #loop through seeds
  set.seed(seeds[s])
  #record a variable to separate samples by seed
  sample <- append(sample,rep(s,batchsize*batchcount*methodcount))
  
  for(b in 1:batchcount){
    #record a variable to separate batches
    batch <- append(batch,rep(b,batchsize))
    
    #first batch must be the SAME for all algorithms, do not sample once per algorithm
    if(b==1){
      for(k in 1:actioncount){
        #identify which population data are the matching action and sample the appropriate amount
        actionrows <- which(popdata$action_pop == k)
        actionsample <- sample(actionrows,
                               size=round(batchsize/actioncount,0),
                               replace=FALSE)
        #append the data once per each method
        for(m in 1:methodcount){
          method <- append(method,rep(methodname[m],round(batchsize/actioncount,0)))
          ad <- append(ad,popdata$ad_pop[actionsample])
          man <- append(man,popdata$man_pop[actionsample])
          y <- append(y,popdata$y_pop[actionsample])
          action <- c(action,popdata$action_pop[actionsample])
        }
      }
      
      #####
      # Plot theta distributions after batch 1
      #####
      
      #All samples are the same for batch 1, only need to calculate and plot once
      thompsontable <- table(action[(sample==s)&(method=="Thompson")],
                             y[(sample==s)&(method=="Thompson")],
                             dnn=list("action",
                                      "response"))
      #get count of successes and failures (+1) for Beta() distribution
      thompsonalphas <- 1+thompsontable[,2]
      thompsonbetas <- 1+thompsontable[,1]
      
      #run customplot
      customplot(prob=prob,
                 alphas=thompsonalphas,
                 betas=thompsonbetas,
                 method="All sampling methods")
      
    }
    
    #####
    # Start batches 2+
    #####
    
    #subsequent batches need a separate loop for each method
    else{
      
      #####
      #Greedy sampling
      #####
      greedytable <- table(action[(sample==s)&(method=="greedy")],
                             y[(sample==s)&(method=="greedy")],
                             dnn=list("action",
                                      "response"))
      #get count of successes and failures (+1) for Beta() distribution
      greedyalphas <- 1+greedytable[,2]
      greedybetas <- 1+greedytable[,1]
      
      #identify the action with the highest expected response rate for sampling
      greedybest <- which.max(greedyalphas/(greedyalphas+greedybetas)) 
      greedyrows <- which(popdata$action_pop == greedybest)
      
      greedysample <- sample(greedyrows,
                             size=batchsize,
                             replace=FALSE)
      
      ad <- append(ad,popdata$ad_pop[greedysample])
      man <- append(man,popdata$man_pop[greedysample])
      y <- append(y,popdata$y_pop[greedysample])
      action <- c(action,popdata$action_pop[greedysample])
      
      #append the sampling method
      method <- append(method,rep("greedy",batchsize))
      
      #Update greedy table
      greedytable <- table(action[(sample==s)&(method=="greedy")],
                           y[(sample==s)&(method=="greedy")],
                           dnn=list("action",
                                    "response"))
      #get count of successes and failures (+1) for Beta() distribution
      greedyalphas <- 1+greedytable[,2]
      greedybetas <- 1+greedytable[,1]
      
      #####
      #Plot theta distributions for greedy batch 2+
      #####
      #run customplot
      customplot(prob=prob,
                 alphas=greedyalphas,
                 betas=greedybetas,
                 method="greedy sampling")
      
      
      #####
      #Thompson sampling
      #####
      thompsontable <- table(action[(sample==s)&(method=="Thompson")],
                             y[(sample==s)&(method=="Thompson")],
                             dnn=list("action",
                                      "response"))
      
      #get count of successes and failures (+1) for Beta() distribution
      thompsonalphas <- 1+thompsontable[,2]
      thompsonbetas <- 1+thompsontable[,1]
      
      #Thompson sampling draws once from each distribution for each observation
      thetahat <- c()
      for(k in 1:actioncount){
        thetas <- rbeta(n=batchsize,
                        shape1=thompsonalphas[k],
                        shape2=thompsonbetas[k])
        thetahat <- append(thetahat, thetas)
      }
      #store the random draws in a matrix
      thompsondraws <- matrix(data=thetahat,
                              nrow=batchsize,
                              ncol=actioncount)
      
      #for each observation, identify the max draw as the action you want to sample
      thompsonaction <- c()
      for(n in 1:batchsize){
        thompsonaction <- append(thompsonaction,
                                 which.max(thompsondraws[n,]))

      }
      #count how many of each action Thompson wants to draw in the next batch
      thompsonactioncount <- table(thompsonaction)
      
      #sample the new data
      
      for(k in 1:actioncount){
        #identify which population data are the matching action and sample the appropriate amount
        actionrows <- which(popdata$action_pop == k)
        actionsample <- sample(actionrows,
                               #Note the size is what changes in Thompson sampling
                               size=thompsonactioncount[k],
                               replace=FALSE)

        ad <- append(ad,popdata$ad_pop[actionsample])
        man <- append(man,popdata$man_pop[actionsample])
        y <- append(y,popdata$y_pop[actionsample])
        action <- c(action,popdata$action_pop[actionsample])
      }
      #append the sampling method
      method <- append(method,rep("Thompson",batchsize))
      
      #####
      # Plot theta distributions for Thompson batch 2+
      #####
      
      #Update Thompson table
      thompsontable <- table(action[(sample==s)&(method=="Thompson")],
                             y[(sample==s)&(method=="Thompson")],
                             dnn=list("action",
                                      "response"))
      #get count of successes and failures (+1) for Beta() distribution
      thompsonalphas <- 1+thompsontable[,2]
      thompsonbetas <- 1+thompsontable[,1]
      
      #run customplot
      customplot(prob=prob,
                 alphas=thompsonalphas,
                 betas=thompsonbetas,
                 method="Thompson sampling")
      
      
      #####
      #RCT sampling
      #####
      
      for(k in 1:actioncount){
        #identify which population data are the matching action and sample the appropriate amount
        actionrows <- which(popdata$action_pop == k)
        actionsample <- sample(actionrows,
                               #Equal size for all actions
                               size=round(batchsize/actioncount,0),
                               replace=FALSE)
        ad <- append(ad,popdata$ad_pop[actionsample])
        man <- append(man,popdata$man_pop[actionsample])
        y <- append(y,popdata$y_pop[actionsample])
        action <- c(action,popdata$action_pop[actionsample])
      }
      #append the method
      method <- append(method,rep("RCT",batchsize))
      
      #Update RCT table
      RCTtable <- table(action[(sample==s)&(method=="RCT")],
                        y[(sample==s)&(method=="RCT")],
                        dnn=list("action",
                                 "response"))
      #get count of successes and failures (+1) for Beta() distribution
      RCTalphas <- 1+RCTtable[,2]
      RCTbetas <- 1+RCTtable[,1]
      
      #####
      # Plot theta distributions for RCT batch 2+
      #####
      
      #run customplot
      customplot(prob=prob,
                 alphas=RCTalphas,
                 betas=RCTbetas,
                 method="RCT sampling")

    }
    
  }
  print(paste("seed",s,"tables and rates"))
  print("greedytable:")
  print(greedytable)
  print(greedytable[,2]/(greedytable[,1]+greedytable[,2]))
  print("thompsontable:")
  print(thompsontable)
  print(thompsontable[,2]/(thompsontable[,1]+thompsontable[,2]))
  print("RCTtable:")
  print(RCTtable)
  print(RCTtable[,2]/(RCTtable[,1]+RCTtable[,2]))
  
  print(paste("seed",s,"cumulative y=1 responses"))
  print("greedy1s")
  print(sum(greedytable[,2]))
  print("thompson1s")
  print(sum(thompsontable[,2]))
  print("RCT1s")
  print(sum(RCTtable[,2]))
}

#Look at cumulative outputs across multiple seeds
allseedstable <- table(action,y,method)





