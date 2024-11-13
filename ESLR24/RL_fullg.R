knitr::opts_chunk$set(echo = TRUE)
source('functions.R')
library(tidyverse)
library(ggplot2)

# Setting our environment

# Session Setup
set.seed(6) 
group_size = 10 
repetition = 10 # number of repetition for a specific combination of parameters
horizon = 100 # total steps of one session


# Bandit Setup
num_options = 2 # the number of arms
mean_list = c(1, 1.5) 
sd_list = c(0.03, 1)



# RL parameter space
Q_initial = 1.25 # prior belief
alpha_list <- c(0.3, 0.5, 0.7)
beta_list <- c(2, 4, 6)
sigma_list <- seq(0, 1, 0.2)
theta_list <- seq(1, 7, 2)


social_learning_model_data <- 
  data.frame(
    
    beta = rep(beta_list, times = 1, each = length(alpha_list)*length(sigma_list)*length(theta_list)*repetition*horizon)
    , alpha = rep(alpha_list, times = length(beta_list), each = length(sigma_list)*length(theta_list)*repetition*horizon)
    , sigma = rep(sigma_list, times = length(beta_list)*length(alpha_list), each = length(theta_list)*repetition*horizon)               
    , theta = rep(theta_list, times = length(beta_list)*length(alpha_list)*length(sigma_list), each = repetition*horizon)               
    , r = rep(1:repetition, times = length(beta_list)*length(alpha_list)*length(sigma_list)*length(theta_list), each = horizon)
    , t = rep(1:horizon, times = length(beta_list)*length(alpha_list)*length(sigma_list)*length(theta_list), each = 1)
    , bestChoiceProb = rep(NA, length(beta_list)*length(alpha_list)*length(sigma_list)*length(theta_list)*repetition*horizon)
    , averageReward = rep(NA, length(beta_list)*length(alpha_list)*length(sigma_list)*length(theta_list)*repetition*horizon)
  )


for (alpha in alpha_list) {
  for (beta in beta_list) {
    for (sigma in sigma_list) {
      for (theta in theta_list) {
        for (r in 1:repetition) {
          
          # Setting individual learning parameters
          alpha_base <- convert_prob_to_base(alpha)
          this_alpha <- (alpha_base + rnorm(group_size, 0, 0.05)) %>% convert_base_to_prob
          this_beta <- beta + rnorm(group_size, 0, 0.05)
          
          # -------- Social learning parameters --------------
          sigma_base <- convert_prob_to_base(sigma)
          this_sigma <- (sigma_base + rnorm(group_size, 0, 0.05)) %>% convert_base_to_prob
          this_theta <- theta + rnorm(group_size, 0, 0.05)
          # --------------------------------------------------
          
          # The initial settings
          Q <- array(dim = c(num_options, horizon, group_size))
          Q[,1,] <- Q_initial # Q_initial = 1.25
          netChoiceProb <- array(dim = c(num_options, horizon, group_size))
          netChoiceProb[,1,] <- 1/num_options
          choices <- matrix(nrow=group_size, ncol=horizon)
          payoffs <- matrix(nrow=group_size, ncol=horizon)
          bestChoiceProb <- matrix(nrow=group_size, ncol=horizon)
          
          # -- A matrix tracking the social frequency information ----
          socialFrequency = matrix(nrow=num_options, ncol=horizon)
          socialFrequency[,] = 1e-1
          # ----------------------------------------------------------
          
          
          for (t in 1:horizon) {
            # each individual chooses one option 
            choices[,t] <- apply(netChoiceProb[,t,]
                                 , MARGIN = 2
                                 , function(netChoiceProbInstance){ 
                                   sample(1:num_options
                                          , 1
                                          , prob=netChoiceProbInstance
                                          , replace=FALSE
                                   ) 
                                 }
            )
            
            # each subject earns payoffs whose amount depends upon choices[,t]
            payoffs[,t] <- mapply(rnorm, 1, mean_list[choices[,t]], sd_list[choices[,t]])
            
            if (t < horizon) {
              # Rescorla-Wagner updating 
              Q <- rescorla_wagner_updating(t
                                            , numoptions
                                            , group_size
                                            , horizon
                                            , this_alpha
                                            , choices[,t]
                                            , payoffs[,t]
                                            , Q
              )
              
              # ------------ Update socialFrequency ----------
              for (k in 1:num_options) {
                if(length(which(names(table(choices[,t]))==k))>0) {
                  socialFrequency[k,t+1] <- 
                    socialFrequency[k,t+1] +
                    table(choices[,t])[which(names(table(choices[,t]))==k)][1]
                }
              }
              # ---------------------------------------------
              
              # Calculate each "exp( beta*Q_i )"
              Q_exp = ( Q[,t+1,] * rep(this_beta, each = num_options) ) %>% exp() 
              # Then, calculating the softmax probability
              # exp( beta*Q_k )/(exp( beta*Q_1 ) + exp( beta*Q_2 )) for each slot
              softmaxMatrix = Q_exp %>% 
                apply(MARGIN=1
                      , divideVector
                      , denominator = apply(Q_exp, MARGIN=2, sum)
                ) %>% 
                t()
              
              # ------ Calculating the social influence --------
              # The function "frequencyDependentCopy" is defined in 
              # "functions.R" file
              freqDepenMatrix <- 
                frequencyDependentCopy(socialFrequency[,t+1]
                                       , choices[,t]
                                       , this_theta
                                       , num_options
                )
              
              netMatrix <- 
                apply(softmaxMatrix, 1, multiplying, A=(1-this_sigma)) %>%
                t() + 
                apply(freqDepenMatrix, 1, multiplying, A=this_sigma) %>% 
                t()
              # ------------------------------------------------
              
              netChoiceProbAperm = aperm(netChoiceProb, c(1,3,2))
              dim(netChoiceProbAperm) = c(num_options*group_size, horizon)
              dim(netMatrix) = c(num_options*group_size, 1)
              netChoiceProbAperm[,t+1] = netMatrix
              dim(netChoiceProbAperm) = c(num_options, group_size, horizon)
              netChoiceProb = aperm(netChoiceProbAperm, c(1,3,2))
            }
          }
          # ==========================================================
          
          for(i in 1:group_size) {
            bestChoiceProb[i,] <- netChoiceProb[2,,i]
          }
          thisGroupPerformance = apply(bestChoiceProb, MARGIN = 2, mean)
          thisGroupPayoff = apply(payoffs, MARGIN = 2, mean)
          # Second, let's record this in the data set you defined 
          # We first find a position where new data should be recorded
          positionOfThisRun = which(social_learning_model_data$beta == beta & 
                                      social_learning_model_data$alpha == alpha & 
                                      social_learning_model_data$sigma == sigma & 
                                      social_learning_model_data$theta == theta & 
                                      social_learning_model_data$r == r) 
          # Then, record it:
          social_learning_model_data$bestChoiceProb[positionOfThisRun] = thisGroupPerformance
          social_learning_model_data$averageReward[positionOfThisRun] = thisGroupPayoff
          # END
          # ==========================================================
           }
        }
      }
     }
   }

## summerize data          
social_learning_model_data_summarised <- social_learning_model_data %>% 
  # the new dataset is copied from the original dataset
  dplyr::group_by(alpha, beta,theta, sigma, t) %>% 
  summarise( # summarising the data set with the group structure considered.
    # the following two become new variables
    bestChoiceProb_globalMean = mean(bestChoiceProb), 
    averageReward_globalMean = mean(averageReward)
  )     


subsetData <- social_learning_model_data_summarised %>%
  filter(sigma == 0 &  theta == 5)
          

###plot
ggplot(subsetData) +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  geom_line(aes(t, bestChoiceProb_globalMean, colour = as.factor(beta))) +
  scale_colour_viridis_d() +
  ylim(c(0,1))+
  labs(x = 'Steps'
       , y = '% Optimal action'
       , colour = 'Inverse\nTemperature'
       , title = 'Social learning model') +
  facet_grid(. ~ alpha) +
  theme_test() 

          
          
          