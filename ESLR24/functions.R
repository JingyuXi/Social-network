
## Functions

# This function is simply dividing two vectors
divideVector = function(numerator, denominator) {
  return( numerator/denominator )
}
# This function is simply multiplying two things
multiplying = function(A, B) {
  return( A*B )
}

# Rescorla-Wagner updating function
rescorla_wagner_updating = function (t
                                     , numoptions
                                     , group_size
                                     , horizon
                                     , alpha
                                     , choices
                                     , payoffs
                                     , Q
                                     ) {
  # updatingPositions tracks which option was chosen by each individual at this trial
  updatingPositions = (choices + num_options*(1:group_size-1))
  
  # Copying the previous Q value to the current Q value
  Q[ , t+1, ] <- Q[ , t, ]
  
  # Q is a three-dimensional matrix, which is a bit tricky to handle in its current form.
  # So, let's reduce its dimension
  QQ <- aperm(Q, c(1,3,2))
  dim(QQ) <- c(num_options*group_size, horizon)
  
  # Here have I defined a new matrix, `QQ`, which is a reduced version of `Q`. 
  # Although its dimension was reduced to 2-D, it contains as same amount of
  # information as the original `Q`.
  
  # Then, Q values at the position `QQ[updatingPositions,t+1]` will be updated 
  # via Rescorla-Wagner rule:
  QQ[updatingPositions,t+1] <- 
    QQ[updatingPositions,t] + 
    alpha * (payoffs - QQ[updatingPositions,t])
  
  # Finally, we have to translate the updated `QQ` back to `Q` value
  dim(QQ) <- c(num_options, group_size, horizon)
  return( aperm(QQ, c(1,3,2)) )
}





multiplyBeta = function(Q, beta) {
    return( beta*Q )
}
expCeiling = function(x) {
    if(length(which(is.infinite(exp(x))))==0) {
        return( exp(x) )
    }else{
        result = exp(x)
        result[which(is.infinite(exp(x)))] = 8.988466e+307
        return( result ) # 2^1023
    }
}
powerTheta = function(f, theta) {
    return( f^theta )
}
frequencyDependentCopy_org = function(F, lastChoice, theta, numOptions) {
    totalN=length(lastChoice)
    f = rep(F, totalN)
    f[lastChoice+((1:totalN)-1)*numOptions] = f[lastChoice+((1:totalN)-1)*numOptions] - 1
    f_matrix = matrix(f, ncol=totalN)
    ftheta = apply(f_matrix, 1, powerTheta, theta=theta) %>% t()
    denom = apply(ftheta,2,sum)
    if(denom == 0) denom = 1
    return ( apply(ftheta, 1, divideVector, denominator=denom) %>% t() )
}
frequencyDependentCopy = function(F, lastChoice, theta, numOptions) {
    totalN=length(lastChoice)
    f = rep(F, totalN)
    f[lastChoice+((1:totalN)-1)*numOptions] = f[lastChoice+((1:totalN)-1)*numOptions] - 1
    f_matrix = matrix(f, ncol=totalN)
    ftheta = f_matrix ^ rep(theta, each= numOptions)
    denom = apply(ftheta,2,sum)
    if(length(which(denom == 0)>0)) denom[which(denom==0)] = 1
    return ( apply(ftheta, 1, divideVector, denominator=denom) %>% t() )
}
myTheme_legend = function() {
  theme(
    #legend.position = 'right',
    legend.title = element_text(size=16, family="Times" ,colour='black'),
    legend.text = element_text(size=15, family="Times" ,colour='black'),
    strip.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
    strip.text = element_text(size=16, family="Times" ), #"Times"
    #axis.line = element_line(colour = "#000000", size = 1, linetype = "solid", lineend = "round"),
    axis.text.x = element_text(size=15, family="Times" ,colour='black'),
    axis.text.y = element_text(size=15, family="Times" ,colour='black'),
    axis.title.x=element_text(size=16, family="Times" ),
    axis.title.y=element_text(size=16, family="Times" ),
    panel.background = element_rect(fill = "white", colour = NA),
    #panel.grid = element_blank(),
    panel.grid = element_line(colour = "black", size = 0.8),
    plot.title = element_text(size=15, family="Times" ),
    plot.background = element_rect(colour = "white")
  )
}
myTheme_gillsansMT = function() {
  theme(
    #legend.position = 'none',
    strip.background = element_rect(fill = NA),
    #panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
    strip.text = element_text(size=16, family="Gill Sans MT" ), #"Gill Sans MT"
    axis.line = element_line(colour = "#000000", size = 0.5, linetype = "solid", lineend = "round"),
    legend.text = element_text(size=15, family="Gill Sans MT" ,colour='black'),
    legend.title = element_text(size=16, family="Gill Sans MT" ,colour='black'),
    axis.text.x = element_text(size=15, family="Gill Sans MT" ,colour='black'),
    axis.text.y = element_text(size=15, family="Gill Sans MT" ,colour='black'),
    axis.title.x=element_text(size=16, family="Gill Sans MT" ),
    axis.title.y=element_text(size=16, family="Gill Sans MT" ),
    plot.title = element_text(size=15, family="Gill Sans MT", hjust=0.5),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    plot.background = element_rect(colour = "white"),
    panel.spacing = unit(0.2, "lines")
  )
}
myTheme_gillsans = function() {
  theme(
    #legend.position = 'none',
    strip.background = element_rect(fill = NA),
    #panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
    strip.text = element_text(size=16, family="Gill Sans" ), #"Gill Sans"
    axis.line = element_line(colour = "#000000", size = 0.5, linetype = "solid", lineend = "round"),
    legend.text = element_text(size=15, family="Gill Sans" ,colour='black'),
    legend.title = element_text(size=16, family="Gill Sans" ,colour='black'),
    axis.text.x = element_text(size=15, family="Gill Sans" ,colour='black'),
    axis.text.y = element_text(size=15, family="Gill Sans" ,colour='black'),
    axis.title.x=element_text(size=16, family="Gill Sans" ),
    axis.title.y=element_text(size=16, family="Gill Sans" ),
    plot.title = element_text(size=15, family="Gill Sans", hjust=0.5),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    plot.background = element_rect(colour = "white"),
    panel.spacing = unit(0.2, "lines")
  )
}

myTheme_Times = function() {
  theme(
    #legend.position = 'none',
    strip.background = element_rect(fill = NA),
    #panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
    strip.text = element_text(size=16, family="Times" ), #"Times"
    axis.line = element_line(colour = "#000000", size = 0.5, linetype = "solid", lineend = "round"),
    legend.text = element_text(size=15, family="Times" ,colour='black'),
    legend.title = element_text(size=16, family="Times" ,colour='black'),
    axis.text.x = element_text(size=15, family="Times" ,colour='black'),
    axis.text.y = element_text(size=15, family="Times" ,colour='black'),
    axis.title.x=element_text(size=16, family="Times" ),
    axis.title.y=element_text(size=16, family="Times" ),
    plot.title = element_text(size=15, family="Times", hjust=0.5),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    plot.background = element_rect(colour = "white"),
    panel.spacing = unit(0.2, "lines")
  )
}

myTheme_Helvetica = function() {
  theme(
    #legend.position = 'none',
    strip.background = element_rect(fill = NA),
    #panel.border = element_rect(fill = NA, colour = "black", size = 0.8),
    strip.text = element_text(size=13, family="Helvetica" ), #"Helvetica"
    axis.line = element_line(colour = "#000000", size = 0.2, linetype = "solid", lineend = "round"),
    legend.text = element_text(size=13, family="Helvetica" ,colour='black'),
    legend.title = element_text(size=13, family="Helvetica" ,colour='black'),
    axis.text.x = element_text(size=13, family="Helvetica" ,colour='black'),
    axis.text.y = element_text(size=13, family="Helvetica" ,colour='black'),
    axis.title.x=element_text(size=13, family="Helvetica" ),
    axis.title.y=element_text(size=13, family="Helvetica" ),
    plot.title = element_text(size=13, family="Helvetica", hjust=0.5),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    plot.background = element_rect(colour = "white"),
    panel.spacing = unit(0.2, "lines")
  )
}

## Converting function
convert_prob_to_base = function (alpha) {
    alpha[which(alpha < 1e-5)] <- 1e-5
    alpha[which(alpha > 1 - 1e-5)] <- 1 - 1e-5
    return (log(-alpha/(alpha-1)))
}

convert_base_to_prob = function (alphaRaw) {
    alphaRaw[which(alphaRaw > 12)] <- 12
    return( 1 / (1 + exp(-alphaRaw)))
}



SN_frequencyDependent = function(F, lastChoice, theta, numOptions) {
  ftheta = socialFrequency[,t+1,] ^ rep(theta, each= numOptions)
  denom = apply(ftheta,2,sum)
  if(length(which(denom == 0)>0)) denom[which(denom==0)] = 1
  return ( apply(ftheta, 1, divideVector, denominator=denom) %>% t() )
}




makeSymm<- function(SN) {
  SN[upper.tri(SN)] <- t(SN)[upper.tri(SN)]
  return(SN)
}


