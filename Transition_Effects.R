EffectsTransition<-function(P_effects, b)
{  
  
  if(P_effects=="test")
  {   
    # here we assume that the probability that a site occupied is being replaced, colonized or extinct by a species is independent of the the other species dynamics
    # POCHERS EXTINCTION independent of the presence of prey
    epsilonAB <- exp(b[13]) / (1+exp(b[13])+exp(b[14])+exp(b[15]))
    trans13 <- exp(b[5]) / (1+exp(b[4])+exp(b[5])+exp(b[6])) # gammaB-gammaAB
    gammaAB <- exp(b[6]) / (1+exp(b[4])+exp(b[5])+exp(b[6])) 
    trans42 <- exp(b[14]) / (1+exp(b[13])+exp(b[14])+exp(b[15]))  # epsilonB-epsilonAB
    trans43 <- exp(b[15]) / (1+exp(b[13])+exp(b[14])+exp(b[15])) 
    trans12 <- exp(b[4]) / (1+exp(b[4])+exp(b[5])+exp(b[6]))
    nuB <- trans42 
    
    omegaAB <- exp(b[8]) / (1+exp(b[15])+exp(b[8])+exp(b[9]))
    etaB <- exp(b[9]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
    nuA <- exp(b[7]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
    omegaBA <- exp(b[11]) / (1+exp(b[14])+exp(b[11])+exp(b[12]))
    etaA <- exp(b[12]) / (1+exp(b[14])+exp(b[11])+exp(b[12]))
    
  }  
  else if(P_effects=="no prey effect")
  {  
    # POCHERS EXTINCTION independent of the presence of prey
   
    trans42 <- exp(b[14]) / (1+exp(b[13])+exp(b[14])+exp(b[15]))  # epsilonB-epsilonAB
    
    
    # POCHERS colonisation independent of the presence of prey
    gammaAB <- exp(b[6]) / (1+exp(b[4])+exp(b[5])+exp(b[6])) 
    trans12 <- exp(b[4]) / (1+exp(b[4])+exp(b[5])+exp(b[6]))  # gammaA-gammaAB
    etaA <- trans12 
    
    # PREYS colonisation and extinction dependent of presence of poacher
    epsilonAB <- exp(b[13]) / (1+exp(b[13])+exp(b[14])+exp(b[15]))
    trans43 <- exp(b[15]) / (1+exp(b[13])+exp(b[14])+exp(b[15])) # epsilonA-epsilonAB
    nuA <-  trans43 
    trans13 <- exp(b[5]) / (1+exp(b[4])+exp(b[5])+ exp(b[6])) # gammaB-gammaAB
    nuB <- exp(b[10]) / (1+exp(b[10])+ exp(b[11])+  exp(b[4]))  # to be checked with Oliver
    omegaBA <- exp(b[11]) / (1+exp(b[10])+exp(b[11]) + exp(b[4])) # to be checked with Oliver
    omegaAB <-  exp(b[8]) / (1+exp(b[9])+ exp(b[8]) + exp(b[15]))
    etaB <- exp(b[9]) / (1+exp(b[15])+exp(b[8])+exp(b[9]))
  
  }
  else if(P_effects=="species interaction")
  { 
    epsilonAB <- exp(b[13]) / (1+exp(b[13])+exp(b[14])+exp(b[15]))
    trans13 <- exp(b[5]) / (1+exp(b[4])+exp(b[5])+exp(b[6])) # gammaB-gammaAB
    gammaAB <- exp(b[6]) / (1+exp(b[4])+exp(b[5])+exp(b[6])) 
    trans42 <- exp(b[14]) / (1+exp(b[13])+exp(b[14])+exp(b[15]))  # epsilonB-epsilonAB
    trans43 <- exp(b[15]) / (1+exp(b[13])+exp(b[14])+exp(b[15])) 
    trans12 <- exp(b[4]) / (1+exp(b[4])+exp(b[5])+exp(b[6]))
    nuA <- exp(b[7]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
    omegaAB <- exp(b[8]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
    etaB <- exp(b[9]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
    nuB <- exp(b[10]) / (1+exp(b[10])+exp(b[11])+exp(b[12]))
    omegaBA <- exp(b[11]) / (1+exp(b[10])+exp(b[11])+exp(b[12]))
    etaA <- exp(b[12]) / (1+exp(b[10])+exp(b[11])+exp(b[12]))
    
    
  }
  else if(P_effects=="poachers performance")
  {
    # colonisation of poachers depend on the presence of prey
    
    trans12 <- exp(b[4]) / (1+exp(b[4])+exp(b[5])+exp(b[6])) # gammaA-gammaAB
    trans13 <- exp(b[5]) / (1+exp(b[4])+exp(b[5])+exp(b[6])) # gammaB-gammaAB
    gammaAB <- exp(b[6]) / (1+exp(b[4])+exp(b[5])+exp(b[6]))
    
    # colonisation of prey independent of poachers
    omegaAB <- exp(b[8]) / (1+exp(b[5])+exp(b[8])+exp(b[14]))
    etaB <- trans13 
    
    nuB <- exp(b[10]) / (1+exp(b[10])+exp(b[11])+exp(b[12]))
    omegaBA <- exp(b[11]) / (1+exp(b[10])+exp(b[11])+exp(b[12]))
    etaA <- exp(b[12]) / (1+exp(b[10])+exp(b[11])+exp(b[12])) 
    # extinction of poachers is independent of prey (poachers are doing a better job than ranger)
    epsilonAB <- exp(b[13]) / (1+exp(b[13])+exp(b[14])+exp(b[15]))
    trans42 <- exp(b[14]) / (1+exp(b[13])+exp(b[14])+exp(b[15])) # epsilonB-epsilonAB
    trans43 <- exp(b[15]) / (1+exp(b[13])+exp(b[14])+exp(b[15])) # epsilonA-epsilonAB
    nuA <-  trans43 
    
    
  }   
  else if(P_effects=="prey counteracting behaviour")
  {
    # POACHERS EXTINCTION dependent of the presence of prey
    epsilonAB <- exp(b[13]) / (1+exp(b[13])+exp(b[14])+exp(b[15]))
    trans42 <- exp(b[14]) / (1+exp(b[13])+exp(b[14])+exp(b[15]))  # epsilonB-epsilonAB
    omegaAB <- exp(b[8]) / (1+exp(b[8])+exp(b[7])+exp(b[9]))
    nuA <- exp(b[7]) / (1+exp(b[8])+exp(b[7])+exp(b[9]))
    
    # POACHERS colonisation independent of the presence of prey
    gammaAB <- exp(b[6]) / (1+exp(b[6])+exp(b[4])+exp(b[5])) 
    trans12 <- exp(b[4]) / (1+exp(b[6])+exp(b[4])+exp(b[5]))  # gammaA-gammaAB
    etaA <-  trans12
    
    # PREYS colonisation dependent of presence of poachers
    trans13 <- exp(b[5]) / (1+exp(b[6])+exp(b[5])+exp(b[4])) # gammaB-gammaAB
    etaB <- exp(b[9]) / (1+exp(b[8])+exp(b[7])+exp(b[9]))
    
    
    # PREYS extinction independent on occupancy of poachers
    trans43 <- exp(b[15]) / (1+exp(b[13])+exp(b[14])+exp(b[15])) # epsilonA-epsilonAB
    nuB <- trans42
    omegaBA <- exp(b[11]) / (1+exp(b[14])+exp(b[11])+exp(b[4]))
    
  }
  if(P_effects=="spatial effects of ranger on poachers")
  {  
    nuA<-matrix(0, nrow=R, ncol=N)
    omegaAB<-matrix(0, nrow=R, ncol=N)
    omegaBA<-matrix(0, nrow=R, ncol=N)
    etaA<-matrix(0, nrow=R, ncol=N)
    nuB<-matrix(0, nrow=R, ncol=N)
    etaB<-matrix(0, nrow=R, ncol=N)
    epsilonAB <- matrix(0, nrow=R, ncol=N)
    trans13 <- matrix(0, nrow=R, ncol=N) # gammaB-gammaAB
    gammaAB <- matrix(0, nrow=R, ncol=N)
    trans42 <- matrix(0, nrow=R, ncol=N) # epsilonB-epsilonAB
    trans43 <- matrix(0, nrow=R, ncol=N)
    trans12 <- matrix(0, nrow=R, ncol=N)
    
    for(s in 1:R)
    {
      epsilonAB[s,] <- exp(b[13] +  b[25]*tempcov[s,]) / (1+exp(b[13] +  b[25]*tempcov[s,])+exp(b[14] + b[26]*tempcov[s,])+exp(b[15]+ b[27]*tempcov[s,])) 
      trans13[s,] <- exp(b[5]) / (1+exp(b[4])+exp(b[5])+exp(b[6])) # gammaB-gammaAB
      gammaAB[s,] <- exp(b[6]) / (1+exp(b[4])+exp(b[5])+exp(b[6])) 
      trans42[s,] <- exp(b[14] + b[26]*tempcov[s,]) / (1+exp(b[13] +  b[25]*tempcov[s,])+exp(b[14] + b[26]*tempcov[s,]) + exp(b[15]+ b[27]*tempcov[s,]))  # epsilonB-epsilonAB
      trans43[s,] <- exp(b[15] + b[27]*tempcov[s,]) / (1+exp(b[13] +  b[25]*tempcov[s,])+exp(b[14] + b[26]*tempcov[s,]) + exp(b[15]+ b[27]*tempcov[s,])) 
      trans12[s,] <- exp(b[4] ) / (1+exp(b[4])+exp(b[5]) + exp(b[6]))
      
      
      nuA[s,] <- exp(b[7]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
      omegaAB[s,] <- exp(b[8]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
      etaB[s,] <- exp(b[9]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
      nuB[s,] <- exp(b[10]+ b[25]*tempcov[s,]) / (1+exp(b[10] + b[25]*tempcov[s,])+exp(b[11] + b[26]*tempcov[s,])+exp(b[12] + b[27]*tempcov[s,]))
      omegaBA[s,] <- exp(b[11] + b[26]*tempcov[s,]) / (1+exp(b[10] + b[25]*tempcov[s,])+exp(b[11]+ b[26]*tempcov[s,])+exp(b[12] + b[27]*tempcov[s,])) 
      etaA[s,] <- exp(b[12] + b[27]*tempcov[s,]) / (1+exp(b[10] + b[25]*tempcov[s,])+exp(b[11]+ b[26]*tempcov[s,])+exp(b[12] + b[27]*tempcov[s,]))
      
    }
  }
  if(P_effects=="temporal effects of ranger on poachers")
  {  
    nuA<-matrix(0, nrow=R, ncol=N)
    omegaAB<-matrix(0, nrow=R, ncol=N)
    omegaBA<-matrix(0, nrow=R, ncol=N)
    etaA<-matrix(0, nrow=R, ncol=N)
    nuB<-matrix(0, nrow=R, ncol=N)
    etaB<-matrix(0, nrow=R, ncol=N)
    epsilonAB <- matrix(0, nrow=R, ncol=N)
    trans13 <- matrix(0, nrow=R, ncol=N) # gammaB-gammaAB
    gammaAB <- matrix(0, nrow=R, ncol=N)
    trans42 <- matrix(0, nrow=R, ncol=N) # epsilonB-epsilonAB
    trans43 <- matrix(0, nrow=R, ncol=N)
    trans12 <- matrix(0, nrow=R, ncol=N)
    
    for(t in 1:(N-1))
    {
      epsilonAB[,t+1] <- exp(b[13] + b[25]*tempcov[,t]) / (1+exp(b[13] +  b[25]*tempcov[,t]) + exp(b[14] + b[26]*tempcov[,t]) + exp(b[15]+ b[27]*tempcov[,t])) 
      trans13[,t+1] <- exp(b[5]) / (1+exp(b[4])+exp(b[5])+exp(b[6])) # gammaB-gammaAB
      gammaAB[,t+1] <- exp(b[6]) / (1+exp(b[4])+exp(b[5])+exp(b[6])) 
      trans42[,t+1] <- exp(b[14] + b[26]*tempcov[,t]) / (1+exp(b[13] +  b[25]*tempcov[,t])+exp(b[14] + b[26]*tempcov[,t]) + exp(b[15]+ b[27]*tempcov[,t]))  # epsilonB-epsilonAB
      trans43[,t+1] <- exp(b[15] + b[27]*tempcov[,t]) / (1+exp(b[13] +  b[25]*tempcov[,t])+exp(b[14] + b[26]*tempcov[,t]) + exp(b[15]+ b[27]*tempcov[,t])) 
      trans12[,t+1] <- exp(b[4] ) / (1+exp(b[4])+exp(b[5]) + exp(b[6]))
      
    
      nuA[,t+1] <- exp(b[7]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
      omegaAB[,t+1] <- exp(b[8]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
      etaB[,t+1] <- exp(b[9]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
      nuB[,t+1] <- exp(b[10]+ b[15]*tempcov[,t]) / (1+exp(b[10] + b[25]*tempcov[,t])+exp(b[11] + b[26]*tempcov[,t])+exp(b[12] + b[27]*tempcov[,t]))
      omegaBA[,t+1] <- exp(b[11] + b[26]*tempcov[,t]) / (1+exp(b[10] + b[25]*tempcov[,t])+exp(b[11]+ b[26]*tempcov[,t])+exp(b[12] + b[27]*tempcov[,t])) 
      etaA[,t+1] <- exp(b[12] + b[27]*tempcov[,t]) / (1+exp(b[10] + b[25]*tempcov[,t])+exp(b[11]+ b[26]*tempcov[,t])+exp(b[12] + b[27]*tempcov[,t]))
      
    }
  }
  if(P_effects=="space")
  {  
    
    nuA<-matrix(0, nrow=R, ncol=N)
    omegaAB<-matrix(0, nrow=R, ncol=N)
    omegaBA<-matrix(0, nrow=R, ncol=N)
    etaA<-matrix(0, nrow=R, ncol=N)
    nuB<-matrix(0, nrow=R, ncol=N)
    etaB<-matrix(0, nrow=R, ncol=N)
    epsilonAB <- matrix(0, nrow=R, ncol=N)
    trans13 <- matrix(0, nrow=R, ncol=N) # gammaB-gammaAB
    gammaAB <- matrix(0, nrow=R, ncol=N)
    trans42 <- matrix(0, nrow=R, ncol=N) # epsilonB-epsilonAB
    trans43 <- matrix(0, nrow=R, ncol=N)
    trans12 <- matrix(0, nrow=R, ncol=N)
    
    for(s in 1:R)
    {
      epsilonAB[s,] <- exp(b[13] +  b[25]*tempcov[s,]) / (1+exp(b[13] +  b[25]*tempcov[s,])+exp(b[14] +  b[26]*tempcov[s,])+exp(b[15] +  b[27]*tempcov[s,]))
      trans13[s,] <- exp(b[5] +  b[26]*tempcov[s,]) / (1+exp(b[4] +  b[25]*tempcov[s,])+exp(b[5] +  b[26]*tempcov[s,])+exp(b[6] +  b[27]*tempcov[s,]))# gammaB-gammaAB
      gammaAB[s,] <- exp(b[6] +  b[27]*tempcov[s,]) / (1+exp(b[4] +  b[25]*tempcov[s,])+exp(b[5] +  b[26]*tempcov[s,])+exp(b[6] +  b[27]*tempcov[s,]))
      trans42[s,] <- exp(b[14] + b[26]*tempcov[s,]) / (1+exp(b[13] +  b[25]*tempcov[s,])+exp(b[14] + b[26]*tempcov[s,])+exp(b[15]+ b[27]*tempcov[s,]))  # epsilonB-epsilonAB
      trans43[s,] <- exp(b[15] + b[27]*tempcov[s,]) / (1+exp(b[13] +  b[25]*tempcov[s,])+exp(b[14] +  b[26]*tempcov[s,]) + exp(b[15] +  b[27]*tempcov[s,])) 
      trans12[s,] <- exp(b[4] +  b[25]*tempcov[s,]) / (1+exp(b[4] +  b[25]*tempcov[s,])+exp(b[5] +  b[26]*tempcov[s,])+exp(b[6] +  b[27]*tempcov[s,]))
      
      nuA[s,] <- exp(b[7] +  b[25]*tempcov[s,]) / (1+exp(b[7] + b[25]*tempcov[s,])+exp(b[8]+ b[26]*tempcov[s,])+exp(b[9] + b[27]*tempcov[s,]))
      omegaAB[s,] <- exp(b[8] +  b[26]*tempcov[s,]) / (1+exp(b[7] + b[25]*tempcov[s,])+exp(b[8]+ b[26]*tempcov[s,])+exp(b[9] + b[27]*tempcov[s,]))
      etaB[s,] <- exp(b[9] + b[27]*tempcov[s,]) / (1+exp(b[7] + b[25]*tempcov[s,])+exp(b[8]+ b[26]*tempcov[s,])+exp(b[9] + b[27]*tempcov[s,]))
      nuB[s,] <- exp(b[10] + b[25]*tempcov[s,]) / (1+exp(b[10] + b[25]*tempcov[s,])+exp(b[11]+ b[26]*tempcov[s,])+exp(b[12] + b[27]*tempcov[s,]))
      omegaBA[s,] <- exp(b[11] + b[26]*tempcov[s,]) / (1+exp(b[10] + b[25]*tempcov[s,])+exp(b[11]+ b[26]*tempcov[s,])+exp(b[12] + b[27]*tempcov[s,]))
      etaA[s,] <- exp(b[12] + b[27]*tempcov[s,]) / (1+exp(b[10] + b[25]*tempcov[s,])+exp(b[11]+ b[26]*tempcov[s,])+exp(b[12] + b[27]*tempcov[s,]))
      
    }
  }
  
  if(P_effects=="time cov")
  {  
    nuA<-matrix(0, nrow=R, ncol=N)
    omegaAB<-matrix(0, nrow=R, ncol=N)
    omegaBA<-matrix(0, nrow=R, ncol=N)
    etaA<-matrix(0, nrow=R, ncol=N)
    nuB<-matrix(0, nrow=R, ncol=N)
    etaB<-matrix(0, nrow=R, ncol=N)
    epsilonAB <- matrix(0, nrow=R, ncol=N)
    trans13 <- matrix(0, nrow=R, ncol=N) # gammaB-gammaAB
    gammaAB <- matrix(0, nrow=R, ncol=N)
    trans42 <- matrix(0, nrow=R, ncol=N) # epsilonB-epsilonAB
    trans43 <- matrix(0, nrow=R, ncol=N)
    trans12 <- matrix(0, nrow=R, ncol=N)
    
    for(t in 1:(N-1))
    {
      epsilonAB[,t+1] <- exp(b[13] +  b[25]*tempcov[,t]) / (1+exp(b[13] +  b[25]*tempcov[,t])+exp(b[14] +  b[26]*tempcov[,t])+exp(b[15] +  b[27]*tempcov[,t]))
      trans13[,t+1] <- exp(b[5] +  b[26]*tempcov[,t]) / (1+exp(b[4] +  b[25]*tempcov[,t])+exp(b[5] +  b[26]*tempcov[,t])+exp(b[6] +  b[27]*tempcov[,t]))# gammaB-gammaAB
      gammaAB[,t+1] <- exp(b[6] +  b[27]*tempcov[,t]) / (1+exp(b[4] +  b[25]*tempcov[,t])+exp(b[5] +  b[26]*tempcov[,t])+exp(b[6] +  b[27]*tempcov[,t]))
      trans42[,t+1] <- exp(b[14] + b[26]*tempcov[,t]) / (1+exp(b[13] +  b[25]*tempcov[,t])+exp(b[14] + b[26]*tempcov[,t])+exp(b[15]+ b[27]*tempcov[,t]))  # epsilonB-epsilonAB
      trans43[,t+1] <- exp(b[15] + b[27]*tempcov[,t]) / (1+exp(b[13] +  b[25]*tempcov[,t])+exp(b[14] +  b[26]*tempcov[,t]) + exp(b[15] +  b[27]*tempcov[,t])) 
      trans12[,t+1] <- exp(b[4] +  b[25]*tempcov[,t]) / (1+exp(b[4] +  b[25]*tempcov[,t])+exp(b[5] +  b[26]*tempcov[,t])+exp(b[6] +  b[27]*tempcov[,t]))
      
      nuA[,t+1] <- exp(b[7] +  b[25]*tempcov[,t]) / (1+exp(b[7] + b[25]*tempcov[,t])+exp(b[8]+ b[26]*tempcov[,t])+exp(b[9] + b[27]*tempcov[,t]))
      omegaAB[,t+1] <- exp(b[8] +  b[26]*tempcov[,t]) / (1+exp(b[7] + b[25]*tempcov[,t])+exp(b[8]+ b[26]*tempcov[,t])+exp(b[9] + b[27]*tempcov[,t]))
      etaB[,t+1] <- exp(b[9] + b[27]*tempcov[,t]) / (1+exp(b[7] + b[25]*tempcov[,t])+exp(b[8]+ b[26]*tempcov[,t])+exp(b[9] + b[27]*tempcov[,t]))
      nuB[,t+1] <- exp(b[10] + b[25]*tempcov[,t]) / (1+exp(b[10] + b[25]*tempcov[,t])+exp(b[11]+ b[26]*tempcov[,t])+exp(b[12] + b[27]*tempcov[,t]))
      omegaBA[,t+1] <- exp(b[11] + b[26]*tempcov[,t]) / (1+exp(b[10] + b[25]*tempcov[,t])+exp(b[11]+ b[26]*tempcov[,t])+exp(b[12] + b[27]*tempcov[,t]))
      etaA[,t+1] <- exp(b[12] + b[27]*tempcov[,t]) / (1+exp(b[10] + b[25]*tempcov[,t])+exp(b[11]+ b[26]*tempcov[,t])+exp(b[12] + b[27]*tempcov[,t]))
    }
  }
  
  if(P_effects=="both delay")
  {  
    nuA<-matrix(0, nrow=R, ncol=N)
    omegaAB<-matrix(0, nrow=R, ncol=N)
    omegaBA<-matrix(0, nrow=R, ncol=N)
    etaA<-matrix(0, nrow=R, ncol=N)
    nuB<-matrix(0, nrow=R, ncol=N)
    etaB<-matrix(0, nrow=R, ncol=N)
    epsilonAB <- matrix(0, nrow=R, ncol=N)
    trans13 <- matrix(0, nrow=R, ncol=N) # gammaB-gammaAB
    gammaAB <- matrix(0, nrow=R, ncol=N)
    trans42 <- matrix(0, nrow=R, ncol=N) # epsilonB-epsilonAB
    trans43 <- matrix(0, nrow=R, ncol=N)
    trans12 <- matrix(0, nrow=R, ncol=N)
    
    
    for(s in 1:R)
    {   for(t in 1:(N-1))
    {
      
      epsilonAB[s,t+1] <- exp(b[13] +  b[28]*tempcov[s,t]) / (1+exp(b[13] +  b[28]*tempcov[s,t])+exp(b[14] +  b[29]*tempcov[s,t])+exp(b[15] +  b[30]*tempcov[s,t]))
      trans13[s,t+1] <- exp(b[5] +  b[32]*tempcov[s,t]) / (1+exp(b[4] +  b[31]*tempcov[s,t])+exp(b[5] +  b[32]*tempcov[s,t])+exp(b[6] +  b[33]*tempcov[s,t]))# gammaB-gammaAB
      gammaAB[s,t+1] <- exp(b[6] +  b[33]*tempcov[s,t]) / (1+exp(b[4] +  b[31]*tempcov[s,t])+exp(b[5] +  b[32]*tempcov[s,t])+exp(b[6] +  b[33]*tempcov[s,t]))
      trans42[s,t+1] <- exp(b[14] + b[29]*tempcov[s,t]) / (1+exp(b[13] +  b[28]*tempcov[s,t])+exp(b[14] + b[29]*tempcov[s,t])+exp(b[15]+ b[30]*tempcov[s,t]))  # epsilonB-epsilonAB
      trans43[s,t+1] <- exp(b[15] + b[30]*tempcov[s,t]) / (1+exp(b[13] +  b[28]*tempcov[s,t])+exp(b[14] +  b[29]*tempcov[s,t]) + exp(b[15] +  b[30]*tempcov[s,t])) 
      trans12[s,t+1] <- exp(b[4] +  b[31]*tempcov[s,t]) / (1+exp(b[4] +  b[31]*tempcov[s,t])+exp(b[5] +  b[32]*tempcov[s,t])+exp(b[6] +  b[33]*tempcov[s,t]))
      nuA[s,t+1] <- exp(b[7] +  b[25]*tempcov[s,t]) / (1+exp(b[7] + b[25]*tempcov[s,t])+exp(b[8]+ b[26]*tempcov[s,t])+exp(b[9] + b[27]*tempcov[s,t]))
      omegaAB[s,t+1] <- exp(b[8] +  b[26]*tempcov[s,t]) / (1+exp(b[7] + b[25]*tempcov[s,t])+exp(b[8]+ b[26]*tempcov[s,t])+exp(b[9] + b[27]*tempcov[s,t]))
      etaB[s,t+1] <- exp(b[9] + b[27]*tempcov[s,t]) / (1+exp(b[7] + b[25]*tempcov[s,t])+exp(b[8]+ b[26]*tempcov[s,t])+exp(b[9] + b[27]*tempcov[s,t]))
      nuB[s,t+1] <- exp(b[10] + b[33]*tempcov[s,t]) / (1+exp(b[10] + b[33]*tempcov[s,t])+exp(b[11]+ b[34]*tempcov[s,t])+exp(b[12] + b[35]*tempcov[s,t]))
      omegaBA[s,t+1] <- exp(b[11] + b[34]*tempcov[s,t]) / (1+exp(b[10] + b[33]*tempcov[s,t])+exp(b[11]+ b[34]*tempcov[s,t])+exp(b[12] + b[35]*tempcov[s,t]))
      etaA[s,t+1] <- exp(b[12] + b[35]*tempcov[s,t]) / (1+exp(b[10] + b[33]*tempcov[s,t])+exp(b[11]+ b[34]*tempcov[s,t])+exp(b[12] + b[35]*tempcov[s,t]))
      
    }
    }
  }
  
  if(P_effects=="both")
  {  
    nuA<-matrix(0, nrow=R, ncol=N)
    omegaAB<-matrix(0, nrow=R, ncol=N)
    omegaBA<-matrix(0, nrow=R, ncol=N)
    etaA<-matrix(0, nrow=R, ncol=N)
    nuB<-matrix(0, nrow=R, ncol=N)
    etaB<-matrix(0, nrow=R, ncol=N)
    epsilonAB <- matrix(0, nrow=R, ncol=N)
    trans13 <- matrix(0, nrow=R, ncol=N) # gammaB-gammaAB
    gammaAB <- matrix(0, nrow=R, ncol=N)
    trans42 <- matrix(0, nrow=R, ncol=N) # epsilonB-epsilonAB
    trans43 <- matrix(0, nrow=R, ncol=N)
    trans12 <- matrix(0, nrow=R, ncol=N)
    
    
    for(s in 1:R)
    {
      for(t in 1: N)
     {
        epsilonAB[s,t] <- exp(b[13] +  b[28]*tempcov[s,t]) / (1+exp(b[13] +  b[28]*tempcov[s,t])+exp(b[14] +  b[29]*tempcov[s,t])+exp(b[15] +  b[30]*tempcov[s,t]))
        trans13[s,t] <- exp(b[5] +  b[32]*tempcov[s,t]) / (1+exp(b[4] +  b[31]*tempcov[s,t])+exp(b[5] +  b[32]*tempcov[s,t])+exp(b[6] +  b[33]*tempcov[s,t]))# gammaB-gammaAB
        gammaAB[s,t] <- exp(b[6] +  b[33]*tempcov[s,t]) / (1+exp(b[4] +  b[31]*tempcov[s,t])+exp(b[5] +  b[32]*tempcov[s,t])+exp(b[6] +  b[33]*tempcov[s,t]))
        trans42[s,t] <- exp(b[14] + b[29]*tempcov[s,t]) / (1+exp(b[13] +  b[28]*tempcov[s,t])+exp(b[14] + b[29]*tempcov[s,t])+exp(b[15]+ b[30]*tempcov[s,t]))  # epsilonB-epsilonAB
        trans43[s,t] <- exp(b[15] + b[30]*tempcov[s,t]) / (1+exp(b[13] +  b[28]*tempcov[s,t])+exp(b[14] +  b[29]*tempcov[s,t]) + exp(b[15] +  b[30]*tempcov[s,t])) 
        trans12[s,t] <- exp(b[4] +  b[31]*tempcov[s,t]) / (1+exp(b[4] +  b[31]*tempcov[s,t])+exp(b[5] +  b[32]*tempcov[s,t])+exp(b[6] +  b[33]*tempcov[s,t]))
        nuA[s,t] <- exp(b[7] +  b[25]*tempcov[s,t]) / (1+exp(b[7] + b[25]*tempcov[s,t])+exp(b[8]+ b[26]*tempcov[s,t])+exp(b[9] + b[27]*tempcov[s,t]))
        omegaAB[s,t] <- exp(b[8] +  b[26]*tempcov[s,t]) / (1+exp(b[7] + b[25]*tempcov[s,t])+exp(b[8]+ b[26]*tempcov[s,t])+exp(b[9] + b[27]*tempcov[s,t]))
        etaB[s,t] <- exp(b[9] + b[27]*tempcov[s,t]) / (1+exp(b[7] + b[25]*tempcov[s,t])+exp(b[8]+ b[26]*tempcov[s,t])+exp(b[9] + b[27]*tempcov[s,t]))
        nuB[s,t] <- exp(b[10] + b[33]*tempcov[s,t]) / (1+exp(b[10] + b[33]*tempcov[s,t])+exp(b[11]+ b[34]*tempcov[s,t])+exp(b[12] + b[35]*tempcov[s,t]))
        omegaBA[s,t] <- exp(b[11] + b[34]*tempcov[s,t]) / (1+exp(b[10] + b[33]*tempcov[s,t])+exp(b[11]+ b[34]*tempcov[s,t])+exp(b[12] + b[35]*tempcov[s,t]))
        etaA[s,t] <- exp(b[12] + b[35]*tempcov[s,t]) / (1+exp(b[10] + b[33]*tempcov[s,t])+exp(b[11]+ b[34]*tempcov[s,t])+exp(b[12] + b[35]*tempcov[s,t]))
        
      }
      
    }
  }
  
 
  
  return(list(epsilonAB, trans42, nuA, omegaAB, trans12, etaA,  gammaAB, trans43, nuB,  omegaBA, trans13, etaB))
  
}
