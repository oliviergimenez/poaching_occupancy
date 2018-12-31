EffectsOccupancy<-function(P_effects, b)
{  
  
  if(P_effects=="no prey effect on occupancy of poachers")
  {  
    # In this section we assume that occupancy of poachers is independent of the occupancy of prey. In other word we describe the presence of poachers in a site regardless of whether or not it is occupied by preys
    # In the most general model occupancy estimates or obtained from the three first regression coefficients on logit scale             so to indicate the absence of effect of the presence of prey on the occupancy of poacher we set the probaility that a site occupied by both species is equal to the one of a poacher occupies a site in absence of prey defined which is defined as a logistic function of the second regression coefficient.  
    psiB <- 1/(1+exp(-(b[1])))
    psiAb <- 1/(1+exp(-(b[2])))
    psiAB <- 1/(1+exp(-(b[2])))
    
    # here we define the matrix with the intial state regarding the occupation of site appropriate to the effect considered
    PI1 <- array(data=NA, dim=c(1,2))
    PI2 <- array(data=NA, dim=c(2,4))
    PI <- array(data=NA, dim=c(1,4))
    
    PI1 <- c(1-psiB,psiB)
    PI2 <- matrix(c(1-psiAb,psiAb,0,0,0,0,1-psiAB,psiAB),nrow=2,byrow=T) # 
    PI <- PI1 %*% PI2
    
  }
  else if(P_effects=="species interaction effects on occupancy")
  {  
    # Here we assume occupancy probability of paochers to vary according to the presence of prey. Therefore we have a regression coefficient associated to each proability of occupancy, the three first ones of the list defined in the general model
    psiB <- 1/(1+exp(-(b[1])))
    psiAb <- 1/(1+exp(-(b[2])))
    psiAB <- 1/(1+exp(-(b[3])))
    
    # here we define the matrix of the intial states occupancy. We set the dimension of the matrix according to the tested hypothesis
    
    PI1 <- array(data=NA, dim=c(1,2))
    PI2 <- array(data=NA, dim=c(2,4))
    PI <- array(data=NA, dim=c(1,4))
    
    PI1 <- c(1-psiB,psiB)
    PI2 <- matrix(c(1-psiAb,psiAb,0,0,0,0,1-psiAB,psiAB),nrow=2,byrow=T) # 
    PI <- PI1 %*% PI2
    
  }
  
  if(P_effects=="spatial effect of patroling effort on occupancy")
  {  
    
    # Here we are looking if the probability of species occupancy is site dependent and moreprecisely vary with the spatial distribution of patrolling effort 
    psiB <-vector(mode="logical", length=R)
    psiAb <-vector(mode="logical", length=R)
    psiAB <-vector(mode="logical", length=R)
    
    # here we define the matrix of the intial states occupancy. We set the dimension of the matrix according to the tested hypothesis. In this case it is a three dimension matrix one for the occupancy of prey the other for occupancy of poacher given the presence absence of prey and the third the site which occupation state depend on the patrolling effort 
    PI1 <- array(data=NA, dim=c(1,2,R))
    PI2 <- array(data=NA, dim=c(2,4,R))
    PI <- array(data=NA, dim=c(1,4,R))
    
    for(s in 1:R)
    {
      psiB[s] <-  1/(1+exp(-(b[1]  + b[22]*mean(tempcov[s,]))))
      psiAb[s] <- 1/(1+exp(-(b[2] + b[23]*mean(tempcov[s,]))))
      psiAB[s] <- 1/(1+exp(-(b[3] + b[24]*mean(tempcov[s,]))))
      
      
      PI1[,,s] <- c(1-psiB[s],psiB[s])
      PI2[,,s] <- matrix(c(1-psiAb[s],psiAb[s],0,0,0,0,1-psiAB[s],psiAB[s]),nrow=2,byrow=T) # 
      PI[,,s] <- PI1[,,s] %*% PI2[,,s]
    }
    
  }
  
  if(P_effects=="spatial effect of patroling effort on poacher occupancy only")
  {
    # Here we are looking if the probability of species occupancy is site dependent and moreprecisely vary with the spatial distribution of patrolling effort 
    psiB <-vector(mode="logical", length=R)
    psiAb <-vector(mode="logical", length=R)
    psiAB <-vector(mode="logical", length=R)
    
    # here we define the matrix of the intial states occupancy. We set the dimension of the matrix according to the tested hypothesis. In this case it is a three dimension matrix one for the occupancy of prey the other for occupancy of poacher given the presence absence of prey and the third the site which occupation state depend on the patrolling effort 
    PI1 <- array(data=NA, dim=c(1,2,R))
    PI2 <- array(data=NA, dim=c(2,4,R))
    PI <- array(data=NA, dim=c(1,4,R))
    
    for(s in 1:R)
    {
      psiB[s] <- 1/(1+exp(-(b[1])))
      psiAb[s] <- 1/(1+exp(-(b[2] + b[22]*mean(tempcov[s,]))))
      psiAB[s] <- 1/(1+exp(-(b[2] + b[22]*mean(tempcov[s,]))))
      
      
      PI1[,,s] <- c(1-psiB[s],psiB[s])
      PI2[,,s] <- matrix(c(1-psiAb[s],psiAb[s],0,0,0,0,1-psiAB[s],psiAB[s]),nrow=2,byrow=T) # 
      PI[,,s] <- PI1[,,s] %*% PI2[,,s]
    }
    
  }
     
  if(P_effects=="spatial effect of patroling effort on poacher occupancy given presence/absence of prey")
  {
    # Here we are looking if the probability of species occupancy is site dependent and moreprecisely vary with the spatial distribution of patrolling effort 
    psiB <-vector(mode="logical", length=R)
    psiAb <-vector(mode="logical", length=R)
    psiAB <-vector(mode="logical", length=R)
    
    # here we define the matrix of the intial states occupancy. We set the dimension of the matrix according to the tested hypothesis. In this case it is a three dimension matrix one for the occupancy of prey the other for occupancy of poacher given the presence absence of prey and the third the site which occupation state depend on the patrolling effort 
    PI1 <- array(data=NA, dim=c(1,2,R))
    PI2 <- array(data=NA, dim=c(2,4,R))
    PI <- array(data=NA, dim=c(1,4,R))
    
    for(s in 1:R)
    {
      psiB[s] <- 1/(1+exp(-(b[1])))
      psiAb[s] <- 1/(1+exp(-(b[2] + b[22]*mean(tempcov[s,]))))
      psiAB[s] <- 1/(1+exp(-(b[3] + b[23]*mean(tempcov[s,]))))
      
      
      PI1[,,s] <- c(1-psiB[s],psiB[s])
      PI2[,,s] <- matrix(c(1-psiAb[s],psiAb[s],0,0,0,0,1-psiAB[s],psiAB[s]),nrow=2,byrow=T) # 
      PI[,,s] <- PI1[,,s] %*% PI2[,,s]
    }
    
  }
  
  if(P_effects=="distance to road on occupancy")
  {  
    
    # Here we are looking if the probability of species occupancy is site dependent and moreprecisely vary with the spatial distribution of patrolling effort 
    psiB <-vector(mode="logical", length=R)
    psiAb <-vector(mode="logical", length=R)
    psiAB <-vector(mode="logical", length=R)
    
    # here we define the matrix of the intial states occupancy. We set the dimension of the matrix according to the tested hypothesis. In this case it is a three dimension matrix one for the occupancy of prey the other for occupancy of poacher given the presence absence of prey and the third the site which occupation state depend on the patrolling effort 
    PI1 <- array(data=NA, dim=c(1,2,R))
    PI2 <- array(data=NA, dim=c(2,4,R))
    PI <- array(data=NA, dim=c(1,4,R))
    
    for(s in 1:R)
    {
      psiB[s] <-  1/(1+exp(-(b[1]  + b[22]*roadcov[s])))
      psiAb[s] <- 1/(1+exp(-(b[2] + b[23]*roadcov[s])))
      psiAB[s] <- 1/(1+exp(-(b[3] + b[24]*roadcov[s])))
      
      
      PI1[,,s] <- c(1-psiB[s],psiB[s])
      PI2[,,s] <- matrix(c(1-psiAb[s],psiAb[s],0,0,0,0,1-psiAB[s],psiAB[s]),nrow=2,byrow=T) # 
      PI[,,s] <- PI1[,,s] %*% PI2[,,s]
    }
    
  }
  
  
  if(P_effects=="distance to station effect on occupancy")
  {  
    
    # Here we are looking if the probability of species occupancy is site dependent and moreprecisely vary with the spatial distribution of patrolling effort 
    psiB <-vector(mode="logical", length=R)
    psiAb <-vector(mode="logical", length=R)
    psiAB <-vector(mode="logical", length=R)
    
    # here we define the matrix of the intial states occupancy. We set the dimension of the matrix according to the tested hypothesis. In this case it is a three dimension matrix one for the occupancy of prey the other for occupancy of poacher given the presence absence of prey and the third the site which occupation state depend on the patrolling effort 
    PI1 <- array(data=NA, dim=c(1,2,R))
    PI2 <- array(data=NA, dim=c(2,4,R))
    PI <- array(data=NA, dim=c(1,4,R))
    
    for(s in 1:R)
    {
      psiB[s] <-  1/(1+exp(-(b[1]  + b[22]*stationcov[s])))
      psiAb[s] <- 1/(1+exp(-(b[2] + b[23]*stationcov[s])))
      psiAB[s] <- 1/(1+exp(-(b[3] + b[24]*stationcov[s])))
      
      
      PI1[,,s] <- c(1-psiB[s],psiB[s])
      PI2[,,s] <- matrix(c(1-psiAb[s],psiAb[s],0,0,0,0,1-psiAB[s],psiAB[s]),nrow=2,byrow=T) # 
      PI[,,s] <- PI1[,,s] %*% PI2[,,s]
    }
    
  }
  
  if(P_effects=="stream length effect on occupancy")
  {  
    
    # Here we are looking if the probability of species occupancy is site dependent and moreprecisely vary with the spatial distribution of patrolling effort 
    psiB <-vector(mode="logical", length=R)
    psiAb <-vector(mode="logical", length=R)
    psiAB <-vector(mode="logical", length=R)
    
    # here we define the matrix of the intial states occupancy. We set the dimension of the matrix according to the tested hypothesis. In this case it is a three dimension matrix one for the occupancy of prey the other for occupancy of poacher given the presence absence of prey and the third the site which occupation state depend on the patrolling effort 
    PI1 <- array(data=NA, dim=c(1,2,R))
    PI2 <- array(data=NA, dim=c(2,4,R))
    PI <- array(data=NA, dim=c(1,4,R))
    
    for(s in 1:R)
    {
      psiB[s] <-  1/(1+exp(-(b[1]  + b[22]*streamcov[s])))
      psiAb[s] <- 1/(1+exp(-(b[2] + b[23]*streamcov[s])))
      psiAB[s] <- 1/(1+exp(-(b[3] + b[24]*streamcov[s])))
      
      
      PI1[,,s] <- c(1-psiB[s],psiB[s])
      PI2[,,s] <- matrix(c(1-psiAb[s],psiAb[s],0,0,0,0,1-psiAB[s],psiAB[s]),nrow=2,byrow=T) # 
      PI[,,s] <- PI1[,,s] %*% PI2[,,s]
    }
    
  }
  
  
  
  return(list(psiB, psiAb, psiAB, PI))
}

