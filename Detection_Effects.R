EffectsDetection<-function(P_effects, b)
{
  # various quantities that will be useful later on

  
  if(P_effects=="cst detection")
  {  
    
    pB <- 1/(1+exp(-(b[16])))
    pAb <- 1/(1+exp(-(b[16])))
    pAB <- 1/(1+exp(-(b[16])))
    
    B <- matrix(c(
      1,0,0,0,
      1-pAb,pAb,0,0,
      1-pB,0,pB,0,
      1-pAB,0,0,pAB),
      nrow = n.states)
    npar<-16 # Define the number of regression coefficient for this effect
  }  
  else if(P_effects=="no prey effects on detection of poachers")
  {  
    
    pB <- 1/(1+exp(-(b[16])))
    pAb <- 1/(1+exp(-(b[17])))
    pAB <- 1/(1+exp(-(b[17])))
    
    B <- matrix(c(
      1,0,0,0,
      1-pAb,pAb,0,0,
      1-pB,0,pB,0,
      1-pAB,0,0,pAB),
      nrow = n.states)
    npar<-17 # Define the number of regression coefficient for this effect
  }  
  else if(P_effects=="no poacher effects on detection of preys")
  {  
    
    pB <- 1/(1+exp(-(b[16])))
    pAb <- 1/(1+exp(-(b[17])))
    pAB <- 1/(1+exp(-(b[16])))
    
    B <- matrix(c(
      1,0,0,0,
      1-pAb,pAb,0,0,
      1-pB,0,pB,0,
      1-pAB,0,0,pAB),
      nrow = n.states)
    npar<-17 # Define the number of regression coefficient for this effect
  }  
  else if(P_effects=="species interaction effects on detection")
  {  
    
    pB <- 1/(1+exp(-(b[16])))
    pAb <- 1/(1+exp(-(b[17])))
    pAB <- 1/(1+exp(-(b[18])))
    
    B <- matrix(c(
      1,0,0,0,
      1-pAb,pAb,0,0,
      1-pB,0,pB,0,
      1-pAB,0,0,pAB),
      nrow = n.states)
    npar<-18 # Define the number of regression coefficient for this effect
  }  
  else if(P_effects=="spatial effect of patroling effort on detection")
  {  
    pB <-vector(mode="logical", length=R)
    pAb <-vector(mode="logical", length=R)
    pAB <-vector(mode="logical", length=R)
    B <-array(0, dim=c(n.states,n.states,R))
    
    for(s in 1:R)
    {    
      pB[s] <- 1/(1+exp(-(b[16] + b[19]*mean(tempcov[s,]))))
      pAb[s] <- 1/(1+exp(-(b[17] + b[20]*mean(tempcov[s,]))))
      pAB[s] <- 1/(1+exp(-(b[18] + b[21]*mean(tempcov[s,]))))
      
      B[,,s] <- matrix(c(1,0,0,0,
                         1-pAb[s],pAb[s],0,0,
                         1-pB[s],0,pB[s],0,
                         1-pAB[s],0,0,pAB[s]),
                       nrow = n.states)
    }
    npar<-21 # Define the number of regression coefficient for this effect
    
  } else if(P_effects=="time effect of patroling effort on detection")
  {  
    pB <-vector(mode="logical", length=N)
    pAb <-vector(mode="logical", length=N)
    pAB <-vector(mode="logical", length=N)
    # obs prob
    B <-array(0, dim=c(n.states,n.states,N))
    
    for(n in 1:N)
    {    
      pB[n] <- 1/(1+exp(-(b[16] + b[19]*mean(tempcov[,n]))))
      pAb[n] <- 1/(1+exp(-(b[17] + b[20]*mean(tempcov[,n]))))
      pAB[n] <- 1/(1+exp(-(b[18] + b[21]*mean(tempcov[,n]))))
      
      
      B[,,n] <- matrix(c(
        1,0,0,0,
        1-pAb[n],pAb[n],0,0,
        1-pB[n],0,pB[n],0,
        1-pAB[n],0,0,pAB[n]),
        nrow = n.states)
      
      npar<-21 # Define the number of regression coefficient for this effect
      
    }
    
    
  } else if(P_effects=="spatio-temporal effect of patroling effort on detection")
  {  
    pB <-matrix(0, nrow=R, ncol=N)
    pAb <-matrix(0, nrow=R, ncol=N)
    pAB <-matrix(0, nrow=R, ncol=N)
    # obs prob
    B <-array(0, dim=c(n.states,n.states,R,N))
    
    for(s in 1:R)
    { 
      for(n in 1:N)
      {    
        pB[s,n] <- 1/(1+exp(-(b[16] + b[19]*tempcov[s,n])))
        pAb[s,n] <- 1/(1+exp(-(b[17] + b[20]*tempcov[s,n])))
        pAB[s,n] <- 1/(1+exp(-(b[18] + b[21]*tempcov[s,n])))
        
        B[,,s,n] <- matrix(c(
          1,0,0,0,
          1-pAb[s,n],pAb[s,n],0,0,
          1-pB[s,n],0,pB[s,n],0,
          1-pAB[s,n],0,0,pAB[s,n]),
          nrow = n.states)
      }
    }  
    npar<-21 # Define the number of regression coefficient for this effect
    
  }  
  else if(P_effects=="interactive effects of years")
  {  
      pB <-vector(mode="logical", length=K)
      pAb <-vector(mode="logical", length=K)
      pAB <-vector(mode="logical", length=K)
      B <-array(0, dim=c(n.states,n.states,K))
      pB[1] <- 1/(1+exp(-(b[16])))
      pAb[1] <- 1/(1+exp(-(b[17])))
      pAB[1] <- 1/(1+exp(-(b[18])))
      
      for(k in 2:K)
      {  
        pB[k] <- 1/(1+exp(-(b[(19 + k)])))
        pAb[k] <- 1/(1+exp(-(b[(20 + (K-1) + k)])))
        pAB[k] <- 1/(1+exp(-(b[(21 + 2*(K-1) + k)])))
        
        B[,,k] <- matrix(c(
          1,0,0,0,
          1-pAb[k],pAb[k],0,0,
          1-pB[k],0,pB[k],0,
          1-pAB[k],0,0,pAB[k]),
          nrow = n.states)
      }
      npar <- 20 + 3*(K) # Define the number of regression coefficient for this effect
      
  } 
  
  else if(P_effects=="additive effects of years")
  {  
    pB <-vector(mode="logical", length=K)
    pAb <-vector(mode="logical", length=K)
    pAB <-vector(mode="logical", length=K)
    B <-array(0, dim=c(n.states,n.states,K))
    pB[1] <- 1/(1+exp(-(b[16])))
    pAb[1] <- 1/(1+exp(-(b[17])))
    pAB[1] <- 1/(1+exp(-(b[18])))
    
    for(k in 2:K)
    {  
      pB[k] <- 1/(1+exp(-(b[(19  + k)])))
      pAb[k] <- 1/(1+exp(-(b[(19 + k)])))
      pAB[k] <- 1/(1+exp(-(b[(19 + k)])))
      
      B[,,k] <- matrix(c(
        1,0,0,0,
        1-pAb[k],pAb[k],0,0,
        1-pB[k],0,pB[k],0,
        1-pAB[k],0,0,pAB[k]),
        nrow = n.states)
    }
    npar <- 20 + 3*(K) # Define the number of regression coefficient for this effect
    
  }
    
  
  else if(P_effects=="distance to road effect on detection")
  {  
    pB <-vector(mode="logical", length=R)
    pAb <-vector(mode="logical", length=R)
    pAB <-vector(mode="logical", length=R)
    B <-array(0, dim=c(n.states,n.states,R))
    
    for(s in 1:R)
    {    
      pB[s] <- 1/(1+exp(-(b[16] + b[19]*roadcov[s])))
      pAb[s] <- 1/(1+exp(-(b[17] + b[20]*roadcov[s])))
      pAB[s] <- 1/(1+exp(-(b[18] + b[21]*roadcov[s])))
      
      B[,,s] <- matrix(c(1,0,0,0,
                         1-pAb[s],pAb[s],0,0,
                         1-pB[s],0,pB[s],0,
                         1-pAB[s],0,0,pAB[s]),
                       nrow = n.states)
    }
    npar<-21 # Define the number of regression coefficient for this effect
    
  }  
  
  else if(P_effects=="stream length effect on detection")
  {  
    pB <-vector(mode="logical", length=R)
    pAb <-vector(mode="logical", length=R)
    pAB <-vector(mode="logical", length=R)
    B <-array(0, dim=c(n.states,n.states,R))
    
    for(s in 1:R)
    {    
      pB[s] <- 1/(1+exp(-(b[16] + b[19]*streamcov[s])))
      pAb[s] <- 1/(1+exp(-(b[17] + b[20]*streamcov[s])))
      pAB[s] <- 1/(1+exp(-(b[18] + b[21]*streamcov[s])))
      
      B[,,s] <- matrix(c(1,0,0,0,
                         1-pAb[s],pAb[s],0,0,
                         1-pB[s],0,pB[s],0,
                         1-pAB[s],0,0,pAB[s]),
                       nrow = n.states)
    }
    npar<-21 # Define the number of regression coefficient for this effect
    
  }  
  
  
  else if(P_effects=="distance to station effect on detection")
  {  
    pB <-vector(mode="logical", length=R)
    pAb <-vector(mode="logical", length=R)
    pAB <-vector(mode="logical", length=R)
    B <-array(0, dim=c(n.states,n.states,R))
    
    for(s in 1:R)
    {    
      pB[s] <- 1/(1+exp(-(b[16] + b[19]*stationcov[s])))
      pAb[s] <- 1/(1+exp(-(b[17] + b[20]*stationcov[s])))
      pAB[s] <- 1/(1+exp(-(b[18] + b[21]*stationcov[s])))
      
      B[,,s] <- matrix(c(1,0,0,0,
                         1-pAb[s],pAb[s],0,0,
                         1-pB[s],0,pB[s],0,
                         1-pAB[s],0,0,pAB[s]),
                       nrow = n.states)
    }
    npar<-21 # Define the number of regression coefficient for this effect
    
  }  
  
  
  return(list(pB, pAb, pAB, B, npar))
  
}  
