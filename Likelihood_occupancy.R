dev_2spocc_dyn2 <- function(b,data,eff, tempcov, roadcov, stationcov, streamcov, garb,nh,primary,secondary,nstates, P_effects){
  
  R <- nh
  # b = vector of parameters on real scale
  # data = site histories
  # eff = nb of sites with that particular history
  # garb = initial states
  # nh = nb of sites
  # nstates = nb of occupancy states
  
  #---------------------------------------------
  # apply multinomial logit and standard logit link
  
  #---------------------------------------------
  # various quantities that will be useful later on
  K <- length(primary)
  J2 <- length(secondary)
  J <- J2/K
  N <- J * K
  #---------------------------------------------
  
  
  #   
  psiB <- EffectsOccupancy(P_effects,b) [[1]]
  psiAb <- EffectsOccupancy(P_effects,b) [[2]]
  psiAB <- EffectsOccupancy(P_effects,b) [[3]]
  PI <-  EffectsOccupancy(P_effects,b) [[4]]
  
  
  
  gammaAB <- exp(b[6]) / (1+exp(b[4])+exp(b[5])+exp(b[6]))
  trans12 <- exp(b[4]) / (1+exp(b[4])+exp(b[5])+exp(b[6])) # gammaA-gammaAB
  trans13 <- exp(b[5]) / (1+exp(b[4])+exp(b[5])+exp(b[6])) # gammaB-gammaAB
  
  
  gammaA = gammaAB + trans12
  gammaB = gammaAB + trans13
  nuA <- exp(b[7]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
  
  omegaAB <- exp(b[8]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
  etaB <- exp(b[9]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
  
  nuB <- exp(b[10]) / (1+exp(b[10])+exp(b[11])+exp(b[12]))
  omegaBA <- exp(b[11]) / (1+exp(b[10])+exp(b[11])+exp(b[12]))
  etaA <- exp(b[12]) / (1+exp(b[10])+exp(b[11])+exp(b[12]))
  
  
  epsilonAB <- exp(b[13]) / (1+exp(b[13])+exp(b[14])+exp(b[15]))
  trans42 <- exp(b[14]) / (1+exp(b[13])+exp(b[14])+exp(b[15])) # epsilonB-epsilonAB
  trans43 <- exp(b[15]) / (1+exp(b[13])+exp(b[14])+exp(b[15])) # epsilonA-epsilonAB
  
  
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
        nrow = nstates)
      
    }
  }  

  
  
  Aprimary <- matrix(c(
    1-trans12-trans13-gammaAB,trans12,trans13,gammaAB,
    nuA,1-nuA-omegaAB-etaB,omegaAB,etaB,
    nuB,omegaBA,1-nuB-omegaBA-etaA,etaA,
    epsilonAB,trans42,trans43,1-trans42-trans43-epsilonAB), nrow = nstates, byrow = TRUE)
  Asecondary <- diag(1,nrow = nstates)
  A <- array(NA,c(nstates,nstates,N))
  for (j in primary) A[1:nstates,1:nstates,j] <- Aprimary
  for (k in secondary[-primary]) A[1:nstates,1:nstates,k] <- Asecondary
  
   

  #---------------------------------------------
  # calculate -log(lik)
  l <- 0
  for (i in 1:nh) # loop on sites
  {
    # Here we define occupancy as a linear function of patrolling effort described by a site-specific covariate 
    # therefore we extract occupancy estimates for each site i and update the intial state PI matrix at   
    oe <- garb[i] + 1 # first obs
    evennt <- data[,i] + 1 # non-det/det -> 1/2
    if(length(dim(PI))==2) 
    {
      ALPHA <- PI* B[oe,,i,1]
    }
    else 
    {
      ALPHA <- PI[,,i]* B[oe,,i,1] #two dimension as no spatio-temporal effects were retained from the smart_analysis0.1and0.1_DetectionEffect
    }
    for (j in 2:N)
    {
      ALPHA <- (ALPHA %*% A[1:nstates,1:nstates,j-1])*B[evennt[j],,i,j]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
}

