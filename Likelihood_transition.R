dev_2spocc_dyn3 <- function(b, data, eff, tempcov, garb, nh, primary, secondary, nstates, P_effects){
  
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
  psiB <- 1/(1+exp(-(b[1])))
  psiAb <- 1/(1+exp(-(b[2])))
  psiAB <- 1/(1+exp(-(b[2])))
  
  
  
  
  ResultEffect <- EffectsTransition(P_effects,b)
  epsilonAB<- ResultEffect[[1]]
  trans42<- ResultEffect[[2]]
  nuA<- ResultEffect[[3]]
  omegaAB<- ResultEffect[[4]]
  
  trans12<- ResultEffect[[5]]
  etaA<- ResultEffect[[6]]
  gammaAB<- ResultEffect[[7]]
  
  trans43<- ResultEffect[[8]]
  nuB<- ResultEffect[[9]]
  omegaBA<- ResultEffect[[10]]
  
  trans13<- ResultEffect[[11]]
  etaB<- ResultEffect[[12]]
  
  
#  gammaA = gammaAB + trans12
#  gammaB = gammaAB + trans13
#  epsilonB = epsilonAB + trans42
#  epsilonA = epsilonAB + trans43
  
  
  
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
  
  
  
  #---------------------------------------------
  # psiB <- 1/(1+exp(-(b[1]  + b[19]* coveff)))
  # psiAb <- 1/(1+exp(-(b[2] + b[20]* coveff)))
  # psiAB <- 1/(1+exp(-(b[3] + b[21]* coveff)))
  # 
  
  
  # transition prob (dynamic model)
  
  # initial states prob
  PI1 <- c(1-psiB,psiB) # 
  PI2 <- matrix(c(1-psiAb,psiAb,0,0,0,0,1-psiAB,psiAB),nrow=2,byrow=T) # 
  PI <- PI1 %*% PI2
  
  A <- array(NA,c(nstates,nstates,R, N))
  Asecondary <- diag(1,nrow = nstates)
  if(length(dim(omegaAB)>1))
  {  
    
    for(z in 1:R)
    {
      for (j in 1:N) {
        A[1:nstates,1:nstates,z,j] <- matrix(c(
          1-trans12[z,j]-trans13[z,j]-gammaAB[z,j],trans12[z,j],trans13[z,j],gammaAB[z,j],
          nuA[z,j],1-nuA[z,j]-omegaAB[z,j]-etaB[z,j],omegaAB[z,j],etaB[z,j],
          nuB[z,j],omegaBA[z,j],1-nuB[z,j]-omegaBA[z,j]-etaA[z,j],etaA[z,j],
          epsilonAB[z,j],trans42[z,j],trans43[z,j],1-trans42[z,j]-trans43[z,j]-epsilonAB[z,j]), nrow = nstates, byrow = TRUE)
      }
    }
    for (k in secondary[-primary]) A[1:nstates,1:nstates,,k] <- Asecondary
  } else {
    Aprimary <- matrix(c(
      1-trans12-trans13-gammaAB,trans12,trans13,gammaAB,
      nuA,1-nuA-omegaAB-etaB,omegaAB,etaB,
      nuB,omegaBA,1-nuB-omegaBA-etaA,etaA,
      epsilonAB,trans42,trans43,1-trans42-trans43-epsilonAB), nrow = nstates, byrow = TRUE)
    
    A <- array(NA,c(nstates,nstates,N))
    for (j in primary) A[1:nstates,1:nstates,j] <- Aprimary
    for (k in secondary[-primary]) A[1:nstates,1:nstates,k] <- Asecondary
  }
  
  
  
  
  # obs prob
  
  
  #---------------------------------------------
  # calculate -log(lik)
  l <- 0
  for (i in 1:nh) # loop on sites
  {
      # Here we define occupancy as a linear function of patrolling effort described by a site-specific covariate 
      # therefore we extract occupancy estimates for each site i and update the intial state PI matrix at   
      oe <- garb[i] + 1 # first obs
      evennt <- data[,i] + 1 # non-det/det -> 1/2
      ALPHA <- PI* B[oe,,i,1] 
      for (j in 2:N)
      {
        ifelse(length(dim(omegaAB)>1),
               ALPHA <- (ALPHA %*% A[1:nstates,1:nstates,i,j-1])*B[evennt[j],,i,j],
               ALPHA <- (ALPHA %*% A[1:nstates,1:nstates,j-1])*B[evennt[j],,i,j])
      }
      l <- l + logprot(sum(ALPHA))*eff[i]
    }
    l <- -l
    l
  }