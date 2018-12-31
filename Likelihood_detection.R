dev_2spocc_dyn <- function(b, data, eff,tempcov, roadcov, stationcov, streamcov, garb, nh, primary, secondary, nstates, P_effects){
  
  
  R<-nh
  # b = vector of parameters on real scale
  # data = site histories
  # eff = nb of sites with that particular history
  # garb = initial states
  # nh = nb of sites
  # nstates = nb of occupancy states
  
  #---------------------------------------------
  # apply multinomial logit and standard logit link
  
  
  
  psiB <- 1/(1+exp(-(b[1])))
  psiAb <- 1/(1+exp(-(b[2])))
  psiAB <- 1/(1+exp(-(b[3])))
  
  
  #---------------------------------------------
  # various quantities that will be useful later on
  K <- length(primary)
  J2 <- length(secondary)
  J <- J2/K
  N <- J * K
  #---------------------------------------------
  
  
  # 
  
  gammaAB <- exp(b[6]) / (1+exp(b[4])+exp(b[5])+exp(b[6]))
  trans12 <- exp(b[4]) / (1+exp(b[4])+exp(b[5])+exp(b[6])) # gammaA-gammaAB
  trans13 <- exp(b[5]) / (1+exp(b[4])+exp(b[5])+exp(b[6])) # gammaB-gammaAB
  
  
  nuA <- exp(b[7]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
  omegaAB <- exp(b[8]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
  etaB <- exp(b[9]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
  
  nuB <- exp(b[10]) / (1+exp(b[10])+exp(b[11])+exp(b[12]))
  omegaBA <- exp(b[11]) / (1+exp(b[10])+exp(b[11])+exp(b[12]))
  etaA <- exp(b[12]) / (1+exp(b[10])+exp(b[11])+exp(b[12]))
  
  
  epsilonAB <- exp(b[13]) / (1+exp(b[13])+exp(b[14])+exp(b[15]))
  trans42 <- exp(b[14]) / (1+exp(b[13])+exp(b[14])+exp(b[15])) # epsilonB-epsilonAB
  trans43 <- exp(b[15]) / (1+exp(b[13])+exp(b[14])+exp(b[15])) # epsilonA-epsilonAB
  
  # 
  result_effect <- EffectsDetection(P_effects,b) 
  pB <- result_effect [[1]]
  pAb <- result_effect [[2]]
  pAB <- result_effect [[3]]
  B <-  result_effect [[4]]
  
  #---------------------------------------------
  # psiB <- 1/(1+exp(-(b[1]  + b[19]* coveff)))
  # psiAb <- 1/(1+exp(-(b[2] + b[20]* coveff)))
  # psiAB <- 1/(1+exp(-(b[3] + b[21]* coveff)))
  # 
  # # initial states prob
  PI1 <- array(data=NA, dim=c(1,2))
  PI2 <- array(data=NA, dim=c(2,4))
  PI <- array(data=NA, dim=c(1,4))
  
  PI1 <- c(1-psiB,psiB)
  PI2 <- matrix(c(1-psiAb,psiAb,0,0,0,0,1-psiAB,psiAB),nrow=2,byrow=T) # 
  PI <- PI1 %*% PI2
  
  
  # transition prob (dynamic model)
  
  # initial states prob
  #PI1 <- c(1-psiB,psiB) # 
  #PI2 <- matrix(c(1-psiAb,psiAb,0,0,0,0,1-psiAB,psiAB),nrow=2,byrow=T) # 
  #PI <- PI1 %*% PI2
  
  
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
    if(length(dim(B))==4) 
    {
      Bi <- B[oe,,i,1]
    }
    else if(length(dim(B))==3) 
    {
      if(dim(B)[3]==R)
        Bi <- B[oe,,i]
      else  
        Bi <- B[oe,,1]
    }
    if(length(dim(B))==2) 
    {
      Bi <- B[oe,]
    }
    ALPHA <- PI*Bi
    for (j in 2:N)
    {
      if(length(dim(B))==4) 
      {
        Blike <- B[evennt[j],,i,j]
      }
      else if(length(dim(B))==3) 
      {
        if(dim(B)[3]==R)
          Blike <- B[evennt[j],,i]
        if(dim(B)[3]==N) 
          Blike <- B[evennt[j],,j]
        if(dim(B)[3]==(K-1)) 
        {
          y<-(j %/% J)+1
          Blike <- B[evennt[j],,y]  
        }
      }
      if(length(dim(B))==2) 
      {
        Blike <- B[evennt[j],]
      }
      ALPHA <- (ALPHA %*% A[1:nstates,1:nstates,j-1])*Blike
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
}
