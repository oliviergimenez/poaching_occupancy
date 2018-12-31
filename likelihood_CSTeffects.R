#the Likelihood Function of a two species dynamic occupancy model with constant parameters


logprot <- function(v){
  # avoid explosion of log(0)
  eps <- 2.2204e-016
  u <- log(eps) * (1+vector(length=length(v)))
  index <- (v>eps)
  u[index] <- log(v[index])
  u
}


#Write down likelihood of a two-species occupancy model: 

dev_2spocc_dyn <- function(b,data,eff,garb,nh,primary,secondary,nstates){
  # b = vector of parameters on real scale
  # data = site histories
  # eff = nb of sites with that particular history
  # garb = initial states
  # nh = nb of sites
  # nstates = nb of occupancy states
  
  #---------------------------------------------
  # apply multinomial logit and standard logit link
  psiB <- 1/(1+exp(-b[1]))
  psiAb <- 1/(1+exp(-b[2]))
  psiAB <- 1/(1+exp(-b[3]))
  
  trans12 <- exp(b[4]) / (1+exp(b[4])+exp(b[5])+exp(b[6])) # gammaA-gammaAB
  trans13 <- exp(b[5]) / (1+exp(b[4])+exp(b[5])+exp(b[6])) # gammaB-gammaAB
  gammaAB <- exp(b[6]) / (1+exp(b[4])+exp(b[5])+exp(b[6]))
  
  nuA <- exp(b[7]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
  omegaAB <- exp(b[8]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
  etaB <- exp(b[9]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
  
  nuB <- exp(b[10]) / (1+exp(b[10])+exp(b[11])+exp(b[12]))
  omegaBA <- exp(b[11]) / (1+exp(b[10])+exp(b[11])+exp(b[12]))
  etaA <- exp(b[12]) / (1+exp(b[10])+exp(b[11])+exp(b[12]))
  
  epsilonAB <- exp(b[13]) / (1+exp(b[13])+exp(b[14])+exp(b[15]))
  trans42 <- exp(b[14]) / (1+exp(b[13])+exp(b[14])+exp(b[15])) # epsilonB-epsilonAB
  trans43 <- exp(b[15]) / (1+exp(b[13])+exp(b[14])+exp(b[15])) # epsilonA-epsilonAB
  
  pB <- 1/(1+exp(-b[16]))
  pAb <- 1/(1+exp(-b[17]))
  pAB <- 1 / (1+exp(-b[18]))
  
  #---------------------------------------------
  
  #---------------------------------------------
  # various quantities that will be useful later on
  K <- length(primary)
  J2 <- length(secondary)
  J <- J2/K
  N <- J * K
  #---------------------------------------------
  
  # initial states prob
  PI1 <- c(1-psiB,psiB) # 
  PI2 <- matrix(c(1-psiAb,psiAb,0,0,0,0,1-psiAB,psiAB),nrow=2,byrow=T) # 
  PI <- PI1 %*% PI2
  
  # transition prob (dyanmic model)
  Aprimary <- matrix(c(
    1-trans12-trans13-gammaAB,trans12,trans13,gammaAB,
    nuA,1-nuA-omegaAB-etaB,omegaAB,etaB,
    nuB,omegaBA,1-nuB-omegaBA-etaA,etaA,
    epsilonAB,trans42,trans43,1-trans42-trans43-epsilonAB), nrow = nstates, byrow = TRUE)
  Asecondary <- diag(1,nrow = nstates)
  A <- array(NA,c(nstates,nstates,N))
  for (j in primary) A[1:nstates,1:nstates,j] <- Aprimary
  for (k in secondary[-primary]) A[1:nstates,1:nstates,k] <- Asecondary
  
  # obs prob
  B <- matrix(c(
    1,0,0,0,
    1-pAb,pAb,0,0,
    1-pB,0,pB,0,
    1-pAB,0,0,pAB),
    nrow = nstates)
  
  #---------------------------------------------
  # calculate -log(lik)
  l <- 0
  for (i in 1:nh) # loop on sites
  {
    oe <- garb[i] + 1 # first obs
    evennt <- data[,i] + 1 # non-det/det -> 1/2
    ALPHA <- PI*B[oe,]
    for (j in 2:N){ # loop on time
      ALPHA <- (ALPHA %*% A[1:nstates,1:nstates,j-1])*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
}

