
res<-NULL

library(HMM)

 
dyn <- matrix(c(
  	1-gammaA-gammaB+gammaAB,gammaA-gammaAB,gammaB-gammaAB,gammaAB,
		nuA,1-nuA-omegaAB-etaB,omegaAB,etaB,
  	nuB,omegaBA,1-nuB-omegaBA-etaA,etaA,
		epsilonAB,epsilonB-epsilonAB,epsilonA-epsilonAB,1-epsilonA-epsilonB+epsilonAB), nrow = n.states, byrow = TRUE)

if(min(dyn) > 0)
{ 
hmm <- initHMM(c("U","OA","OB","OAB"), c("1","2","3","4"), 
c((1-psiB)*(1-psiAb),(1-psiB)*psiAB, psiB*(1-psiAB),psiB*psiAB), 
dyn,
matrix(c(
		1,0,0,0,
		0,1,0,0,
		0,0,1,0,
		0,0,0,1), nrow = n.states, byrow = TRUE))



simHMM(hmm,n.primary)

  
#Simulate a dataset of nb.ind individuals

dat = NULL
dat_st = NULL
i = 0
while (i < R){
	temp = simHMM(hmm,n.primary)
	current_ind = temp$observation # this might be a row of non-detections
	current_ind_st = temp$states
	dat = rbind(dat, current_ind)
	dat_st = rbind(dat_st, current_ind_st)
	i = i+1
}
events = matrix(as.numeric(dat),byrow=F,ncol=n.primary)
data = events - 1
states = dat_st
head(data)
head(dat_st)


#Second, add on top of the states the observation process. Note that states (dat_st) and observations (data) coincide for now because detection is perfect. Let us use the observations as they are in numeric format. As a reminder:

# 0 is for U  
# 1 is for OA  
# 2 is for OB  
# 3 is for OAB  


obs = array(NA,dim=c(R,n.primary,n.secondary))
for (i in 1:nrow(data)){
  current_site = data[i,]
  for (j in 1:n.primary){
    if (current_site[j]==0) obs[i,j,1:n.secondary] = 0 # if site unoccupied, species not seen whatever secondary occasion
    if (current_site[j]==1){    # if site occupied by species A only, species A seen with proba pAb
      for (k in 1:n.secondary){
        obs[i,j,k] = rbinom(1,1,pAb)
      }
    }
    if (current_site[j]==2){    # if site occupied by species B only, species B seen with proba pB
      for (k in 1:n.secondary){
        obs[i,j,k] = 2 * rbinom(1,1,pB)
      }
    }
    if (current_site[j]==3){    # if site occupied by both species, both species are seen with proba pAB
      for (k in 1:n.secondary){
        obs[i,j,k] = 3 * rbinom(1,1,pAB)
      }
    }    
  }
}


#Now reshape array in matrix with secondary occ within primary occ:

obs2 = matrix(NA,nrow=R,ncol=n.primary*n.secondary)
count = 1
for (i in 1:n.primary){
  for (j in 1:n.secondary){
    obs2[,count] <- obs[,i,j]
    count = count+1
  }
}


## Model fitting


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
      	l <- l + log(sum(ALPHA))*eff[i]
   }
l <- -l
l
}


#Various quantities

J = n.secondary
K = n.primary
R = dim(obs2)[1] # nb sites
eff = rep(1,R) # nb sites with this particular history
garb = obs2[,1] # initial states
data = t(obs2) # transpose

# primary and secondary occasions
primary = seq(J,J*K,by=J)
secondary = 1:(J*K)

#Run optimisation
#binit = runif(18)


binit = runif(18) # generate initial values for parameters
res = optim(binit,dev_2spocc_dyn,NULL,hessian=FALSE,data,eff,garb,R,primary,secondary,n.states,method="BFGS",control=list(trace=0, REPORT=1))

#dev <- res$value

}



#Get estimates and back-transform

b = res$par
mle_psiB <- 1/(1+exp(-b[1]))
mle_psiAb <- 1/(1+exp(-b[2]))
mle_psiAB <- 1/(1+exp(-b[3]))
mle_gammaAB <- exp(b[6]) / (1+exp(b[4])+exp(b[5])+exp(b[6]))
trans12 <- exp(b[4]) / (1+exp(b[4])+exp(b[5])+exp(b[6])) # gammaA-gammaAB
trans13 <- exp(b[5]) / (1+exp(b[4])+exp(b[5])+exp(b[6])) # gammaB-gammaAB
mle_gammaA = gammaAB + trans12
mle_gammaB = gammaAB + trans13
mle_nuA <- exp(b[7]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
mle_omegaAB <- exp(b[8]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
mle_etaB <- exp(b[9]) / (1+exp(b[7])+exp(b[8])+exp(b[9]))
mle_nuB <- exp(b[10]) / (1+exp(b[10])+exp(b[11])+exp(b[12]))
mle_omegaBA <- exp(b[11]) / (1+exp(b[10])+exp(b[11])+exp(b[12]))
mle_etaA <- exp(b[12]) / (1+exp(b[10])+exp(b[11])+exp(b[12]))
mle_epsilonAB <- exp(b[13]) / (1+exp(b[13])+exp(b[14])+exp(b[15]))
trans42 <- exp(b[14]) / (1+exp(b[13])+exp(b[14])+exp(b[15])) # epsilonB-epsilonAB
trans43 <- exp(b[15]) / (1+exp(b[13])+exp(b[14])+exp(b[15])) # epsilonA-epsilonAB
mle_epsilonB = epsilonAB + trans42
mle_epsilonA = epsilonAB + trans43
mle_pB <- 1/(1+exp(-b[16]))
mle_pAb <- 1/(1+exp(-b[17]))
mle_pAB <- 1 / (1+exp(-b[18]))


#Compare true parameter to their estimates:

res=data.frame(param=c('psiB', 
                       'psiAb', 
                       'psiAB',
                       'gammaAB', 
                       'gammaA', 
                       'gammaB', 
                       'nuA',
                       'omegaAB', 
                       'etaB',
                       'nuB',
                       'omegaBA', 
                       'etaA',
                       'epsilonAB', 
                       'epsilonB',
                       'epsilonA',
                       'pB',
                       'pAb', 
                       'pAB'),truth=c(psiB, 
                                      psiAb, 
                                      psiAB ,
                                      gammaAB, 
                                      gammaA, 
                                      gammaB, 
                                      nuA ,
                                      omegaAB, 
                                      etaB ,
                                      nuB,
                                      omegaBA, 
                                      etaA ,
                                      epsilonAB, 
                                      epsilonB ,
                                      epsilonA ,
                                      pB,
                                      pAb, 
                                      pAB),estimates=c(mle_psiB, 
                                                       mle_psiAb, 
                                                       mle_psiAB ,
                                                       mle_gammaAB, 
                                                       mle_gammaA, 
                                                       mle_gammaB, 
                                                       mle_nuA ,
                                                       mle_omegaAB, 
                                                       mle_etaB ,
                                                       mle_nuB,
                                                       mle_omegaBA, 
                                                       mle_etaA ,
                                                       mle_epsilonAB, 
                                                       mle_epsilonB ,
                                                       mle_epsilonA ,
                                                       mle_pB ,
                                                       mle_pAb, 
                                                       mle_pAB))


