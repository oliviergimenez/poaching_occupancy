

# Lucile Marescot 12/30/2018

# 1 Simulation study of two-species occupancy data focusing on the occupancy design   
# In both simulation study we used the same transition parameter as the one obatined from the best model of the SMART data analysis
# here we vary occupancy design (nb of sites, nb of primary occasions, nb of secondary occasions) 
# (35 sites with 12 repeated survey during 6 sampling sessions)


rm(list=ls()) # remove all objects/memory
gc()




setwd("C:/Users/Marescot/Documents/GitHub/sdp_occ_wwf/Last files")



source('Likelihood_CSTeffects.R')




# a) Simulate the data sets 
library(HMM)


# transition parameters (see Table 8.7 page 245 in MacKenzie et al. 2006)
epsilonAB = 0.0001 # 
epsilonA = 0.20 # 
epsilonB = 0.79 # 
nuA = 0.14
nuB = 0.41
gammaAB = 0.001 
gammaA = 0.01 
gammaB = 0.28 
etaA = 0.01
etaB = 0.05
omegaAB = 0.12 
omegaBA = 0.43

n.obs <- 4 # nb events Unoccupied, occupied by poachers only, by prey only and by both
n.states <- 4 # nb states



# The first simulation study consists in fixing the parameters of occupancy and detection and varying the occupancy design

sites = c(35,100,250, 500) # range of number of sites 



primary.occasion <- c(3,5,10) # number of primary occasions
sec <- c(3, 6, 12) # number of secondary occasions
MCiter<-1000 # number of iterations

tab<-array(data=0, dim=c(18, length(sites), length(primary.occasion),length(sec), MCiter)) # table with all possible combinations of occupancy design
dfplus<-NULL # data frame to save the results of simulations


for(it in 1:MCiter)
{ 
  inc1 <-0 # increment for the loop on site
  for(v in sites)
  {
    inc1 <- inc1 + 1 
    inc2 <-0 # increment for the loop on primary occasions
    for(p in primary.occasion)
    {  
      inc2 <- inc2 + 1  
      
      inc3 <-0 # increment for the loop on secondary occasions
      for(s in sec)
      {
        inc3 <- inc3 + 1  
        R = v
        n.secondary <- s # nb secondary occasions
        n.primary <- p   # nb primary occasions    
        
        #source('Multispecies_occupancy_dynamic simul_fitting.R')
        # occupancy in Waddle parameterization
        psiAB = 0.06 # probability that species A is present, given that species B is present
        psiAb = 0.06 # probability that species A is present, given that species B is absent
        psiB = 0.46 # probability that species B is present, regardless species A
        
        
        # SMART
        pB = 0.10 # probability of detecting B, Regardless of occurence of A
        pAb = 0.17 # probability of detecting only species A, given that B is absent
        pAB = 0.50 # probability of detecting only species A, given that both are present
        
        
          dyn <- matrix(c(
            1-gammaA-gammaB+gammaAB,gammaA-gammaAB,gammaB-gammaAB,gammaAB,
            nuA,1-nuA-omegaAB-etaB,omegaAB,etaB,
            nuB,omegaBA,1-nuB-omegaBA-etaA,etaA,
            epsilonAB,epsilonB-epsilonAB,epsilonA-epsilonAB,1-epsilonA-epsilonB+epsilonAB), nrow = n.states, byrow = TRUE)
          
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
          
          
          # b) model fitting  
          
          #Run optimisation
          binit = runif(18) # generate initial values for parameters
          res = optim(binit,dev_2spocc_dyn,NULL,hessian=FALSE,data,eff,garb,R,primary,secondary,n.states,method="BFGS",control=list(trace=0, REPORT=1))
          
          
          
          
          
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
          
          resdf=data.frame(param=c('psiB', 
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
          
        
          
          
          
          tab[,inc1,inc2,inc3,it]<- (resdf[,3] -  resdf[,2])
         }
      }
  }
} 



# c) calculate the mean bias and minimum square error 

library(plyr)  

df1<-apply(tab,5, mean)
df1<-adply(tab, c(2,3,4))

df1[,1]<-sites[df1[,1]]
df1[,2]<-primary.occasion[df1[,2]]
df1[,3]<-sec[df1[,3]]

param=c('psiB', 
        'psiAb', 
        'psiAB',
        'gammaA', 
        'gammaB', 
        'gammaAB', 
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
        'pAB')

df1<-t(as.matrix(df1))

#colnames(df1)[c(1,2,3)]<- c("sites", "primary.occasion", "secondary.occasion")
#colnames(df1)[4: ncol(df1)]<- param




write.table(df1, "biasSmartCellNumber.txt")

save(df1, tab, file="ResultsSimu_SebCluster_12SEP.RData")
