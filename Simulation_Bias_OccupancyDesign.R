

# Lucile Marescot 12/30/2018

# 1 Simulation study of two-species occupancy data focusing on the occupancy design   
# In both simulation study we used the same transition parameter as the one obatined from the best model of the SMART data analysis
# here we vary occupancy design (nb of sites, nb of primary occasions, nb of secondary occasions) 
# (35 sites with 12 repeated survey during 6 sampling sessions)


rm(list=ls()) # remove all objects/memory
gc()




setwd("C:/Users/Marescot/Documents/GitHub/sdp_occ_wwf/Last files")




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
psiAB = 0.06 # probability that species A is present, given that species B is present
psiAb = 0.06 # probability that species A is present, given that species B is absent
psiB = 0.46 # probability that species B is present, regardless species A


# SMART
pB = 0.10 # probability of detecting B, Regardless of occurence of A
pAb = 0.17 # probability of detecting only species A, given that B is absent
pAB = 0.50 # probability of detecting only species A, given that both are present




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
        
        source('Multispecies_occupancy_dynamic_simul_fitting.R')
        
          tab[,inc1,inc2,inc3,it]<- (res[,3] -  res[,2])
      }
      }
  }
} 



# c) calculate the mean bias and minimum square error 

library(plyr)  

dfmean<-apply(tab,c(1,2,3,4), mean)
df1<-adply(dfmean, c(2,3,4))

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


colnames(df1)[c(1,2,3)]<- c("sites", "primary.occasion", "secondary.occasion")
colnames(df1)[4: ncol(df1)]<- param




write.table(df1, "biasSmartCellNumber.txt")

save(df1, tab, file="ResultsSimu_SebCluster_12SEP.RData")
