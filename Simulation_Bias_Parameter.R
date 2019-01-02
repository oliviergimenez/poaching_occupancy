# Lucile Marescot 12/30/2018

# 1 Simulation study of two-species occupancy data focusing on the occupancy design   
# In both simulation study we used the same transition parameter as the one obatined from the best model of the SMART data analysis
# here we vary detection and occupancy probabilities and fixed the study design




rm(list=ls()) # remove all objects/memory
gc()


setwd("C:/Users/Marescot/Documents/GitHub/sdp_occ_wwf/Simulation")


# occupancy in Waddle parameterization
psiABvec = seq(0.1, 0.5, 0.15) # probability that species A is present, given that species B is present
psiAbvec = seq(0.1, 0.5, 0.15) # probability that species A is present, given that species B is absent
psiBvec = seq(0.1, 0.5, 0.15) # probability that species B is present, regardless species A

# SMART
pBvec = seq(0.1, 0.5, 0.15) # probability of detecting B, Regardless of occurence of A
pAbvec = seq(0.1, 0.5, 0.15) # probability of detecting only species A, given that B is absent
pABvec = seq(0.1, 0.5, 0.15) # probability of detecting only species A, given that both are present


epsilonAB = 0.1 # 
epsilonA = 0.4 # 
epsilonB = 0.3 # 
nuA = 0.1
nuB = 0.1
gammaAB = 0.1 
gammaA = 0.2 
gammaB = 0.5 
etaA = 0.6
etaB = 0.4
omegaAB = 0.1 
omegaBA = 0.2

# occupancy design fixed
R = 35 # number of sites 
n.states <- 4 # nb states
n.obs <- 4 # nb events
n.primary <- 4 # nb primary occasions
n.secondary <- 12 # nb secondary occasions




dfplus<-NULL

MCiter<-2


tab<-array(data=0, dim=c(18, length(psiAbvec), length(psiABvec),length(psiBvec),length(pBvec), MCiter))
dfmean<-array(data=0, dim=c(18, length(psiAbvec), length(psiABvec),length(psiBvec),length(pBvec)))
inc1 <-0
for(Aocc in psiABvec)
{
  psiAB <- Aocc   
  inc1 <- inc1 + 1 # increment on the loop of the occupancy vector of species A given presence of B
  inc2 <-0
  for(Abocc in psiAbvec)
  {
    psiAb <- Abocc # increment on the loop of the occupancy vector of species B given presence of A
    inc2 <- inc2 + 1  

    inc3 <-0 
    for(Bocc in psiBvec)
    {
      psiB <- Bocc # increment on the loop of the occupancy vector of species B 
      inc3 <- inc3 + 1  
      
      inc4 <-0
      for(pBocc in pBvec)
      {
        pB <- pAB <- pAb <- pBocc
        inc4 <- inc4 + 1  # increment on the loop of the detection probability
        incmean<-rep(0, 18) 
      for(it in 1:MCiter)
      { 
        source('Multispecies_occupancy_dynamic_simul_fitting.R')
        
        tab[,inc1,inc2,inc3,inc4, it]<- (res[,3] - res[,2])
        incmean<- incmean + (res[,3] - res[,2])
        
      }
      dfmean[,inc1,inc2,inc3,inc4] <- incmean / MCiter
      
      # for(it in 1:MCiter)
      # { 
      #   incsum<-tab[,inc1,inc2,inc3,it]-dfmean[,inc1,inc2,inc3]
      #   MSE[,inc1,inc2,inc3]<- MSE[,inc1,inc2,inc3] + incsum 
      # }
      #   MSE[,inc1,inc2,inc3]<- MSE[,inc1,inc2,inc3]/MCiter
      }
    }
  }
}


library(plyr)  

df<-adply(dfmean, c(2,3,4,5))

df[,1]<-psiABvec[df[,1]]
df[,2]<-psiAbvec[df[,2]]
df[,3]<-psiBvec[df[,3]]
df[,4]<-pBvec[df[,4]]




param<-c(psiB, psiAb, 
               psiAB,
               gammaA, 
               gammaB, 
               gammaAB, 
               nuA,
               omegaAB, 
               etaB,
               nuB,
               omegaBA, 
               etaA,
               epsilonAB, 
               epsilonB,
               epsilonA,
               pB,
               pAb, 
               pAB)
names(param)=c('psiB', 
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

colnames(df)[c(1,2,3,4)]<- c("psiAB", "psiAb", "psiB", "pB")
#colnames(df)[4: ncol(df)]<- names(param)





write.table(df, "biasSmartPARAM.txt")
write.table(dfMSE, "MSESmart.txt")



save(tab, dfmean, file = "SimulationBiasSmartPARAMDetection.RData")

save()
install.packages("gridExtra")
require(gridExtra) # also loads grid
# # 
# # 
layout()
# # 
colnames(df)[3]<-"occupancyA"
p<-list()
# # 
library(lattice)
inc<-0
for(i in unique(df[,"pB"]))
{
  inc<-inc+1
  select<-df[df[,"pB"]==i,] 
  p[[inc]]<-levelplot(psiAb ~ psiAbvec * psiABvec | occupancyA , data=select, at=seq(-0.3,0.3,0.01), main= "", xlab=expression(psi["P/A"]), ylab=expression(psi["P/noA"]),  zlab="",  colorkey = TRUE , region = TRUE, contour = FALSE, pretty=TRUE, col.regions= (grey.colors(100)))  
  #if(inc >3) p[[i]]<-levelplot(df[,i] ~ psiAbvec * psiABvec | occupancyA , data=df, at=seq(-0.1,0.1,0.02), main= paste (i), xlab=expression(psi["P/A"]), ylab=expression(psi["P/noA"]),  zlab="",  colorkey = TRUE , region = TRUE, contour = FALSE, pretty=TRUE, col.regions= (cm.colors(100)))  
  # print(p, split = c(1, 1, 2, 2), more = TRUE)
}
im1<-grid.arrange(p[[1]], p[[2]], p[[3]], nrow=3, ncol=1)


#save(df,  im1, im2, im3, im4, im5, file = "WWFSimulationBiasSmartV2.RData")
#save(df, tab, file = "WWFSimulationBiasSmartPARAdetection.RData")
#save(tab, file = "WWFSimulationBiasSmartPARA.RData")



