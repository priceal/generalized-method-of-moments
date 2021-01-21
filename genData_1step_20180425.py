#########################################################
# 
# Script for running GMM 1 step method on simulated data
#
# This script is used to perform calculations reported in
# "Analyzing Dwell Times with the Generalized Method of Moments",
# by Sadie Piatt and Allen Price, submitted to PLOS ONE.
#
########################################################
# The script will test  GMM for one step process run all in a loop so don't
# have to store large amounts of simulated data
#
#  A  - >  B
#
# basic parameters are
# N[] = list of number of samples 
# tau1 = decay time 1
# order[] = list of orders of method
# trials = number of trials computed for each set of paramters (N,order)
#
# results will be in 2D arrays
# e.g., meanA[j,k] where 
# j = 0,1,...  index for N (number of samples)
# k = 0,1,...  index for order

# keep decay time 1 constant
tau1 = 10.0

# try for many values of N = number of samples
N = [5, 10, 20, 50, 100, 200, 500, 1000]

# try for different orders
order = [1,2,3,4]

# for each set of parameters (tau2,N,order) compute trails:
trials = 1000
tauRange = [[1.0],[10.0],[100.0],[1000.0]]

# initialize arrays to hold results
# these hold results for identity matrix
tempAi = zeros(trials)
meanAi = zeros([len(N),len(order)])
meandevAi = zeros([len(N),len(order)])
stdAi = zeros([len(N),len(order)])
seAi = zeros([len(N),len(order)])

# these hold results for diagonal jackknife
tempAd = zeros(trials)
meanAd = zeros([len(N),len(order)])
meandevAd = zeros([len(N),len(order)])
stdAd = zeros([len(N),len(order)])
seAd = zeros([len(N),len(order)])

# these hold results for full jackknife
tempAf = zeros(trials)
meanAf = zeros([len(N),len(order)])
meandevAf = zeros([len(N),len(order)])
stdAf = zeros([len(N),len(order)])
seAf = zeros([len(N),len(order)])

print "now calculating tau1 =",tau1
for j in range(len(N)): # loop through numbers of samples
    print "   ...N =", N[j]
    for k in range(len(order)): # loop through orders
        print "      ...order =", order[k]
        for n in range(trials): 
            times = sim.multi_poissonN([tau1],N[j])
            junki=gmm.gmmG(times,tauRange,n=order[k],diag=True,weight='iden',bc=True )
            tempAi[n]=junki  # this is smallest decay time
           
            junkd=gmm.gmmG(times,tauRange,n=order[k],diag=True,weight='jack',bc=True )
            tempAd[n]=junkd  # this is smallest decay time

            junkf=gmm.gmmG(times,tauRange,n=order[k],diag=False,weight='jack',bc=True )
            tempAf[n]=junkf  # this is smallest decay time

        meanAi[j,k],meandevAi[j,k],stdAi[j,k],seAi[j,k] = cu.stats(tempAi)
        meanAd[j,k],meandevAd[j,k],stdAd[j,k],seAd[j,k] = cu.stats(tempAd)
        meanAf[j,k],meandevAf[j,k],stdAf[j,k],seAf[j,k] = cu.stats(tempAf)
          




