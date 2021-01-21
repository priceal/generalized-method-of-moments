#########################################################
# 
# Script for running GMM 2 step/1 pass method on simulated data
#
# This script is used to perform calculations reported in
# "Analyzing Dwell Times with the Generalized Method of Moments",
# by Sadie Piatt and Allen Price, submitted to PLOS ONE.
#
########################################################
# The script will test  GMM for 2 step process run all in a loop so don't
# have to store large amounts of simulated data.  This took about 34 hours
# to run on 3 GHz dual core machine.
#
# basic parameters are
# N[] = list of number of samples (molecules)
# tau1 = decay time 1
# tau2[] = list of decay times 2 
# trials = number of trials computed for each set of paramters (tau2,N,order)
# order[] = list of orders of method
#
# results will be in 3D arrays
# e.g., meanA[i,j,k] where 
# i = 0,1,...  index for tau2
# j = 0,1,...  index for N (number of samples)
# k = 0,1,...  index for order

# keep decay time 1 constant, and vary time 2
tau1 = 10.0
tau2 = [10.0, 20.0, 50.0, 100.0 ]

# try for many values of N = number of samples
N = [5, 10, 20, 50, 100, 200, 500, 1000]


# try for different orders
order = [ 2,3,4 ]

# for each set of parameters (tau2,N,order) compute trA1fls:
trials = 1000
tauRange = [[1.0,10.0],[1.0,100.0],[1.0,1000.0],[10.0,100.0],[10.0,1000.0],[100.0,1000.0]]


# initialize arrays to hold results
# these hold results for tauA, diagonal jackknife
tempAd = zeros(trials)
meanAd = zeros([len(tau2),len(N),len(order)])
meandevAd = zeros([len(tau2),len(N),len(order)])
stdAd = zeros([len(tau2),len(N),len(order)])
seAd = zeros([len(tau2),len(N),len(order)])
# these hold results for tauB, diagonal jackknife
tempBd = zeros(trials)
meanBd = zeros([len(tau2),len(N),len(order)])
meandevBd = zeros([len(tau2),len(N),len(order)])
stdBd = zeros([len(tau2),len(N),len(order)])
seBd = zeros([len(tau2),len(N),len(order)])

# these hold results for tauA, full jackknife
tempAf = zeros(trials)
meanAf = zeros([len(tau2),len(N),len(order)])
meandevAf = zeros([len(tau2),len(N),len(order)])
stdAf = zeros([len(tau2),len(N),len(order)])
seAf = zeros([len(tau2),len(N),len(order)])
# these hold results for tauB, full jackknife
tempBf = zeros(trials)
meanBf = zeros([len(tau2),len(N),len(order)])
meandevBf = zeros([len(tau2),len(N),len(order)])
stdBf = zeros([len(tau2),len(N),len(order)])
seBf = zeros([len(tau2),len(N),len(order)])

# these hold results for tauA, identity weighting
tempAi = zeros(trials)
meanAi = zeros([len(tau2),len(N),len(order)])
meandevAi = zeros([len(tau2),len(N),len(order)])
stdAi = zeros([len(tau2),len(N),len(order)])
seAi = zeros([len(tau2),len(N),len(order)])
# these hold results for tauB, identity weighting
tempBi = zeros(trials)
meanBi = zeros([len(tau2),len(N),len(order)])
meandevBi = zeros([len(tau2),len(N),len(order)])
stdBi = zeros([len(tau2),len(N),len(order)])
seBi = zeros([len(tau2),len(N),len(order)])

for i in range(len(tau2)): # loop through tau2 values
    print "now calculating tau2 =",tau2[i]
    for j in range(len(N)): # loop through numbers of samples
        print "   ...N =", N[j]
        for k in range(len(order)):
            print "      ...order =", order[k]
            for n in range(trials): 
                times = sim.multi_poissonN([tau1,tau2[i]],N[j])

                junkd=gmm.gmmG(times,tauRange,n=order[k],diag=True,weight='jack',bc=True )
                tempAd[n]=junkd[0]  # this is smallest decay time
                tempBd[n]=junkd[1]  # this is largest decay time

                junkf=gmm.gmmG(times,tauRange,n=order[k],diag=False,weight='jack',bc=True )
                tempAf[n]=junkf[0]  # this is smallest decay time
                tempBf[n]=junkf[1]  # this is largest decay time

                junki=gmm.gmmG(times,tauRange,n=order[k],diag=True,weight='iden',bc=True )
                tempAi[n]=junki[0]  # this is smallest decay time
                tempBi[n]=junki[1]  # this is largest decay time

            meanAd[i,j,k],meandevAd[i,j,k],stdAd[i,j,k],seAd[i,j,k] = cu.stats(tempAd)
            meanBd[i,j,k],meandevBd[i,j,k],stdBd[i,j,k],seBd[i,j,k] = cu.stats(tempBd)

            meanAf[i,j,k],meandevAf[i,j,k],stdAf[i,j,k],seAf[i,j,k] = cu.stats(tempAf)
            meanBf[i,j,k],meandevBf[i,j,k],stdBf[i,j,k],seBf[i,j,k] = cu.stats(tempBf)

            meanAi[i,j,k],meandevAi[i,j,k],stdAi[i,j,k],seAi[i,j,k] = cu.stats(tempAi)
            meanBi[i,j,k],meandevBi[i,j,k],stdBi[i,j,k],seBi[i,j,k] = cu.stats(tempBi)






