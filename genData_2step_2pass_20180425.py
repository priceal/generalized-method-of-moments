#########################################################
# 
# Script for running GMM 2 step/2 pass method on simulated data
#
# This script is used to perform calculations reported in
# "Analyzing Dwell Times with the Generalized Method of Moments",
# by Sadie Piatt and Allen Price, submitted to PLOS ONE.
#
########################################################
# The script will test  GMM for 2 step process run all in a loop so don't
# have to store large amounts of simulated data.  This took about 5 hours
# to run on 3 GHz dual core machine.
#
# Runs first iteration as second order diagonal jackknife
# Tries 3rd order and 4th order with full interpolated matrix
#
# basic parameters are
# N[] = list of number of samples (molecules)
# tau1 = decay time 1
# tau2[] = list of decay times 2 
#
# results will be in 3D arrays
# e.g., meanA2[i,j] where 
# i = 0,1,...  index for tau2
# j = 0,1,...  index for N (number of samples)

# keep decay time 1 constant, and vary time 2
tau1 = 10.0
tau2 = [10.0, 20.0, 50.0, 100.0 ]

# try for many values of N = number of samples
N = [5, 10, 20, 50, 100, 200, 500, 1000]

# for each set of parameters (tau2,N,order) compute trA1fls:
trials = 1000
tauRange = [[1.0,10.0],[1.0,100.0],[1.0,1000.0],[10.0,100.0],[10.0,1000.0],[100.0,1000.0]]

# initialize arrays to hold results
# these hold results for tauA, 1 pass
tempA1 =zeros(trials)
meanA1 =zeros([len(tau2),len(N)])
meandevA1 =zeros([len(tau2),len(N)])
stdA1 =zeros([len(tau2),len(N)])
seA1 =zeros([len(tau2),len(N)])
# these hold results for tauB, 1 pass
tempB1 =zeros(trials)
meanB1 =zeros([len(tau2),len(N)])
meandevB1 =zeros([len(tau2),len(N)])
stdB1 =zeros([len(tau2),len(N)])
seB1 =zeros([len(tau2),len(N)])

# these hold results for tauA, 2nd pass, 3rd order
tempA2 =zeros(trials)
meanA2 =zeros([len(tau2),len(N)])
meandevA2 =zeros([len(tau2),len(N)])
stdA2 =zeros([len(tau2),len(N)])
seA2 =zeros([len(tau2),len(N)])
# these hold results for tauB, 2nd pass, 3rd order
tempB2 =zeros(trials)
meanB2 =zeros([len(tau2),len(N)])
meandevB2 =zeros([len(tau2),len(N)])
stdB2 =zeros([len(tau2),len(N)])
seB2 =zeros([len(tau2),len(N)])

# these hold results for tauA, 2nd pass, 4th order
tempA3 =zeros(trials)
meanA3 =zeros([len(tau2),len(N)])
meandevA3 =zeros([len(tau2),len(N)])
stdA3 =zeros([len(tau2),len(N)])
seA3 =zeros([len(tau2),len(N)])
# these hold results for tauB, 2nd pass, 4th order
tempB3 =zeros(trials)
meanB3 =zeros([len(tau2),len(N)])
meandevB3 =zeros([len(tau2),len(N)])
stdB3 =zeros([len(tau2),len(N)])
seB3 =zeros([len(tau2),len(N)])

for i in range(len(tau2)): # loop through tau2 values
    print "now calculating tau2 =",tau2[i]
    for j in range(len(N)): # loop through numbers of samples
        print "   ...N =", N[j]
        for n in range(trials): 
            times = sim.multi_poissonN([tau1,tau2[i]],N[j])

            junk1=gmm.gmmG(times,tauRange,n=2,diag=True,weight='jack',bc=True )
            tempA1[n]=junk1[0]  # this is smallest decay time
            tempB1[n]=junk1[1]  # this is largest decay time

            junk2=gmm.gmm(times,[tempA1[n],tempB1[n]],n=3,diag=False,weight='int',bc=True )
            tempA2[n]=junk2[0]  # this is smallest decay time
            tempB2[n]=junk2[1]  # this is largest decay time

            junk3=gmm.gmm(times,[tempA1[n],tempB1[n]],n=4,diag=False,weight='int',bc=True )
            tempA3[n]=junk3[0]  # this is smallest decay time
            tempB3[n]=junk3[1]  # this is largest decay time

        meanA1[i,j],meandevA1[i,j],stdA1[i,j],seA1[i,j] = cu.stats(tempA1)
        meanB1[i,j],meandevB1[i,j],stdB1[i,j],seB1[i,j] = cu.stats(tempB1)

        meanA2[i,j],meandevA2[i,j],stdA2[i,j],seA2[i,j] = cu.stats(tempA2)
        meanB2[i,j],meandevB2[i,j],stdB2[i,j],seB2[i,j] = cu.stats(tempB2)

        meanA3[i,j],meandevA3[i,j],stdA3[i,j],seA3[i,j] = cu.stats(tempA3)
        meanB3[i,j],meandevB3[i,j],stdB3[i,j],seB3[i,j] = cu.stats(tempB3)
 

