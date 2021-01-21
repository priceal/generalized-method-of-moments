#########################################################
# 
# Script for running GMM 3 step/1 pass method on simulated data
#
# This script is used to perform calculations reported in
# "Analyzing Dwell Times with the Generalized Method of Moments",
# by Sadie Piatt and Allen Price, submitted to PLOS ONE.
#
########################################################
# The script will test  GMM for 3 step process run all in a loop so don't
# have to store large amounts of simulated data.  
#
# basic parameters are
# N[] = list of number of samples (molecules)
# tau1 = decay time 1
# tau[] = list of decay times 2 
# trials = number of trials computed for each set of paramters (tau2,N,order)
# order[] = list of orders of method
#
# results will be in 3D arrays
# e.g., meanA[i,j,k] where 
# i = 0,1,...  index for tau2
# j = 0,1,...  index for N (number of samplep
# k = 0,1,...  index for order

# decay constants
tau = [ [10.0, 10.0, 10.0], [10.0, 10.0, 50.0], [10.0, 30.0, 100.0] ]

# try for many values of N = number of samples
N = [5, 10, 20, 50, 100, 200, 500, 1000]

# try for different orders
order = [ 3,4 ]

# for each set of parameters (tau2,N,order) compute trA1fls:
trials = 1000
initialtau = [ [ 1.0, 1.0, 1.0 ], [ 1.0, 1.0, 10.0 ], [ 1.0, 1.0, 100.0 ], [ 1.0, 1.0, 1000.0 ], [ 1.0, 10.0, 10.0 ], [ 1.0, 10.0, 100.0 ], [ 1.0, 10.0, 1000.0 ], [ 1.0, 100.0, 100.0 ], [ 1.0, 100.0, 1000.0 ], [ 1.0, 1000.0, 1000.0 ], [ 10.0, 10.0, 10.0 ], [ 10.0, 10.0, 100.0 ], [ 10.0, 10.0, 1000.0 ], [ 10.0, 100.0, 100.0 ], [ 10.0, 100.0, 1000.0 ], [ 10.0, 1000.0, 1000.0 ], [ 100.0, 100.0, 100.0 ], [ 100.0, 100.0, 1000.0 ], [ 100.0, 1000.0, 1000.0 ], [ 1000.0, 1000.0, 1000.0 ]]

# initialize arrays to hold results
# these hold results for tauA
tempA = zeros(trials)
meanA = zeros([len(tau),len(N),len(order)])
meandevA = zeros([len(tau),len(N),len(order)])
stdA = zeros([len(tau),len(N),len(order)])
seA = zeros([len(tau),len(N),len(order)])
# these hold results for tauB
tempB = zeros(trials)
meanB = zeros([len(tau),len(N),len(order)])
meandevB = zeros([len(tau),len(N),len(order)])
stdB = zeros([len(tau),len(N),len(order)])
seB = zeros([len(tau),len(N),len(order)])
# these hold results for tauC
tempC = zeros(trials)
meanC = zeros([len(tau),len(N),len(order)])
meandevC = zeros([len(tau),len(N),len(order)])
stdC = zeros([len(tau),len(N),len(order)])
seC = zeros([len(tau),len(N),len(order)])

for i in range(len(tau)): # loop through tau values
    print "now calculating tau2 =",tau[i]
    for j in range(len(N)): # loop through numbers of samples
        print "   ...N =", N[j]
        for k in range(len(order)):
            print "      ...order =", order[k]
            for n in range(trials): 
                times = sim.multi_poissonN(tau[i],N[j])

                junk=gmm.gmmG(times,initialtau,n=order[k],diag=True,weight='jack',bc=True )
                tempA[n]=junk[0]  # this is smallest decay time
                tempB[n]=junk[1]  # this is middle decay time
                tempC[n]=junk[2]  # this is largest decay time

            meanA[i,j,k],meandevA[i,j,k],stdA[i,j,k],seA[i,j,k] = cu.stats(tempA)
            meanB[i,j,k],meandevB[i,j,k],stdB[i,j,k],seB[i,j,k] = cu.stats(tempB)
            meanC[i,j,k],meandevC[i,j,k],stdC[i,j,k],seC[i,j,k] = cu.stats(tempC)
 









