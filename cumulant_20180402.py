############################################################
############################################################
#############   cumulants.py ##########################
##############################################################
############################################################
#
# Cumulants module v20170926
#
# This module is used to perform calculations reported in
# "Analyzing Dwell Times with the Generalized Method of Moments",
# by Sadie Piatt and Allen Price, submitted to PLOS ONE.
#
############################################################
############################################################

import numpy as np
import sim
import cPickle as pickle

############################################################
############# define global variable
##################################################

GrandMatrix = np.zeros([1,1,1])
NumSamples = [1]
TauTwo = [1.0]

##############################################################3
############## read in the grand matrix
###############################################################
#
def initialize(filename="grandmatrix.pkl"):
    """initialize the grand matrix"""

    global GrandMatrix
    global NumSamples
    global TauTwo

    junk = pickle.load( open( filename, "rb" ) )

    NumSamples = junk[0]
    TauTwo = junk[1]
    GrandMatrix = junk[2]

    return True


def printGrandMatrix(index=0):

    print "NumSamples =", NumSamples
    print "TauTwo = ", TauTwo
    print "size of grand matrix = ", len(GrandMatrix), "X", len(GrandMatrix[0]), "X", len(GrandMatrix[0,0]), "X", len(GrandMatrix[0,0,0])

    if index != 0:
        print "covariance ", index[0], index[1]
        print "N, tauB =", NumSamples[index[0]], TauTwo[index[1]]
        print GrandMatrix[index[0],index[1]]

    return True


##############################################################
#############   stats()   ##########################
##############################################################
#
# calculate simple statistics from a list of wait times
#
# input:
# times		array or list of wait times
#
# output:
# average, mean deviation, standard deviation, standard error in mean
#
def stats(junk):

    times=np.array(junk)
    n = len(times)
    average = np.mean(times)
    standard = np.std(times)
    deviation = np.abs(times-average)
    meandev = np.mean(deviation)
    stderr = standard/np.sqrt(n)

    return average, meandev, standard, stderr

##############################################################
#############   cumulants()   ##########################
##############################################################
#
# calculate the (biased) cumulants from a sample
#
# input:
# sample	array or list of wait times
# n		number of orders, up to 6 (default=4)
# jack		determine jackknife estimate of uncertaity? (default:False)
#
# output:
# k		if jack=False
# k, kcov	if jack=True
#
# k = [mean, 2nd, 3rd, ... up to 6th cumulant]
# kcov = covariance in [mean, 2nd, 3rd, ... up to 6th cumulant]
#
def cumulants(sample,n=4,jack=False,bc=True):

#   process arguments
    if type(n) == int:
        n = np.concatenate( [np.ones(n),np.zeros(6-n)] )
    n = np.array([bool(n[0]),bool(n[1]),bool(n[2]),bool(n[3]),bool(n[4]),bool(n[5])])
    N = len(sample)

#   calculate cumulants
    rm1 = np.mean(sample)  # raw first moment
    deviation = sample - rm1  
    cm2 = np.mean(deviation*deviation)
    cm3 = np.mean(deviation*deviation*deviation)
    cm4 = np.mean(deviation*deviation*deviation*deviation)
    cm5 = np.mean(deviation*deviation*deviation*deviation*deviation)
    cm6 = np.mean(deviation*deviation*deviation*deviation*deviation*deviation)
    
    k = np.array([rm1,0,0,0,0,0])
    if bc==True:
        k[1] = cm2*N/(N-1.0)
        k[2] = cm3*N*N/(N-1.0)/(N-2.0)
        k[3] = ( cm4*(N+1.0) - 3.0*cm2*cm2*(N-1.0) )*N*N/(N-1.0)/(N-2.0)/(N-3.0)
        k[4] = cm5 - 10.0*cm3*cm2       
        k[5] = cm6 - 15.0*cm4*cm2 - 10.0*cm3*cm3 + 30.0*cm2*cm2*cm2    
        if jack==True:
            # initialize the jackknife cumulants
            jackknife_k1 = np.zeros(N)
            jackknife_k2 = np.zeros(N)
            jackknife_k3 = np.zeros(N)
            jackknife_k4 = np.zeros(N)
            jackknife_k5 = np.zeros(N)
            jackknife_k6 = np.zeros(N)

            # calculate the jackknife sample cumulants out to order 6
            for i in range(N):
                jackknife_sample = np.delete(sample,i)
                rm1 = np.mean(jackknife_sample)  
                deviation = jackknife_sample - rm1
                cm2 = np.mean(deviation*deviation)
                cm3 = np.mean(deviation*deviation*deviation)
                cm4 = np.mean(deviation*deviation*deviation*deviation)
                cm5 = np.mean(deviation*deviation*deviation*deviation*deviation)
                cm6 = np.mean(deviation*deviation*deviation*deviation*deviation*deviation)
                jackknife_k1[i] = rm1
                jackknife_k2[i] = cm2*N/(N-1.0)
                jackknife_k3[i] = cm3*N*N/(N-1.0)/(N-2.0)
                jackknife_k4[i] = ( cm4*(N+1.0) - 3.0*cm2*cm2*(N-1.0) )*N*N/(N-1.0)/(N-2.0)/(N-3.0)
                jackknife_k5[i] = cm5 - 10.0*cm3*cm2      
                jackknife_k6[i] = cm6 - 15.0*cm4*cm2 - 10.0*cm3*cm3 + 30.0*cm2*cm2*cm2    

            # calculate the jackknife estimate of the covariance matrix
            kcov = np.cov([jackknife_k1,jackknife_k2,jackknife_k3,jackknife_k4,jackknife_k5,jackknife_k6],bias=True)[n][:,n]
            k    = k[n], kcov  # using bias=True gives 1/N normalization for jackknife estimate
        else:
            k = k[n]

    else:
        k[1] = cm2
        k[2] = cm3
        k[3] = cm4 - 3.0*cm2*cm2
        k[4] = cm5 - 10.0*cm3*cm2
        k[5] = cm6 - 15.0*cm4*cm2 - 10.0*cm3*cm3 + 30.0*cm2*cm2*cm2

        if jack==True:
            # initialize the jackknife cumulants
            jackknife_k1 = np.zeros(N)
            jackknife_k2 = np.zeros(N)
            jackknife_k3 = np.zeros(N)
            jackknife_k4 = np.zeros(N)
            jackknife_k5 = np.zeros(N)
            jackknife_k6 = np.zeros(N)

            # calculate the jackknife sample cumulants out to order 6
            for i in range(N):
                jackknife_sample = np.delete(sample,i)
                rm1 = np.mean(jackknife_sample)  
                deviation = jackknife_sample - rm1
                cm2 = np.mean(deviation*deviation)
                cm3 = np.mean(deviation*deviation*deviation)
                cm4 = np.mean(deviation*deviation*deviation*deviation)
                cm5 = np.mean(deviation*deviation*deviation*deviation*deviation)
                cm6 = np.mean(deviation*deviation*deviation*deviation*deviation*deviation)
                jackknife_k1[i] = rm1
                jackknife_k2[i] = cm2
                jackknife_k3[i] = cm3
                jackknife_k4[i] = cm4 - 3.0*cm2*cm2
                jackknife_k5[i] = cm5 - 10.0*cm3*cm2      
                jackknife_k6[i] = cm6 - 15.0*cm4*cm2 - 10.0*cm3*cm3 + 30.0*cm2*cm2*cm2    

            # calculate the jackknife estimate of the covariance matrix
            kcov = np.cov([jackknife_k1,jackknife_k2,jackknife_k3,jackknife_k4,jackknife_k5,jackknife_k6],bias=True)[n][:,n]
            k    = k[n], kcov  # using bias=True gives 1/N normalization for jackknife estimate
        else:
            k = k[n]
 
    return k

##############################################################
#############   kcov()   ##########################
##############################################################
#
# calculate the theoretical covariance matrix of the
# sample cumulants using monte-carlo method.  Used in GMM
#
# input:
# tau		[tau1, tau2, ...]
# N		sample size
# trials	number of trials in calculation (default 500)
# n		order of cumulants requested (default 4)
#
# output:
# kcov		the cumulant-covariance matrix
#
def kcov(tau,N,trials=500,n=4):

#   process arguments
    if type(n) == int:
        n = np.concatenate( [np.ones(n),np.zeros(4-n)] )
    n = np.array([bool(n[0]),bool(n[1]),bool(n[2]),bool(n[3])])
    tau = np.array(tau)

    k = np.zeros([trials,4])
    for m in range(trials):
        times = sim.multi_poissonN(tau,N)
        k[m] = cumulants(times,n=4,jack=False,bc=True)

    return np.cov(k.transpose(),bias=False)[n][:,n]

##############################################################
#############   kcovint()   ##########################
##############################################################
#
# calculate the theoretical covariance matrix of the
# sample cumulants using interpolation method.  
# only works for two step reaction scheme, and for specific 
# values of N
#
# input:
# tau		[tau1, tau2]
# N		sample size
# n		order of cumulants requested (default 4)
#
# output:
# kcov		the cumulant-covariance matrix
#
def kcovint(tauIn,N,n=4):

#   process arguments
    if type(n) == int:
        n = np.concatenate( [np.ones(n),np.zeros(4-n)] )
    n = np.array([bool(n[0]),bool(n[1]),bool(n[2]),bool(n[3])])
    tau = np.zeros(2)
    tau[0] = np.min(tauIn)
    tau[1] = np.max(tauIn)
    s = tau[1]/tau[0]	#rescaled decay times
    
    # first convert N value to an index
    Nindex = NumSamples.index(N)
    
    # now convert s to a pair of indices spanning s
    # assume grandmatrix has s = 1,2,3,4,etc...
    tIndexHi = int(np.floor(s))
    tIndexLow = tIndexHi - 1
    weightHi = s-np.floor(s)
    weightLow = 1.0 - weightHi
    if tIndexLow < 0:
        tIndexHi = 0
        tIndexLow = 0
    if tIndexHi > (len(TauTwo)-1) :
        tIndexHi = len(TauTwo)-1
        tIndexLow = len(TauTwo)-1
       
    # here is the linear interpolation (maybe geometric interpolation better?)
    kcovRAW = weightLow*GrandMatrix[Nindex,tIndexLow] + weightHi*GrandMatrix[Nindex,tIndexHi]
    
    # now calculate scale factors
    factor = np.zeros([4,4])
    for i in range(4):
        for j in range(4):
            factor[i,j]=tau[0]**(i+j+2)

    return (kcovRAW * factor)[n][:,n]

