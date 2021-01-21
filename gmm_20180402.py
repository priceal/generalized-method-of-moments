############################################################
############################################################
#############   gmm.py   ##########################
##############################################################
############################################################
#
# Generalized Method of Moments module v20170926
#
# This module is used to perform calculations reported in
# "Analyzing Dwell Times with the Generalized Method of Moments",
# by Sadie Piatt and Allen Price, submitted to PLOS ONE.
#
############################################################
############################################################

import numpy as np
from scipy.optimize import leastsq, minimize
import cumulant as cu

##############################################################
#############   residual()   ##########################
##############################################################
#
# returns residuals of measured cumulants versus theoretical
# from supplied decay times for n step process:
#
# A -> B -> C -> ...
# tA   tB   tC   ...
#
# order is determined by number of cumulants supplied in k
# number of steps is determined by number of paramters in tau
#
# input:
# tau	[tauA, tauB, ...] = lifetimes
# k	[cm1, cm2, cm3, ..., cmN], the sample cumulants up to order N (max 6)
# n	[a,b,c,d,e,f] where a,b,.. are each boolean values indicating if the
#        corresponding cumulant is to be used.  Note that a+b+... must equal
#        the length of k or error will occur.
#
# output:
# residual		vector of residuals
#
def residual(tau, k, n):
    """returns residuals up order O(6)"""
        
    tau = np.array(tau)
    k = np.array(k)
    n = np.array(n)

    # calculate the theoretical cumulants from tau 
    cm1 = tau.sum()    				# <T> = tA + tB + ...
    cm2 = (tau*tau).sum()			# <DT^2> = tA^2 + tB^2 + ...
    cm3 = 2.0*(tau*tau*tau).sum()		# <DT^3> = 2tA^3 + 2tB^3 + ...
    cm4 = 6.0*(tau*tau*tau*tau).sum()		# <DT^4> = 6tA^4 + 6tB^4 + ...
    cm5 = 24.0*(tau*tau*tau*tau*tau).sum()	# <DT^5> = 24tA^5 + 24tB^5 + ...
    cm6 = 120.0*(tau*tau*tau*tau*tau*tau).sum()	# <DT^6> = 120tA^6 + 120tB^6 + ...

    theory = np.array([cm1,cm2,cm3,cm4,cm5,cm6])

    return theory[n]-k

##############################################################
#############   cost()   ##########################
##############################################################
#
# returns GMM cost function of measured cumulants versus theoretical
# for n step process:
#
# A -> B -> C ...
# tA   tB   tC  ...
#
# order is determined by number of cumulants supplied in k
# number of steps is determined by number of paramters in tau
# note: order must be > or = to number of steps!
#
# input:
# tau 	lifetimes, in form [tauA,tauB,...]
# k	array of measured cumulants up to order 6
# w	symmetric, positive definite weight matrix, must have dimensions NxN,
#        where N = order of k
# n	[a,b,c,d,...] where a,b,.. are each boolean values indicating if the
#        corresponding cumulant is to be used.  Note that a+b+c+... must equal
#        the length of k or error will occur.
#
# output:
# cost		scalar, non-negative cost
#
def cost(tau,k,w,n):
    """GMM cost function for n step process"""

    tau = np.array(tau)	# decay times
    k = np.array(k)    	# measured cumulants from sample
    w = np.array(w)    	# positive definite weight matrix
    n = np.array(n)  	# boolean values indicating which cumulants to use

    g = residual(tau,k,n)  # residuals

    return np.dot(g,np.dot(w,g))  # cost = gWg

##############################################################
#############   dcost()   ##########################
##############################################################
#
# returns derivate of GMM cost function above.
# input is same
#
# output:
# dcost		vector derivative O(n)
#
def dcost(tau,k,w,n):
    """derivative of GMM cost function for n step process"""

    tau = np.array(tau)	# decay times
    k = np.array(k)    	# measured cumulants from sample
    w = np.array(w)    	# positive definite weight matrix
    n = np.array(n)

    g = residual(tau,k,n)  # residuals

    # this is the derivative of the residuals = n!t^(n-1)
    dg = np.ones([6,len(tau)])
    dg[1] = 2.0*tau           		# second order
    dg[2] = 6.0*tau*tau       		# third
    dg[3] = 24.0*tau*tau*tau  		# fourth
    dg[4] = 120.0*tau*tau*tau*tau  	# fifth
    dg[5] = 720.0*tau*tau*tau*tau*tau  	# sixth
    dg = dg[n]
    
    return 2.0*np.dot(g,np.dot(w,dg)) # dcost = 2gWg'

##############################################################
#############   gmm()   #####################################
##############################################################
#
def gmm(t,tau0,n=1,diag=False,weight='jack',verbose=False,bc=True):
    """GMM for N step process, with a weight matrix, number of
    steps is determined by length of tau0=[tau10,tau20, ...]

    input: 
    t        array of sample times
    tau0     initial guess time constants [tau10, tau20, ...] 
    n        order of method up to 6 (default = number of steps)
             or an array 1/0 values indicating which cumulants to use.
             for example, [0,1,1,0,0,0] indicates use cumulants 2 and 3
    diag     if True, use diagonalized covariance matrix
    weight   'jack': estimate covariance with jackknife method
             'mc': use monte-carlo method to calculate covariance
             'int': use interpolation method to calculate covariance
             'iden': sets weight = identity matrix
    verbose  set to True if you want cost function value returned also
             this is necessary for global searches
    output:
    tau      estimates of decay times [tau1, tau2, ...]"""
  
#   process inputs
    times = np.array(t)
    tau0 = np.array(tau0) 
    numparams = len(tau0) 
    if type(n) == int:
        if n<numparams:
            n = np.concatenate( [np.ones(numparams),np.zeros(6-numparams)] )
        else:
            nmax = min(n,6)
            n = np.concatenate( [np.ones(nmax),np.zeros(6-nmax)] )
    else:
        n = np.array(n)
    order = n.sum()
    if order<numparams:
        n = np.concatenate( [np.ones(numparams),np.zeros(6-numparams)] )

#   calculate cumulants and weights
    if weight=='jack':  # use jackknife estimate
        k, kcov = cu.cumulants(times,n=n,jack=True,bc=bc)
    elif weight=='mc':  # use montecarlo method
        k = cu.cumulants(times,n=n,jack=False,bc=bc)
        kcov = cu.kcov(tau0,len(t),trials=500,n=n)
    elif weight=='int': # use interpolation method
        k = cu.cumulants(times,n=n,jack=False,bc=bc)
        kcov = cu.kcovint(tau0,len(t),n=n)
    elif weight=='iden': # set weight=identity matrix
        k = cu.cumulants(times,n=n,jack=False,bc=bc)
        kcov = np.identity(order)

    if diag==True:
        w = np.diag(1.0/np.diag(kcov))
    elif diag==False:
        w = np.linalg.inv(kcov)

#   perform minimization
    n = [bool(n[0]),bool(n[1]),bool(n[2]),bool(n[3]),bool(n[4]),bool(n[5])]  # set boolean values before minimization
    result = minimize(cost,tau0,args=(k,w,n),method='BFGS',jac=dcost,options={'gtol': 1e-8, 'disp': False})

#   return result and cost function value minimum if needed
    if verbose == True:
        return np.sort(result['x']),result['fun']
    else:
        return np.sort(result['x'])

##############################################################
#############   gmmG()   #####################################
##############################################################
#
def gmmG(t,taulist,n=1,diag=False,weight='jack',bc=True):
    """a global search wrapper for gmm.  see gmm.gmm() for complete
    list of input.  Here, input taulist is a list of initial values
    to try.  Number of steps is determined from list.  

    input: 
    taulist	list of initial values to try.  For example, for 
		a two step could be something like:
		[ [10,10],[10,20],[10,30],[20,20],[20,30],[30,0] ]
		this example would perform 6 minimizations and
		return the results that yield the overall minimum
    output:
    tau     	estimates of decay times [tau1, tau2] which
		minimizes cost function in region specified"""

    taulist = np.array(taulist)
    numpoints = len(taulist)
    nsteps = len(taulist[0])
    decay = np.zeros([numpoints,nsteps])
    value = np.zeros([numpoints])
    i = 0
    for tau in taulist:
        result=gmm(t,tau,n=n,diag=diag,weight=weight,verbose=True,bc=bc) 
	result[0].sort()
        decay[i] = result[0]
        value[i] = result[1]
        i = i + 1
            
    minindex = np.argmin(value)
    return decay[minindex]



















