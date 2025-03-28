import re
import os
import shutil
import subprocess
import tinygp
import numpy as np
import matplotlib.pyplot as plt
import jax.numpy as jnp
import jax.random as jr
# Modified statistical distribution
def stat_dist(alpha,l):
    return (2 *l + 1)*np.exp(alpha*l)
def statlinear_dist(l):
    return (2 *l + 1)
# Vogel's quadratic distribution
def quad_dist(a,b,l):
    return 1 + a*l + b*l**2

def log_dist(A,B,x):
    return A*np.log(B*(x + np.exp(1/A)/B))
# Routine to normalise the input distribution
# Currently not in use
def norm_func(dist, alpha, npoints, X):

    sum = 0
    if(dist == "stat"):
        for i in range(0,npoints):
            sum += stat_dist(alpha, X[i])
    elif(dist == "quad"):
        for i in range(0,npoints):
            sum += quad_dist(a, b, X[i])
    return sum

# This routine takes optimised parameter(s) and produces input
# files using GP samples. At the moment, the kernel is fixed
def run_gaussian_process(*params, element, nsamples, l_distribution, signal, scale, nmax):
    key = jr.PRNGKey(123)

    # print(signal, scale, l_distribution)

    # Choose an appropriate kernel
    # k = 0.01 * tinygp.kernels.ExpSquared(scale=0.01, distance=tinygp.kernels.L2Distance())
    # k = 1e-15 * tinygp.kernels.RationalQuadratic(scale=1e-15, distance=tinygp.kernels.L2Distance())
    k =float(signal) * tinygp.kernels.ExpSquared(scale=float(scale))
    # k = float(signal)*tinygp.kernels.Matern32(scale=scale)

    # Make the l axis
    npoints  = nmax
    X_p = np.linspace(0,nmax-1, npoints)

    # Build our prior GP
    prior_gp = tinygp.GaussianProcess(k, X_p)

    # Statistical
    if(int(l_distribution)==0):
        alpha = params[0]
        distribution = "stat"
	# Quadratic
    elif(int(l_distribution) == 2):
        a = params[0]
        b = params[1]
        distribution = "quad"
    elif(int(l_distribution) == 4):
        distribution = "statlinear"
    elif(int(l_distribution) == -1):
        distribution = "custom"
        custom_distribution = params[0]
    elif(int(l_distribution) == 16):
        a = params[0]
        b = params[1]
        distribution = "log"

    # Obtain all of our samples
    samples = prior_gp.sample(key=key, shape=(nsamples,)).T
    mean_vals = np.zeros(npoints)

    # normalisation = norm_func(distribution, alpha, npoints, X_p)

    # For now, just keep it unnormalised
    normalisation = 1

    # This will contain our initial population for each sample
    val = np.zeros((nsamples, npoints))

    # Loop over number of samples
    for j in range(nsamples):

        # Loop over each l value
        for i in range(0,nmax):

            # Choose distribution
            if(distribution == "quad"):

                # If our population goes below zero, force it to zero
                if((samples[i,j] + quad_dist(a,b,X_p[i])/normalisation) < 0):
                   val[j,i] = 0
                
                # Otherwise, add the sample to function evaluation
                else:
                   val[j,i] = (samples[i,j] + quad_dist(a,b,X_p[i])/normalisation)
            elif(distribution == "log"):

                # If our population goes below zero, force it to zero
                if((samples[i,j] + log_dist(a,b,X_p[i])/normalisation) < 0):
                   val[j,i] = 0
                
                # Otherwise, add the sample to function evaluation
                else:
                   val[j,i] = (samples[i,j] + log_dist(a,b,X_p[i])/normalisation)
									
            elif(distribution == "stat"):

                # If our population goes below zero, force it to zero
                if((samples[i,j] + stat_dist(alpha,X_p[i])/normalisation) < 0):
                   val[j,i] = 0
                
                # Otherwise, add the sample to function evaluation
                else:
                   val[j,i] = (samples[i,j] + stat_dist(alpha,X_p[i])/normalisation)
            elif(distribution == "statlinear"):

                # If our population goes below zero, force it to zero
                if((samples[i,j] + statlinear_dist(X_p[i])/normalisation) < 0):
                   val[j,i] = 0
                
                # Otherwise, add the sample to function evaluation
                else:
                   val[j,i] = (samples[i,j] + statlinear_dist(X_p[i])/normalisation)
            elif(distribution == "custom"):

                # If our population goes below zero, force it to zero
                if((samples[i,j] + custom_distribution[i]/normalisation) < 0):
                   val[j,i] = 0
                
                # Otherwise, add the sample to function evaluation
                else:
                   val[j,i] = (samples[i,j] + custom_distribution[i]/normalisation)

    # f = open("statsamples.dat", "w")
    # for i in range(20):
    #     line = ' '.join(map(str, val[:,i]))
    #     f.write(str(X_p[i])+" "+line+"\n")
    # plt.show()
    # Loop over samples again
    for i in range(nsamples):

        # Make a new input file for the GP sample 
        shutil.copyfile('basic_input', "input_"+element+"_"+str(i+1)+"_"+str(signal)+"_"+str(scale))
        f = open("input_"+element+"_"+str(i+1)+"_"+str(signal)+"_"+str(scale), "a")
        
        # Write the population to input file
        for j in range(nmax):
            f.write("PL    FR   "+str(j)+" "+ str(val[i,j])+"\n")

        # Finish up
        f.close()
        
