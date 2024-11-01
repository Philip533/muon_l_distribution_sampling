#!/usr/env/python3
import sys
import seaborn as sns
import re
import math
import matplotlib.pyplot as plt
import shutil
import os
import numpy as np
import subprocess
import scipy.stats as st

# Included in the repository
from element_information import return_elements
from gaussianprocess import run_gaussian_process

# 2l + 1
def statlinear_dist(l):
    return (2 *l + 1)

# (2l+1)exp(alpha * l)
def stat_dist(alpha,l):
    return (2 *l + 1)*np.exp(alpha*l)

# 1 + al + bl^2
def quad_dist(a,b,l):
    return 1 + a*l + b*l**2

# Routine to read the experimental data from custom 
# defined format as show in exp_intens.dat
# We also return the number of s transitions
# i.e n->n-1
def read_exp_data(element_name):

    # Open file and read all the lines
    f = open("exp_intens.dat", "r")
    lines = f.readlines()

    # Start counting each s transition
    s_trans = 0

    # This dictionary contains all transitions from the key state
    intensity_dict = {}
    counter = 0

    # Loop over each line
    for line in lines:

		# Split into each transition
        splitline = line.split()

        # Find the appropriate element and count how many starting
        # energy levels we have
        if(splitline[0] == element_name):
          counter += 1

		  # Start point is the energy level in which we start
          start_point = splitline[1]

		  # Put the entire line of transitions into the dictionary
          intensity_dict[start_point] = splitline[2:]

          # Count if it's an s transition
          if(str(splitline[2:][int(start_point)-2])) != "-2+":
              s_trans += 1
                         
    return intensity_dict, s_trans

# Return an array containing a normal distribution for a given
# experimental mean and std
def generate_experimental_gaussian(mean, std, linspace):
    
    # Values of Gaussian
    gaussian = []

    # Loop over the intensity (x) axis and calculate the value of the Gaussian
    for i in linspace:
        gaussian.append(1/(math.sqrt(2*np.pi*std**2)) * np.exp((-(i - mean)**2)/(2*std**2)))

    return np.array(gaussian)

# Get all of the transitions and put them in a dictionary
def parse_akylas_transitions(out_name):

    # This dictionary holds all transitions starting at i, in the key i
    akylas_intensities = {}

    # We want to loop over our entire set of transitions given by Akylas
    # which is from 2 to 20
    for i in range(2,21):

        # Make an empty list in the dictionary for the starting point
        akylas_intensities[i] = []

        # Now loop over the end points
        for j in range(1,i):

            # Here we need a space between the = and number for both states
            if ((i <= 9) and (j <= 9)):

                # Searches for the appropriate transition in the output file
                line = subprocess.check_output(["grep", "N1= "+str(i)+", N2= "+str(j)+",", out_name])

                # Obtain the intensity
                float_line = float(line.decode("utf-8").split()[7].split("(")[0])

            # Here there is a space only for the end state
            elif ((i > 9) and (j <= 9)):

                # Searches for the appropriate transition in the output file
                line = subprocess.check_output(["grep", "N1="+str(i)+", N2= "+str(j)+",", out_name])

                # Obtain the intensity
                float_line = float(line.decode("utf-8").split()[6].split("(")[0])

            # Neither state has a gap between = and the number
            elif ((i > 9) and (j > 9)):

                # Searches for the appropriate transition in the output file
                line = subprocess.check_output(["grep", "N1="+str(i)+", N2="+str(j)+",", out_name])

                # Obtain the intensity
                float_line = float(line.decode("utf-8").split()[5].split("(")[0])

            # Add each transition intensity to the appropriate starting state entry
            akylas_intensities[i].append(float_line)

    return akylas_intensities
    
# If no element is provided then kill the program
if(len(sys.argv) == 1):
    print("Fatal error: Element must be provided at command line")
    sys.exit()

# For now we are keeping nmax at 20
nmax = 20

# Get the element and l-distribution from command line
element = sys.argv[1]
l_distribution = sys.argv[2]

# Get all the experimental data, but it's not split 
# into intensities and errors yet
exp_intens = {}
exp_intens, num_s_transitions = read_exp_data(element)

# First let's generate an array containing the mean population
distribution = np.zeros(nmax)

# Custom distribution used
# Unused behaviour for now
if(int(l_distribution) == -1):
    print ("Custom l-distribution chosen.")
    nop_line = "NOP   OR   "+str(l_distribution)
    experiment_name = element+"custom"+alpha

# Modified statistical distribution
elif(int(l_distribution) == 0):
    print ("Modified statistical distribution chosen")
    if(len(sys.argv) <= 3):
        print ("Fatal error: Exponent must be provided")
        sys.exit()
    else:
        alpha = sys.argv[3]
        print ("Exponent alpha = ", alpha)
        nop_line = "NOP   OR   "+str(l_distribution)
        experiment_name = element+"_stat_"+alpha
        alpha = float(alpha)

        # Generates a mean distribution using the given alpha
        for i in range(nmax):
            distribution[i] = stat_dist(alpha, i)

# Quadratic distribution
elif(int(l_distribution) == 2):
    print ("Quadratic distribution chosen")
    if (len(sys.argv) <= 4):
        print ("Fatal error: Two coefficients must be provided")
        sys.exit()
    else:
        a, b = sys.argv[3:]
        print ("Coefficients a and b = ", a, b)
        nop_line = "NOP   OR   "+str(l_distribution)
        experiment_name = element+"_quad_"+a+"_"+b
        a = float(a)
        b = float(b)

        # Generate the mean distribution
        for i in range(nmax):
            distribution[i] = quad_dist(a,b,i)

# Regular statistical distribution
elif(int(l_distribution) == 4):
    print ("Statistical distribution chosen")
    print ("Parameter free")
    nop_line = "NOP   OR   "+str(l_distribution)
    nmx_line = "NMX   OR   "+str(nmax)+" 0"
    experiment_name = element+"_statlinear"
    for i in range(nmax):
        distribution[i] = statlinear_dist(i)

else:
    print("No valid distribution chosen")
    sys.exit()

# Name of the folder which will contain the results
os.makedirs(experiment_name, exist_ok=True)

# This gets all our element specific parameters for Akylas
# from data tabulated in the electron_binding_energies folder
z_line, zs_line, a_line, be_line = (return_elements(element))

# Compile the code here
os.system("gfortran -g --std=legacy -fd-lines-as-comments -ffpe-trap=zero,invalid,denorm,underflow -fbacktrace -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8 -funsafe-math-optimizations -fno-align-commons -fmax-errors=1 muon.f")

# More necessary lines for the input file
# which are required for every run
stop_line = "STO(P)NE"
xeq_line = "XEQ   NE"

# Move basic file to directory
shutil.copyfile("basic_input", experiment_name+"/basic_input")

# Move into the proper directory
os.chdir(experiment_name)

# Number of samples to use in each degree of freedom
# for the sensitivity analysis
iters = 100

# This dictionary will contain all possible n level transitions
# from Akylas
intensity_table = []

# This contains a dictionary for each sample of the parameters
list_of_dicts = []

# Counter for keeping track of the samples
counter = 0

# Loop over a set of bounds from the user specified exponent
if(int(l_distribution) == 0):

	# Set the range we are looking at
    sensitivity_range = 3.0
    upper_bound = alpha * ( 1 + sensitivity_range)
    lower_bound = alpha * ( 1 - sensitivity_range)

	# Zero the alpha array
    alpha_vals = np.zeros(iters)

	# Loop over our range of parameter space
    for i in np.arange(lower_bound, upper_bound, (upper_bound - lower_bound)/iters):

		# Round for naming the files
        j = np.round(i,3)
        name = "input_"+experiment_name+"_"+str(j)
        out_name = "output_"+experiment_name+"_"+str(j)
        nmx_line = "NMX   OR   "+str(nmax)+" "+str(i)

		# Copy the basic input and append all of the relevant
		# input lines to it
        shutil.copyfile("basic_input", name)
        f = open(name, "a")
        f.write(z_line+"\n")
        f.write(zs_line+"\n")
        f.write(a_line+"\n")
        f.write(be_line+"\n")
        f.write(nop_line+"\n")
        f.write(nmx_line+"\n")
        f.write(xeq_line+"\n")
        f.write(stop_line+"\n")
        f.close()
		
		# Run the code with our chosen distribution
        os.system("../a.out < "+name+" >"+out_name)
		
		# Collect all of the intensities
        intensity_table = parse_akylas_transitions(out_name)
		  
        # Keep track of the intensities and parameter value
        alpha_vals[counter] = j
        list_of_dicts.append(intensity_table)
        counter += 1
 
elif(int(l_distribution) == 2 or int(l_distribution) == 16):

  sensitivity_range = 0.8
  a = float(a)
  b = float(b)
  upper_bound1 = a * ( 1 + sensitivity_range)
  lower_bound1 = a * ( 1 - sensitivity_range)
  upper_bound2 = b * ( 1 + sensitivity_range)
  lower_bound2 = b * ( 1 - sensitivity_range)

  a_vals = np.zeros(iters**2)
  b_vals = np.zeros(iters**2)
  i_count = 0
  j_count = 0

  # Loop over a values
  for i in np.arange(lower_bound1, upper_bound1, (upper_bound1 - lower_bound1)/iters):

	# Loop over b values
    for j in np.arange(lower_bound2, upper_bound2, (upper_bound2 - lower_bound2)/iters):

	    # Round for nicer file names
        k = np.round(i,3)
        l = np.round(j,3)

        # Name everythin
        name = "input_"+experiment_name+"_"+str(k)+"_"+str(l)
        out_name = "output_"+experiment_name+"_"+str(k)+"_"+str(l)


        # Set alpha to 0
        nmx_line = "NMX   OR   "+str(nmax)+" 0"

        # Write the quadratic parameters
        cl_line = "CL   DE   "+str(i)+" "+str(j)
					
        shutil.copyfile("basic_input", name)
        f = open(name, "a")
        f.write(z_line+"\n")
        f.write(zs_line+"\n")
        f.write(a_line+"\n")
        f.write(be_line+"\n")
        f.write(nop_line+"\n")
        f.write(nmx_line+"\n")
        f.write(cl_line+"\n")
        f.write(xeq_line+"\n")
        f.write(stop_line+"\n")
        f.close()

        # Run the code
        os.system("../a.out < "+name+" >"+out_name)

        # Obtain all of our intensities
        intensity_table = parse_akylas_transitions(out_name)
        a_vals[counter] = i
        b_vals[counter] = j
        list_of_dicts.append(intensity_table)
        counter += 1
# Linear statistical distribution
elif(int(l_distribution) == 4):
    print("No sensitivity analysis required")
# Square number of iterations for the 2 parameter fit
if(int(l_distribution) == 2 or int(l_distribution) == 16):
    iters = iters**2

if(int(l_distribution) != 4):
    # Loss as a function of parameter space
    loss_array = np.zeros(iters)

    # Loop over each of our samples
    for l in range(iters):
        loss_func = 0

        # Loop over each experimental start point
        for i,j in zip(exp_intens.values(), exp_intens.keys()):

            counter = 1

            # Loop over each transition 
            for k in i:

                # Split up into intensity and error
                intens = k.split("+")[0]
                err = k.split("+")[1]

                # If we reached a -2, then we can move on to the next line as
                # no experimental transition exists
                if(float(intens) < -1):
                    continue

                # Add the contribution to the loss function, scaled with the 
                # experimental uncertainty
                loss_func += ((float(list_of_dicts[l][int(j)][counter-1]) - float(intens))/float(err))**2
                counter += 1
            loss_array[l] = loss_func

# Find the fitted parameter which corresponds to the minimum
# value of the loss function
if(int(l_distribution) == 0):
	fitted_alpha = alpha_vals[np.argmin(loss_array)]
	print("Optimal parameter = ", fitted_alpha)
elif(int(l_distribution) == 2 or int(l_distribution) == 16):
	# Find both of the fitted parameters
    fitted_a= a_vals[np.argmin(loss_array)]
    fitted_b= b_vals[np.argmin(loss_array)]
    print("Optimal parameters = ", fitted_a, fitted_b)


# Make a new directory to contain data we will
# want to plot with GNUplot
figure_data = "figure_data"
os.makedirs(figure_data, exist_ok=True)

# Name of the lossplot file
loss_data = "lossplot_"+experiment_name+".dat"

if(int(l_distribution) != 4):
    # Write the loss function to file
    f = open(loss_data, "w")
    if(int(l_distribution) == 0):
        for i in range(len(alpha_vals)):
            f.write(str(alpha_vals[i])+" "+ str(loss_array[i])+"\n")
    elif(int(l_distribution) == 2 or int(l_distribution) == 16):
        for i in range(len(a_vals)):
            # for j in range(len(b_vals)):
            f.write(str(a_vals[i])+" "+str(b_vals[i]) +" "+str(loss_array[i])+"\n")
    f.close()

    # Move the loss plot to a folder
    shutil.move(loss_data, figure_data+"/"+loss_data)

# Obtained a fitted parameter through least squares
# Next we want to do gaussian process
# Name the file 
if(int(l_distribution) == 0):
    gaussian_process_name = "gp_"+element+"_"+str(fitted_alpha)
if(int(l_distribution) == 2 or int(l_distribution) == 16):
    gaussian_process_name = "gp_"+element+"_"+str(fitted_a)+"_"+str(fitted_b)
if(int(l_distribution) == 4):
    gaussian_process_name = "gp_"+element

# Name of the folder which will contain the gaussian process
os.makedirs(gaussian_process_name, exist_ok=True)

# Copy our basic input file
shutil.copyfile("basic_input", gaussian_process_name+"/basic_input")

# Move into the gaussian process directory
os.chdir(gaussian_process_name)

# Number of GP samples to run 
num_samples = 100

signal_axis = np.zeros(6)
scale_axis = np.zeros(6)
entropy_values = np.zeros((6,6))
signal_count = 0
# for signal_val in range(1,6):
# for signal_power in -2,-1,0, 1, 2, 3:
for signal_power in range(1):
    signal_axis[signal_count] = float(signal_power)
    scale_count = 0
    signal_val = 10**signal_power
    # for scale_val in range(1,6):
    # for scale_power in -2,-1,0, 1, 2, 3:
    for scale_power in range(2,3):
        scale_axis[scale_count] = float(scale_power)
        scale_val = 10**scale_power
        # We obtain all of our input files here
        if(int(l_distribution) == 0):
            # Stat distribution
            run_gaussian_process(fitted_alpha, element=element, nsamples=num_samples, l_distribution=l_distribution, signal=signal_val, scale=scale_val)
        elif(int(l_distribution) == 2 or int(l_distribution) == 16):
            # Quad distribution
            run_gaussian_process(fitted_a, fitted_b, element=element, nsamples=num_samples, l_distribution=l_distribution, signal=signal_val,scale=scale_val)
        elif(int(l_distribution) == 4):
            # Statistical distribution
            run_gaussian_process(element=element, nsamples=num_samples, l_distribution=l_distribution, signal=signal_val,scale=scale_val)

        intensities = np.zeros(num_samples)
        list_of_gp_intensities = []

        # Loop over all our samples
        for i in range(1,num_samples+1):
            name = "input_"+element+"_"+str(i)+"_"+str(signal_val)+"_"+str(scale_val)
            out_name = "output_"+element+"_"+str(i)+"_"+str(signal_val)+"_"+str(scale_val)
            nop_line = "NOP   OR   -1"
            f = open(name, "a")
            f.write(z_line+"\n")
            f.write(zs_line+"\n")
            f.write(a_line+"\n")
            f.write(be_line+"\n")
            f.write(nop_line+"\n")
            f.write(nmx_line+"\n")
            f.write(xeq_line+"\n")
            f.write(stop_line+"\n")
            f.close()

            # Run the code and parse the results, same as with the
            # parameter sweep
            os.system("../../a.out < "+name+" >"+out_name)
            gp_intensities = parse_akylas_transitions(out_name)
            list_of_gp_intensities.append(gp_intensities)


        # This is the metric used for GP testing
        kl_div_sum = 0

        # Loop over 2 to 6 as exp data doesn't go beyond this
        for j in range(2,len(exp_intens)+2):

            # We want to skip out the line if there is no experimental
            # data for it
            exp_intensity = float(exp_intens[str(j)][j-2].split("+")[0])
            if (exp_intensity < -1):
                continue

            intensity_file_name = "intensities"+str(j)+str(j-1)+".dat"

            # We will only be looking at the s transitions
            s_intensities = []

            # Loop over each sample
            for i in range(num_samples):

                # Get the correct intensity, n -> n-1
                s_intensities.append(list_of_gp_intensities[i][j][j-2])
            
            # Make a KDE plot, and store the line so we can use it
            # in a KL divergence calculation
            plot_kde = sns.kdeplot(s_intensities, warn_singular=False)
            points = plot_kde.get_lines()[0].get_data()
            x = points[0]
            y = points[1]

            # Find the experimental error
            exp_std = float(exp_intens[str(j)][j-2].split("+")[1])

            # We use this to plot the Gaussian. It needs a different range to the calculation
            # as it may be very far away
            x2 = np.linspace(exp_intensity-3*exp_std, exp_intensity+3*exp_std, 100)

            # Generate the data for the experimental Gaussian to be plotted
            experimental_gaussian = generate_experimental_gaussian(exp_intensity, exp_std, x2)

            # and we use this to calculate the entropy, as it needs to be the same
            # domain as the calculation
            experimental_gaussian_entropy = generate_experimental_gaussian(exp_intensity, exp_std, x)

            # Name of the file
            gaussian_file_name = "gaussian"+str(j)+str(j-1)+".dat"

            # Write experimental Gaussian to file
            f = open(gaussian_file_name, "w")
            f2 = open(intensity_file_name, "w")
            for i in range(len(experimental_gaussian)):
                f.write(str(x2[i])+" "+str(experimental_gaussian[i])+"\n")
            for i in range(len(y)):
                f2.write(str(x[i])+" "+str(y[i])+"\n")
            f.close()
            f2.close()

            # Move file to folder
            shutil.copyfile(intensity_file_name, "../"+figure_data+"/"+intensity_file_name)
            shutil.copyfile(gaussian_file_name, "../"+figure_data+"/"+gaussian_file_name)

            # Find the surprise of using our GP distribution instead of experimental
            print(st.entropy(experimental_gaussian_entropy, y), "ENTROPY for ",j,"->",j-1 )
            kl_div_sum += st.entropy(experimental_gaussian_entropy, y)

            # Plot the experimental gaussian, which will be on the same axis
            # as the kde plot
            plt.cla()
            # plt.plot(x2,experimental_gaussian, label="Experimental data")
            # plt.show()
            
        # Make a nice GNUplot figure. This only works if S transitions up to n = 6 are 
        # present in the experimental data
        os.chdir("../figure_data")

        # Get the correct script depending on how many experimental transitions
        # we actually have
        if (num_s_transitions >= 4):
            shutil.copyfile("../../multiplot4.gpi", "./multiplot4.gpi")
            subprocess.run(["gnuplot", "multiplot4.gpi"])
        if (num_s_transitions == 3):
            shutil.copyfile("../../multiplot3.gpi", "./multiplot3.gpi")
            subprocess.run(["gnuplot", "multiplot3.gpi"])
        if (num_s_transitions == 2):
            shutil.copyfile("../../multiplot2.gpi", "./multiplot2.gpi")
            subprocess.run(["gnuplot", "multiplot2.gpi"])
        os.chdir("../"+gaussian_process_name)

        print("Sum of entropy = ", kl_div_sum)
        entropy_values[signal_count, scale_count] = kl_div_sum
        scale_count += 1
    signal_count += 1
grid1, grid2 = np.meshgrid(signal_axis, scale_axis)
print(grid1, grid2)
print(entropy_values)
a, b = (np.argmin(entropy_values,0), np.argmin(entropy_values,1))
print(np.unravel_index(entropy_values.argmin(), entropy_values.shape))
# plt.contour(grid1, grid2, entropy_values)
# plt.colorbar()
# plt.show() 
# plt.show()
plt.xlabel('scale (10^x)', fontsize=20)
plt.ylabel('signal (10^y)', fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.imshow(entropy_values, extent=[-2,3,-2,3], origin='lower', cmap="coolwarm")
plt.colorbar().ax.tick_params(labelsize=18)
plt.show()
