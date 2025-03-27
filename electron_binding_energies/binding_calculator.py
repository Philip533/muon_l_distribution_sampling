#!/usr/bin/python3
# THIS CODE WILL RETURN BINDING ENERGIES FOR Z-1 ATOM, SO IF Au IS REQUESTED
# THEN THE BINDING ENERGIES FOR Pt ARE GIVEN
import sys
import re
import numpy as np

# Binding energies taken from https://xdb.lbl.gov/Section1/Table_1-1.pdf
with open('./binding.dat') as f:
    lines = f.read().splitlines()

element = []
done = False

# Get the appropriate element
for i in range(len(lines)):
    split = []
    split = (lines[i].split())
    if (split[0] == "#"):
        continue
    if (split[2] == sys.argv[1]):
        element = lines[i-1].split()
        print(element)
        done = True
        break

# Catch incorrect name
if (done == False):
    print("Error in element symbol ", sys.argv[1], ". Program terminating")
    sys.exit()

# All elements have a K binding energy
k_energy = element[3]

# Convert to floats
element_float =  [float(i) for i in element[4:]]

# Here we have the full set of L and M energies
if(len(element[4:]) >= 8):
    l_energy = (sum(element_float[:3])/3)
    m_energy = (sum(element_float[3:8])/5)

# Here we have full L but non full M
elif(3 < len(element[4:]) < 8):
    l_energy = (sum(element_float[:3])/3)
    m_energy = (sum(element_float[3:])/len(element_float[3:]))

# Here we catch the case of full L but empty M
elif(len(element[4:]) == 3):
    l_energy = (sum(element_float[:3])/3)
    m_energy = 0

# Here we have empty M and non full L
elif(0 < len(element[4:]) < 3):
    l_energy = (sum(element_float)/len(element_float))
    m_energy = 0

# No L energies at all
elif(len(element[4:]) == 0):
    l_energy = 0
    m_energy = 0

# Atomic charge
charge = int(element[1])+1

# Now get atomic mass table
with open('./atomic_masses.dat') as f:
    lines = f.read().splitlines()

# Get the appropriate atomic mass
for i in range(len(lines)):
    split = []
    split = lines[i].split()
    if (split[1] == sys.argv[1]):
        mass = split[3].split("(")[0]
        break


# Effective charges for K, L, M shell
print("ZS    NE   ", charge-1, charge-3, charge -9)
print("A     OR   ", mass)
print("BE    NE   ", k_energy,l_energy, m_energy)

