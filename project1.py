## Author:  Zander Mausolf
## Purpose: Monte Carlo simulation of neutrons through limestone to search for oil
## Compilation: python project1
## Inputs: none
## Output:  % of neutrons reflected back
## Notes:  Reuqires numpy, matplotlib (for plotting) 

import numpy as np
import matplotlib.pyplot as plt


#Variables
num_src_neut = 10
abs_xs_lime = 0.1 #micro Units: b
abs_xs_oil = 0.2 #micro Units: b
total_xs_limestone = 0.3 #micro Units: b
total_xs_oil = 0.4 #micro Units: b
y=0
density = 2.3*(1 + 0.1*y)
intial_energy = 2.5 #Units: MeV
N_A = 6.002 * np.power(10,13)
A_mass_lime = 100.0869 #Units: g/mol [src: PNNL]
A_mass_oil = 50  #Units: g/mol [src: PNNL]
oil_dist = 20 #Units: meters
scattered_back = 0
tot_absorb = 0
tot_absorb_ar = [] #array to hold absorbed neutrons
tot_reflected = 0
tot_reflected_ar = [] #array to hold neutrons reflected back to the Earth
trans = 0
trans_ar = [] # array to hold neutrons transmitted past pt of interest (ie 25cm)

#create neutron class
class Neutron(object):
    energy= 2.5
    direction = 1
    ther_fast = False #False = fast // True = thermal
    distance = 0
    def __init__(self,energy,direction,ther_fast,distance):
        self.energy = energy
        self.direction = direction
        self.ther_fast = ther_fast
        self.distance = distance

#function to determine which atom a neutron will interact with in limestone
#returns microscopic xs of selected material
#Determined from atomic fraction given in PNNL data
def which_mat_lime():
    xs = 0
    #pick random #
    rand = np.random.rand()
    if rand < 0.2 :
        #select carbon
        xs = 1 #fill in cs for carbon
    else:
        if rand < 0.4:
            #select calcium
            xs =2 # fill in xs for calcium
        else:
            #select Oxygen
            xs = 3
    return xs
#function to determine which atom a neutron will interact with in Oil
def which_mat_oil():
    xs=0
    #pick random #
    rand = np.random.rand()
    if rand < 0.002815 :
    #select Sulfur
        xs = 1
    else:
        if rand < 0.002578:
            #select Nitrogen
            xs = 2
        else:
            if rand < 0.36522:
                #select Carbon
                xs = 3
            else:
                #select Hydrogen
                xs = 4
            
    return xs

## Calculated macroscopic xs
def calc_macro_xs(N_A,A_mass,p,micro):
    macro = 1
    macro = (N_A*p*micro)/A_mass
    return macro
    

#example creation of a neutron

i = 0
while i < num_src_neut:
   total_dist = 0  #
   absorbed = 0 #flag if absorbed
   #initally travel normal to Earths Surface
   ##determine initial cross section
   curr_dist = 1
   n = Neutron(2.5,0,False,0)
   while (absorbed == 0):  
       #particle not absorbed keep going
       abs_rand = np.random.rand() 
       ###in Oil
       if n.distance>oil_dist: # in oil
           print 'in oil'
           #test if  absorbed in oil
           if abs_rand < (abs_xs_oil/total_xs_oil): 
                print 'absorbed in oil'
                tot_absorb = tot_absorb + 1
                tot_absorb_ar.append(n)
                absorbed = 1
           else:
               #scattering
               print 'scattering'
               ##Determine which element to scatter with
               xs = which_mat_lime()
               #calculate macroscopic xs
               macro =calc_macro_xs(N_A,A_mass_oil,p,xs)
                    #calculate distance to travel
                    #calculate angle
                    #calculate new energy
                    #add to total distance
               if n.distance < 0 :
                    #neutron reflected back to surface
                    tot_reflected = tot_reflected + 1
                    tot_reflected_ar.append(n)
                    absorbed = 1 #flag to get out of loop
               elif n.distance > 25 :
                    #neutron reflected past point of interest
                    absorbed = 1
                    trans = trans +1
                    trans_ar.append(n)
       ####in limestone
       if n.distance<oil_dist:  
            print 'in limestone'
            #check if absorbed in limestone
            if abs_rand < (abs_xs_lime/total_xs_limestone):
                print 'absorbed in limestone'
                tot_absorb = tot_absorb + 1
            else:
                #scattering
                print 'scattering in limestone'
                ##Determine which element to scatter with
                xs = which_mat_oil()
                    #calculate distance to travel
                    #calculate angle
                    #calculate new energy
                    #add to total distance
                if n.distance < 0 :
                    #neutron reflected back to surface
                    tot_reflected = tot_reflected + 1
                    tot_reflected_ar.append(n)
                    absorbed = 1 #flag to get out of loop
                elif n.distance > 25 :
                    #neutron reflected past point of interest
                    absorbed = 1
                    trans = trans +1
                    trans_ar.append(n)


        
        
    

