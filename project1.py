## Author:  Zander Mausolf
## Purpose: Monte Carlo simulation of neutrons through limestone to search for oil
## Compilation: python project1
## Inputs: none
## Output:  % of neutrons reflected back
## Notes:  Reuqires numpy, matplotlib (for plotting) 

import numpy as np
import matplotlib.pyplot as plt

####------------Variables----------------###
num_src_neut = 1000
abs_xs_lime = 0.1 #micro Units: b
abs_xs_oil = 0.2 #micro Units: b
total_xs_limestone = 0.3 #micro Units: b
total_xs_oil = 0.4 #micro Units: b
y=0
density = 2.3*(1 + 0.1*y)
intial_energy = 2.5 #Units: MeV
Avagadro = 6.002 * np.power(10,13)
A_mass_lime = 100.0869 #Units: g/mol [src: PNNL]
A_mass_oil = 50  #Units: g/mol [src: PNNL]
oil_dist = 20 #Units: meters
y_max = 25
a=2.3
b=0.1
###-------Cross Sections----------------###
##thermal##
calcium_xs_total_thermal = 13.29
calcium_xs_elastic_thermal = 8.0
carbon_xs_total_thermal= 21.08
carbon_xs_elastic_thermal = 21.035
oxygen_xs_total_thermal = 15.0037
oxygen_xs_elastic_thermal = 15.00
hydrogen_xs_total_thermal = 301.1
hydrogen_xs_total_thermal = 296.86
sulfur_xs_total_thermal = 9.58
sulfur_xs_elastic_thermal = 2.80
nitro_xs_total_thermal = 65.38
nitro_xs_elastic_thermal = 41.05
##fast##
calcium_xs_total_fast = 3.51
calcium_xs_elastic_fast = 3.43
carbon_xs_total_fast = 3.55
carbon_xs_elastic_fast = 3.55
oxygen_xs_total_fast = 3.263
oxygen_xs_elastic_fast = 3.263
hydrogen_xs_total_fast = 18.27
hydrogen_xs_elastic_fast = 18.22
sulfur_xs_total_fast = 3.23
sulfur_xs_elastic_fast = 3.19
nitro_xs_total_fast = 3.83
nitro_xs_elastic_fast = 3.70

##Avg Number of collions to thermalize##
calcium_collision = 375.2984654	
carbon_collision = 116.8557575	
oxygen_collision = 153.5684244
hydrogen_collision = 18.42386871
sulfur_collision = 301.502774	
nitro_collision = 135.2254046

###--------------------------------------###
tally={'absorbed':0,'scattered':0,'leaked':0,'transmitted':0}

##calculate the density based on given linear density variation##
def calc_density(a,b,position):
    density_max = a*(1+b*position)
    return density_max

## Calculated macroscopic xs
def calc_macro_xs(N_A,A_mass,density,micro_xs):
    macro_xs = 1
    macro_xs = (N_A*density*micro_xs)/A_mass
    return macro_xs

##create neutron object
class Neutron(object):
    def __init__(self,direction=1,group='fast',distance=0,sigma_a=1,sigma_s=1,result='scattered'):
        self.direction = direction
        self.group = group
        self.distance = distance
        self.sigma_a = sigma_a
        self.sigma_s = sigma_s
        self.result = result

    def absorbed(self):
        self.result = "absorbed"

    def scattered(self):
        self.result = "scattered"
        
    def leaked(self):
        self.result = "leaked"

    def transmitted(self):
        self.transmitted = "transmitted"

##-----Goes through the logic to see what interaction happens, how neutron transports through !!-----##
    def transport(self):
        #Determine if absorbed

        if self.distance<oil_dist:  
            #in limestone
            self.lime_stone()
        else:
            #in oil
            self.in_oil()
            
    def scatter(self):
        #counter = 0
        #total_dist=0
        A_mass=1
        micro_xs = .5
        previous_pt = self.distance
        density_max = calc_density(a,b,y_max)
        ## Find max macroscopic xs
        macro_xs_max = calc_macro_xs(Avagadro,A_mass,density_max,self.sigma_s)
        
        ##--get new path length--##
        path_length = np.divide(-1,macro_xs_max) * np.log(np.random.rand())

        ##Sample new angle
        self.direction = np.random.uniform(-1,1)

        ##--Calculate new point--##
        new_pt = previous_pt + path_length*self.direction
        ##--Find second random number to determine if path length should be accepted##
        
        rand_2 = np.random.rand()
        den_at_pt = calc_density(a,b,new_pt)
        xs_at_pt = calc_macro_xs(Avagadro,A_mass,den_at_pt,micro_xs)
        
        #print 'xs/max xs', np.divide(xs_at_pt,macro_xs_max)
        #print 'rand 2 ', rand_2
        if rand_2 < np.divide(xs_at_pt,macro_xs_max):
            ##pt is accpeted
            self.distance = previous_pt + new_pt        

    #function to determine which atom a neutron will interact with in limestone
    #Determined from atomic fraction given in PNNL data
    def which_mat_lime(self):
        xs = 0
        #pick random #
        rand = np.random.rand()
        if rand < 0.2 :
            #select carbon
            self.sigma_s = 1 #fill in cs for carbon
        else:
            if rand < 0.4:
                #select calcium
                self.sigma_s  = 2 # fill in xs for calcium
            else:
                #select Oxygen
                self.sigma_s  = 3

    #function to determine which atom a neutron will interact with in Oil
    def which_mat_oil(self):
        xs=0
        #pick random #
        rand = np.random.rand()
        if rand < 0.002815 :
        #select Sulfur
            self.sigma_s  = 1
        else:
            if rand < 0.002578:
                #select Nitrogen
                self.sigma_s  = 2
            else:
                if rand < 0.36522:
                    #select Carbon
                    self.sigma_s  = 3
                else:
                    #select Hydrogen
                    self.sigma_s  = 4

    def lime_stone(self):
        ####----in limestone----####
        abs_rand = np.random.rand()
        if abs_rand < (abs_xs_lime/total_xs_limestone):
            self.absorbed()
        else:
            #scattering
            ##Determine which element to scatter with
            self.which_mat_oil()
            self.scatter()

            if self.distance < 0 :
                #neutron reflected back to surface
                self.leaked()
            elif self.distance > 25 :
                #neutron reflected past point of interest
                self.transmitted()

###--------in Oil---------###
    def in_oil(self):
        abs_rand = np.random.rand()
        #test if  absorbed in oil
        if abs_rand < (abs_xs_oil/total_xs_oil): 
                self.absorbed()
        else:
            #scattering
            ##Determine which element to scatter with
            self.which_mat_lime()
            self.scatter()

            if self.distance < 0 :
                    #neutron reflected back to surface
                    self.leaked()
            elif self.distance > 25 :
                    #neutron reflected past point of interest
                    self.transmitted()

#example creation of a neutron
def run():
    i = 0
    for i in range(num_src_neut):
        #total_dist = 0  #
        #absorbed = 0 #flag if absorbed

        #initally travel normal to Earths Surface
        ##determine initial cross section
        n = Neutron(1,'fast',0,1,1,'scattered')
        while (n.result == 'scattered'):  
            #particle not absorbed or scattered out of system keep going
            n.transport()
        tally[n.result] +=1

run()

print tally