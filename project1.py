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

y=0
density = 2.3*(1 + 0.1*y)
intial_energy = 2.5 #Units: MeV
Avagadro = 6.002 * np.power(10,13)

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
hydrogen_xs_elastic_thermal = 296.86
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

##total cross sections
scatter_xs_lime_thermal = calcium_xs_elastic_thermal + carbon_xs_elastic_thermal + oxygen_xs_elastic_thermal
scatter_xs_lime_fast = calcium_xs_elastic_fast + carbon_xs_elastic_fast + oxygen_xs_elastic_fast

scatter_xs_oil_thermal = carbon_xs_elastic_thermal + hydrogen_xs_elastic_thermal + nitro_xs_elastic_thermal + sulfur_xs_total_thermal
scatter_xs_oil_fast = carbon_xs_elastic_fast + hydrogen_xs_elastic_fast + nitro_xs_elastic_fast + sulfur_xs_total_fast

total_xs_lime_thermal = calcium_xs_total_thermal + carbon_xs_total_thermal + oxygen_xs_total_thermal 
total_xs_lime_fast = calcium_xs_total_fast + carbon_xs_total_fast + oxygen_xs_total_fast
abs_xs_lime_thermal = (calcium_xs_total_thermal - calcium_xs_elastic_thermal) + (carbon_xs_total_thermal-carbon_xs_elastic_thermal) + (oxygen_xs_total_thermal -oxygen_xs_elastic_thermal)
abs_xs_lime_fast = (calcium_xs_total_fast - calcium_xs_elastic_fast) + (carbon_xs_total_fast - carbon_xs_elastic_fast) + (oxygen_xs_total_fast-oxygen_xs_elastic_fast)

total_xs_oil_thermal = carbon_xs_total_thermal + hydrogen_xs_total_thermal + nitro_xs_total_thermal + sulfur_xs_total_thermal
total_xs_oil_fast = carbon_xs_total_fast + hydrogen_xs_total_fast + nitro_xs_total_fast + sulfur_xs_total_fast
abs_xs_oil_thermal = (carbon_xs_total_thermal-carbon_xs_elastic_thermal) + (hydrogen_xs_total_thermal - hydrogen_xs_elastic_thermal) + (nitro_xs_total_thermal - nitro_xs_elastic_thermal) + (sulfur_xs_total_thermal- sulfur_xs_elastic_thermal)
abs_xs_oil_fast = (carbon_xs_total_fast - carbon_xs_elastic_fast) + (hydrogen_xs_total_fast - hydrogen_xs_elastic_fast) + (nitro_xs_total_fast - nitro_xs_total_fast) + (sulfur_xs_total_fast - sulfur_xs_elastic_fast)


##Avg Number of collions to thermalize##
calcium_collision = 375.2984654	
carbon_collision = 116.8557575	
oxygen_collision = 153.5684244
hydrogen_collision = 18.42386871
sulfur_collision = 301.502774	
nitro_collision = 135.2254046

##Mass Numbers##
calcium_mass	= 40.078	
carbon_mass	= 12.0107
oxygen_mass	= 15.9994
hydrogen_mass   = 1.00794
sulfur_mass	= 32.065
nitro_mass   = 14.0067

limestone_mass = 100.0869 #Units: g/mol [src: PNNL]
oil_mass = 50  #Units: g/mol [src: PNNL]			

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
    def __init__(self,direction=1,group='fast',distance=0,sigma_a=1,sigma_s=1,result='scattered',mass=1,collisions=1):
        self.direction = direction
        self.group = group
        self.distance = distance
        self.sigma_a = sigma_a
        self.sigma_s = sigma_s
        self.result = result
        self.mass = mass
        self.collisions = collisions
    def absorbed(self):
        self.result = "absorbed"

    def scattered(self):
        self.result = "scattered"
        
    def leaked(self):
        self.result = "leaked"

    def transmitted(self):
        self.transmitted = "transmitted"
        
    def down_scatter(self):
       rand_down = np.random.rand()
       if rand_down < self.collisions : 
            self.group = 'thermal'
       else:
            self.group = 'fast'            
    
    def scatter(self):
        ##decide which material to scatter with ##        
        #self.which_mat_lime()
        ##decide if the scatter thermalizes neutron##

        previous_pt = self.distance
        density_max = calc_density(a,b,y_max)
        ## Find max macroscopic xs
        macro_xs_max = calc_macro_xs(Avagadro,self.mass,density_max,self.sigma_s)
        
        ##--get new path length--##
        path_length = np.divide(-1,macro_xs_max) * np.log(np.random.rand())

        ##Sample new angle
        self.direction = np.random.uniform(-1,1)

        ##--Calculate new point--##
        new_pt = previous_pt + path_length*self.direction
        ##--Find second random number to determine if path length should be accepted##
        
        rand_2 = np.random.rand()
        den_at_pt = calc_density(a,b,new_pt)
        xs_at_pt = calc_macro_xs(Avagadro,self.mass,den_at_pt,self.sigma_s)
        
        #print 'xs/max xs', np.divide(xs_at_pt,macro_xs_max)
        #print 'rand 2 ', rand_2
        if rand_2 < np.divide(xs_at_pt,macro_xs_max):
            ##pt is accpeted
            self.distance = previous_pt + new_pt        

####-------in limestone------####
    def in_lime_stone(self):
        abs_rand = np.random.rand()
        self.mass = limestone_mass
        #check group
        if self.group == 'fast':
            xs_comp = abs_xs_lime_fast /total_xs_lime_fast
            self.sigma_s = scatter_xs_lime_fast
        else:
            xs_comp = abs_xs_lime_thermal/total_xs_lime_thermal
            self.sigma_s = scatter_xs_lime_thermal
        #test if absorbed in oil
        if abs_rand < xs_comp:
            self.absorbed()
        else:
            #scattering
            self.scatter()
            if self.distance < 0 :
                #neutron reflected back to surface
                self.leaked()
            elif self.distance > oil_dist:
                #neutron in oil
                self.in_oil()
            elif self.distance > 25 :
                #neutron reflected past point of interest           
                self.transmitted()

###--------in Oil---------###
    def in_oil(self):
        abs_rand = np.random.rand()
        self.mass = oil_mass
        #check group
        if self.group == 'fast':
            xs_comp = abs_xs_oil_fast /total_xs_oil_fast
            self.sigma_s = scatter_xs_oil_fast
        else:
            xs_comp = abs_xs_oil_thermal/total_xs_oil_thermal
            self.sigma_s = scatter_xs_oil_thermal
        #test if absorbed in oil
        if abs_rand < xs_comp :
            self.absorbed()
        else:
            #scattering
            ##Determine which element to scatter with
            #self.which_mat_lime()
            self.scatter()

            if self.distance < 0 :
                #neutron reflected back to surface
                self.leaked()
            elif self.distance > 25 :
                #neutron reflected past point of interest
                self.transmitted()
                    
##-----Goes through the logic to see what interaction happens, how neutron transports through !!-----##
    def transport(self):
        #determine if in oil or limestone
        if self.distance < oil_dist:  
            #in limestone
            self.in_lime_stone()
        else:
            #in oil
            self.in_oil()
            
#example creation of a neutron
def run():
    i = 0
    for i in range(num_src_neut):
        #initally travel normal to Earths Surface
        n = Neutron()
        while (n.result == 'scattered'):  
            #particle not absorbed or scattered out of system keep going
            n.transport()
        tally[n.result] +=1

run()

print tally