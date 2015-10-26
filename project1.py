# -*- coding: utf-8 -*-

## Author:  Zander Mausolf
## Purpose: Monte Carlo simulation of neutrons through limestone to search for oil
## Compilation: python project1
## Inputs: none
## Output:  % of neutrons reflected back
## Notes:  Reuqires numpy, matplotlib (for plotting) 

import numpy as np
import matplotlib.pyplot as plt


####------------Variables----------------###
num_src_neut = input('Enter the number of source neutrons: ')
y=0
#density = 2.3*(1 + 0.1*y) #Units: g/cm^3
density_oil = .875 #Units g/cm^3
Avagadro = 6.0022e23
##Random number seeds##
np.random.seed(1)

oil_dist = 10 #distance the oil is from the groud Units: cm
oil_depth = 5#Units: cm
y_max = 25 #Units: cm
a=2.3 # parameters for variable density
b=0.1 # parameters for variable density
###-------Cross Sections----------------###
###Units: Barns (10^-24 cm)
##--------thermal-----------------------###
calcium_xs_total_thermal = 3.45901
calcium_xs_elastic_thermal = 3.01605
carbon_xs_total_thermal= 4.94262 
carbon_xs_elastic_thermal = 4.93925 
oxygen_xs_total_thermal = 3.96192 * 3 
oxygen_xs_elastic_thermal = 3.96173 * 3
hydrogen_xs_total_thermal = 30.4007 * 5
hydrogen_xs_elastic_thermal = 30.0687 * 5
sulfur_xs_total_thermal = 1.52499
sulfur_xs_elastic_thermal = 0.995363
#nitro_xs_total_thermal = 12.3647
#nitro_xs_elastic_thermal = 10.3596
##--------fast-----------------------###
calcium_xs_total_fast = 2.24281
calcium_xs_elastic_fast = 0.874423
carbon_xs_total_fast = 2.38502 
carbon_xs_elastic_fast = 2.37137 
oxygen_xs_total_fast = 2.79387 * 3
oxygen_xs_elastic_fast = 2.78041 * 3
hydrogen_xs_total_fast = 3.98955 * 5
hydrogen_xs_elastic_fast = 3.98952 * 5
sulfur_xs_total_fast = 2.67586
sulfur_xs_elastic_fast = 2.38327
#nitro_xs_total_fast = 1.9781
#nitro_xs_elastic_fast = 1.84764

##--------total cross sections ------###
scatter_xs_lime_thermal = calcium_xs_elastic_thermal + carbon_xs_elastic_thermal + oxygen_xs_elastic_thermal
scatter_xs_lime_fast = calcium_xs_elastic_fast + carbon_xs_elastic_fast + oxygen_xs_elastic_fast

scatter_xs_oil_thermal = carbon_xs_elastic_thermal*12 + hydrogen_xs_elastic_thermal  + sulfur_xs_total_thermal
scatter_xs_oil_fast = carbon_xs_elastic_fast*12 + hydrogen_xs_elastic_fast  + sulfur_xs_total_fast

total_xs_lime_thermal = calcium_xs_total_thermal + carbon_xs_total_thermal + oxygen_xs_total_thermal 
total_xs_lime_fast = calcium_xs_total_fast + carbon_xs_total_fast + oxygen_xs_total_fast
abs_xs_lime_thermal = (calcium_xs_total_thermal - calcium_xs_elastic_thermal) + (carbon_xs_total_thermal-carbon_xs_elastic_thermal) + (oxygen_xs_total_thermal -oxygen_xs_elastic_thermal)
abs_xs_lime_fast = (calcium_xs_total_fast - calcium_xs_elastic_fast) + (carbon_xs_total_fast - carbon_xs_elastic_fast) + (oxygen_xs_total_fast-oxygen_xs_elastic_fast)

total_xs_oil_thermal = carbon_xs_total_thermal*12 + hydrogen_xs_total_thermal + sulfur_xs_total_thermal
total_xs_oil_fast = carbon_xs_total_fast*12 + hydrogen_xs_total_fast + sulfur_xs_total_fast
abs_xs_oil_thermal = (carbon_xs_total_thermal*12-carbon_xs_elastic_thermal*12) + (hydrogen_xs_total_thermal - hydrogen_xs_elastic_thermal) + (sulfur_xs_total_thermal- sulfur_xs_elastic_thermal)
abs_xs_oil_fast = (carbon_xs_total_fast*12 - carbon_xs_elastic_fast*12) + (hydrogen_xs_total_fast - hydrogen_xs_elastic_fast)  + (sulfur_xs_total_fast - sulfur_xs_elastic_fast)

###-----------From PNNL---------------###
### Oil, crude (Mexican)
#Element         weight          Atom        Atom
#                fraction        fraction    density
#H 1001  1000    0.104039        0.587752    0.060606
#C 6000  6000    0.853733        0.404749    0.041736
#S 16000 16000   0.042228        0.007499    0.000773
H_frac = 0.587752
C_frac = 0.404749
S_frac = 0.007499

#                             Weight   Atom     Atom
#Element Neutron              Fraction Fraction Density
#C 6000 6000                  0.120003 0.200002 0.016847
#O 8016 8000                  0.479554 0.599991 0.050541
#Ca 20000 20000               0.400443 0.200007 0.016848
#Total                        1.000000 1.000000 0.084236
C_lime_frac = 0.200002
O_frac = 0.599991
Ca_frac = 0.200007

##Avg Number of collions to thermalize##
calcium_collision = 375.2984654	
carbon_collision = 116.8557575	
oxygen_collision = 153.5684244
hydrogen_collision = 18.42386871
sulfur_collision = 301.502774	
#nitro_collision = 135.2254046

limestone_collision = (calcium_collision*Ca_frac + carbon_collision*C_lime_frac + oxygen_collision*O_frac*3)
oil_collision = (carbon_collision*C_frac + hydrogen_collision*H_frac + sulfur_collision*S_frac)


##Mass Numbers##
calcium_mass	= 40.078	
carbon_mass	= 12.0107
oxygen_mass	= 15.9994
hydrogen_mass   = 1.00794
sulfur_mass	= 32.065
#nitro_mass   = 14.0067

limestone_mass = 100.0869 #Units: g/mol [src: PNNL]
oil_mass = 72.093  #Units: g/mol [src: PNNL]			

##arrays for plotting
dist_array = []
abs_array = []
leaked_array = []
scat_array = []
###--------------------------------------###
tally={'absorbed':0,'leaked':0,'transmitted':0}
tally_group = {'fast':0,'thermal':0}
##calculate the density based on given linear density variation##
def calc_density(a,b,position):
    density_max = a*(1+b*position)
    return density_max

## Calculated macroscopic xs
def calc_macro_xs(N_A,A_mass,density,micro_xs):
    barn = 1e-24
    macro_xs = ((N_A*density)/A_mass)*(micro_xs*barn)
    return macro_xs

##create neutron object
class Neutron(object):
    def __init__(self,direction=1,group='fast',distance=0,sigma_a=1,sigma_s=1,result='scattered',mass=1,collisions=1,material = 'limestone'):
        self.direction = direction
        self.group = group
        self.distance = distance
        self.sigma_a = sigma_a
        self.sigma_s = sigma_s
        self.result = result
        self.mass = mass
        self.collisions = collisions
        self.material = material
    def absorbed(self):
        self.result = "absorbed"
        abs_array.append(float(self.distance))
    def scattered(self):
        self.result = "scattered"
        scat_array.append(self.distance)      
    def leaked(self):
        self.result = "leaked"
        leaked_array.append(self.distance)
    def transmitted(self):
        self.result = "transmitted"
        
    def down_scatter(self):
       rand_down = np.random.rand()
       if rand_down < (1/self.collisions) : 
            self.group = 'thermal'
       else:
            self.group = 'fast'

    def scatter_limestone(self):
        self.scattered()
        flag = 'False'
        previous_pt = self.distance

        while flag == 'False':
            if self.material == 'limestone':
                density_max = calc_density(a,b,y_max)
                macro_xs_max = calc_macro_xs(Avagadro,self.mass,density_max,self.sigma_s)
                ##--get new path length--##
                path_length = np.divide(-1,macro_xs_max) * np.log(np.random.rand(1))
                ##Sample new angle
                self.direction = np.random.uniform(-1,1)
                ##--Calculate new point--##
                new_pt = previous_pt - path_length*self.direction
                ##--Find second random number to determine if path length should be accepted##     
                rand_2 = np.random.rand()
                den_at_pt = calc_density(a,b,new_pt)            
                xs_at_pt = calc_macro_xs(Avagadro,self.mass,den_at_pt,self.sigma_s)
                print 'LIMESTONE den_at_pt',den_at_pt
                print ' LIMESTONE xs_at_pt',xs_at_pt
                print ' LIMESTONE new pt', new_pt
                print ' LIMESTONE max xs', macro_xs_max
                
                if rand_2 < np.divide(xs_at_pt,macro_xs_max):
                        ##pt is accpeted
                        self.distance = previous_pt + new_pt 
                        flag = 'True'
                        return

    def scatter_oil(self):
        self.scattered()
        previous_pt = self.distance
        density_max = density_oil
            ## Find max macroscopic xs for oil## assuming constant
        macro_xs_max = calc_macro_xs(Avagadro,self.mass,density_max,self.sigma_s)
        print 'OIL previous pt',previous_pt
        print 'OIL density',density_oil
        print 'OIL XS', macro_xs_max
        path_length = np.divide(-1,macro_xs_max) * np.log(np.random.rand(1))
                ##Sample new angle
        self.direction = np.random.uniform(-1,1)
                ##--Calculate new point--##
        self.distance =  previous_pt + path_length*self.direction
                                
          
####-------in limestone------####
    def in_lime_stone(self):
        abs_rand = np.random.rand(1)
        self.mass = limestone_mass
        self.collisions = limestone_collision
        self.material = 'limestone'
        self.down_scatter()
        
        #check group
        if self.group == 'fast':
            xs_comp = abs_xs_lime_fast/total_xs_lime_fast
            self.sigma_s = scatter_xs_lime_fast
        else:
            xs_comp = abs_xs_lime_thermal/total_xs_lime_thermal
            self.sigma_s = scatter_xs_lime_thermal
        #test if absorbed in limestone
        if abs_rand < xs_comp:
            self.absorbed()
            return
        else:
            #scattering
            self.scatter_limestone()
            if self.distance < 0 :
                #neutron reflected back to surface
                self.leaked()
                return
            if oil_dist < self.distance < oil_dist +5:
                #neutron in oil
                self.in_oil()
                return
            if self.distance < oil_dist:
                self.in_lime_stone()
                return
            if self.distance > 25 :
                #neutron reflected past point of interest           
                self.transmitted()
                return
        return
###--------in Oil---------###
    def in_oil(self):
        abs_rand = np.random.rand(1)
        self.mass = oil_mass
        self.collisions = oil_collision
        self.material = 'oil'
        self.down_scatter()
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
            return
        else:
            #scattering
            self.scatter_oil()
            if self.distance < 0 :
                #neutron reflected back to surface
                self.leaked()
                return
            if self.distance < oil_dist:
                self.in_lime_stone()
                return
            if oil_dist < self.distance < oil_dist + oil_depth:
                self.in_oil
                return
            if self.distance > y_max:
                self.transmitted()
                return
        return       
##-----Goes through the logic to see what interaction happens, how neutron transports through !!-----##
    def transport(self):
        #determine if in oil or limestone
        if self.distance < oil_dist or self.distance > oil_dist + oil_depth:  
            #in limestone
            self.in_lime_stone()
            return
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

       # dist_array.append(n)
        tally[n.result] += 1
        tally_group[n.group] += 1

run()

arraytest=[]
for i in range(0,10):
    arraytest.append(i)
hist_scat,bins = np.histogram(scat_array,25)
#hist_scat,bins = np.histogram(scat_array,100)
#hist_abs, bins = np.histogram(dist_array,100)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2

plt.hist(scat_array, bins=100)
#plt.bar(center, hist_scat, align='center', width=width)
#plt.bar(center, hist_abs, align='center', width=width)

plt.show()

print tally
print tally_group