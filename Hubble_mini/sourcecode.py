# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 20:13:21 2019

@author: Haomin(frank) Zhang
"""
#### import packages####
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
import glob as gb 
from scipy.optimize import curve_fit


good_data=[]
bad_data=[]
observation_num=[]
####identify good data ####
file_list=gb.glob('Halpha_spectral_data/Halpha_spectral_data*.csv')
for i in file_list:
    with open(i,'r') as file:
        line1=file.readline()#start inspecting datas
    good="Good"
    line1_split=line1.split(',')#define line_1_split for extracting the observation number.
    if good in line1:
        good_data.append(i)#pick the good datas
        observation_num.append(int((line1_split[2]).split()[1]))#split line 1, then chose element2, then with in the 'observation：xxxxx',chose'xxxxx'.
print(good_data)
print(observation_num)
print(line1_split[2].split()[1])
#### we exclude data NO.2,5,12,13,22,24 as they are bad.
##Print out assigned variables to check if the code works.

        
#%%
#### Play around with good data####
##straight and Gaussian###
def fit_func(x,a,mu,sig,m,c):
    gaussian=a*sp.exp(-(x-mu)**2/(2*sig**2))
    linear=m*x+c
    return gaussian+linear
## guess until converging to the answer##
#%%
def guess(x,y):
    fit_1=sp.polyfit(x,y,1)
    fit_values=sp.poly1d(fit_1)
    grad=fit_1[0]
    y_int=fit_1[1]
    residuals=y - fit_values(x)
    max_res=max(residuals)
    max_res=float(max_res)
    index=sp.where(residuals==max_res)
    u=y[index[0]]
    u=float(u)
    v=x[index[0]]
    v=float(v)
    std = sp.std(y-((grad)*x+y_int))
    return[u,v,std,grad,y_int]
#%%
wavelengths=[]
for i in good_data:
    x,y=sp.loadtxt(i,skiprows=2,delimiter=',',usecols=range(2),unpack=True)
    initial_guess=guess(x,y)
    po,po_cov=sp.optimize.curve_fit(fit_func,x,y,initial_guess)
    po=sp.array(po)
    plt.plot(x,y)
    plt.plot(x,fit_func(x,po[0],po[1],po[2],po[3],po[4]))
    plt.show()
    print('The wavelengths is:')
    print("wavelengths=%.1f+/- %.1f"%(po[1],sp.sqrt(po_cov[1,1]))+'nm')
    wavelengths.append((po[1]))
print(wavelengths)



#%%Set up function to calculate Velocity
def vel(x):
    E=656.28
    c=2.9979*(10**8)
    v=((c*(c**2-E**2))/(x**2+E**2))
    return v
    
#%% calculate V for each corresponding Wavelength
Velocity=[]
for i in range(0,24):
    velocity=vel(wavelengths[i])/1000#km/s
    Velocity.insert(i,velocity)
print('Velocity in km/s')
print(Velocity)

#%%# set up for loops to match the datas form Distance to the datas of Observation
distances=[]
data=sp.loadtxt('Distance_Mpc.csv',skiprows=1,delimiter=',',unpack=True)
for i in range(0,24):
    for j in range(0,30):
        if observation_num[i]==data[0,j]:
            distances.append(data[1,j])
print(distances)

#%%#Now plot the graph to find Hubble Constant

h,cov_h = sp.polyfit(distances,Velocity,1,cov=True)
fit=sp.poly1d(h)
hub_uncertain= sp.sqrt(cov_h[0,0])
plt.grid()
plt.ylabel('Redshift_V: km/s')
plt.xlabel('Distances: Mpc')
plt.title('Redshift Velocity against distance')
plt.plot(distances,fit(distances))
plt.plot(distances,Velocity,'+')
plt.show()
print('Hubble;s Constant:')
print("%.1f km/s/Mpc +/- %.1f km/s/Mpc" % (h[0],hub_uncertain))




    

        



    
    
    

   
    
    

    
        
    
  



    