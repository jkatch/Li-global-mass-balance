'''
Li isotope mass balance.
    
Python code to run Li isotope mass balance model for Kalderon, Katchinoff,
 et al., 2021 Nature.
'''

#importing necessary libraries
import scipy as sci 
import pandas as pd
import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
import math
from matplotlib.ticker import (MultipleLocator, ScalarFormatter)


'''Setting up initial model parameters'''
# time at model start and end in Ga
AGE_OLD = 3
AGE_YOUNG = 0.001

#time step in yrs
t_interval = 1E6

#number of resampling
num_monte = 1000

#time length of model to run
t_end = (AGE_OLD-AGE_YOUNG) *1E9

#Re-organizing time to start from 3 Ga to present
t_change = np.arange(0, (((AGE_OLD - AGE_YOUNG) * 1E9) + 1), t_interval)


# Modern Li mass balance parameters, all fluxes in mol/yr, ratios in permil
# Fluxes and per mil from Misra and Froelich, 2012; Li and West, 2014; Coogan et al, 2017; Caves et al., 2019
Li_riv = 10e9  # riverine Li flux (mol/yr)
R_Li_riv = 23 # riverine isotopic value (permil)


Li_ht = 5.2e9 # high-temperature hydrothermal Li flux (mol/yr)
R_Li_ht = 6.3 # high-temperature hydrothermal isotopic value (permil)

Li_lt = 15.2e9 * 0.9 # low-temperature hydrothermal Li flux (mol/yr); equal to 90% of total input to acheive modern mass balance
R_Li_lt = 13 # low-temperature hydrothermal isotopic value (permil)

Li_maac = 15.2e9 * 0.1 # marine authigenic clay Li flux (mol/yr); equal to 10% of total input to acheive modern mass balance
R_Li_maac = 20 # marine authigenic clay isotopic value (permil)


dLi_rock = 1.5 # isotope value of silicate rocks (permil)
cap_d = -17 #fractionation factor from primary rocks to secondary minerals (permil)
a = 4.5 # dimensionless constants used to fit curve to the data, Caves et al., 2019 Nature
b = 0.0575 # dimensionless constants used to fit curve to the data, Caves et al., 2019 Nature


#Load in data and read in data from excel
df_dLi_ocean = pd.read_excel("Lowess_lowest10_Age_d7Li.xlsx")
age_dLi_ocean = df_dLi_ocean['Age (Ga)']
dLi_ocean_record = df_dLi_ocean['dLi_0.6']

#add carbonate fractionation factor
dLi_ocean_record = dLi_ocean_record + 4


#interpolate Li isotope records
age_dLi_ocean = (AGE_OLD - age_dLi_ocean)
dLi_ocean_interp = sci.interpolate.interp1d(age_dLi_ocean, dLi_ocean_record, kind = "cubic")



#intial Li isotope value (i.e. 3 Gya Li isotope value)
dLi_ocean_now = dLi_ocean_interp(age_dLi_ocean[0])



#set up ages for interpolation and plot
age_all = np.linspace((age_dLi_ocean.iloc[-1]), age_dLi_ocean.iloc[0], int(t_end/t_interval)+1)
age_all_plot =  (AGE_OLD - age_all)
dLi_plot_real = dLi_ocean_interp(age_all)


#add upper and lower filtering bounds on Li record
upper_dLi_plot_real = dLi_plot_real + 5
lower_dLi_plot_real = dLi_plot_real - 4

#set up dataframes/arrays to capture results
success = pd.DataFrame()
Li_success_results = np.full((t_change.size, num_monte), np.nan)
K_riv_success_results = np.full((t_change.size, num_monte), np.nan)
K_R_riv_success_results = np.full((t_change.size, num_monte), np.nan)
K_R_riv_ini_success_results = np.full((t_change.size, num_monte), np.nan)
K_R_lt_success_results = np.full((t_change.size, num_monte), np.nan)
K_R_maac_success_results = np.full((t_change.size, num_monte), np.nan)
lt_contr_success_results = np.full((t_change.size, num_monte), np.nan)
maac_contr_success_results = np.full((t_change.size, num_monte), np.nan)
Li_result = np.zeros((t_change.size, num_monte))


#exponential decay function for modelling long-term outgassing history
def outgas(N, tx, tau):
    """Generate exponential decay curve."""
    return N * np.exp(-tx/tau) + 1

t2 = np.linspace(0, 3, 50)
y2 = outgas(1, t2, 0.7)

nrow = len(t2)

df_degass = np.zeros(shape=(nrow,num_monte))

for i in range(num_monte):
    df_degass[:,i] = outgas(np.random.uniform(0,5), t2, 0.7)
    
degassing_all_mean = np.mean(df_degass,1)
degassing_all_std = np.std(df_degass,1)
degass_fun = sci.interpolate.interp1d(t2, degassing_all_mean)

################################################################
#Model initialization, i.e. solve fluxes at time zero 
t = 0

'''Calculate Li isotope mass balance'''
#For loop that calculates Li isotope mass balance 1000x for each time-step    
for i in range(0, t_change.size):

   
   for j in range(0, num_monte):
              
       #generate random uniform distributions of key parameters of size num_monte
       Li_riv_monte = np.round(np.random.uniform(2E9, 50E9, num_monte), 2)
       K_R_riv_monte = np.round(np.random.uniform(0, 30, num_monte), 2)
       K_R_riv_monte_ini = np.round(np.random.uniform(0, 30, num_monte), 2)
       K_R_lt_monte = np.round(np.random.uniform(10, 25, num_monte), 2)
       K_R_maac_monte = np.round(np.random.uniform(1, 25, num_monte), 2)
       lt_contr_monte = np.round(np.random.uniform(0, 1, num_monte), 3)
       maac_contr_monte = 1-lt_contr_monte
       
       
       #Li isotope mass balance equation
       Li_burial_temp = ((Li_riv_monte[j]*(1 + (math.exp((K_R_riv_monte[j] - 1.5)/-17) - math.exp((K_R_riv_monte_ini[j] - 1.5)/-17)))) + (degass_fun(t) * Li_ht))
       
       
       Li_result[i,j] = ((Li_riv_monte[j]*(1 + (math.exp((K_R_riv_monte[j] - 1.5)/-17) - math.exp((K_R_riv_monte_ini[j] - 1.5)/-17))) * K_R_riv_monte[j] + (degass_fun(t) * Li_ht * R_Li_ht))/(Li_burial_temp)) + (K_R_lt_monte[j] * lt_contr_monte[j] + K_R_maac_monte[j] * maac_contr_monte[j])
       
       #filtering clause - Li isotope mass balance calculated above must be between upper and lower bounds
       if (Li_result[i,j] <= dLi_ocean_interp(t) + 5) and (Li_result[i,j] >= dLi_ocean_interp(t) - 4):
           Li_success_results[i,j] = Li_result[i,j] #stores successful Li seawater value
           
           #stores successful K (i.e., perturbation) values
           K_riv_success_results[i,j] = (Li_riv_monte[j]*(1 +(math.exp((K_R_riv_monte[j] - 1.5)/-17) - math.exp((K_R_riv_monte_ini[j] - 1.5)/-17))))
           K_R_riv_success_results[i,j] = K_R_riv_monte[j]
           K_R_riv_ini_success_results[i,j] = K_R_riv_monte_ini[j]
           K_R_lt_success_results[i,j] = K_R_lt_monte[j]
           K_R_maac_success_results[i,j] = K_R_maac_monte[j]
           lt_contr_success_results[i,j] = lt_contr_monte[j]
           maac_contr_success_results[i,j] = maac_contr_monte[j]
  

   t = t + (age_all[1]-age_all[0])
   
   
#replace all nan's with negative 1
Li_success_wo_nan = np.where(np.isfinite(Li_success_results), Li_success_results, -1)
K_R_riv_wo_nan = np.where(np.isfinite(K_R_riv_success_results), K_R_riv_success_results, -1)
K_R_riv_ini_wo_nan = np.where(np.isfinite(K_R_riv_ini_success_results), K_R_riv_ini_success_results, -1)
K_R_lt_wo_nan = np.where(np.isfinite(K_R_lt_success_results), K_R_lt_success_results, -1)
K_R_maac_wo_nan = np.where(np.isfinite(K_R_maac_success_results), K_R_maac_success_results, -1)
K_riv_wo_nan = np.where(np.isfinite(K_riv_success_results), K_riv_success_results, -1)
lt_contr_wo_nan = np.where(np.isfinite(lt_contr_success_results), lt_contr_success_results, -1)
maac_contr_wo_nan = np.where(np.isfinite(maac_contr_success_results), maac_contr_success_results, -1)

#replace all negative 1 with empty spaces and make a masked array
Li_success_wo = np.ma.masked_equal(Li_success_wo_nan,-1)
K_R_riv_wo = np.ma.masked_equal(K_R_riv_wo_nan,-1)
K_R_riv_ini_wo = np.ma.masked_equal(K_R_riv_ini_wo_nan,-1)
K_R_lt_wo = np.ma.masked_equal(K_R_lt_wo_nan,-1)
K_R_maac_wo = np.ma.masked_equal(K_R_maac_wo_nan,-1)
K_riv_wo = np.ma.masked_equal(K_riv_wo_nan,-1)
lt_contr_wo = np.ma.masked_equal(lt_contr_wo_nan, -1)
maac_contr_wo = np.ma.masked_equal(maac_contr_wo_nan, -1)

#Flattens array into 1D
Li_results_ravel = Li_success_wo.flatten()
K_R_riv_results_ravel = K_R_riv_wo.flatten()
K_R_riv_ini_results_ravel = K_R_riv_ini_wo.flatten()
K_R_lt_results_ravel = K_R_lt_wo.flatten()
K_R_maac_results_ravel = K_R_maac_wo.flatten()
K_riv_results_ravel = K_riv_wo.flatten()
lt_contr_success_ravel = lt_contr_wo.flatten()
maac_contr_success_ravel = maac_contr_wo.flatten()

########################################################
'''Parameters for plot'''

#Initialize array with first set of ages for plotting
#This first array has the age 3 Ga
age_test = np.full(np.shape(Li_success_wo[0,:]), 3)

#initial age for for loop
z = 3
#for loop that makes age plotting array
ages_array = np.full(np.shape(Li_success_wo[0,:]), z)

for i in range (0, t_change.size -1):
#    z = z - np.round(age_all[1]-age_all[0], 2)
    z = z - (age_all[1]-age_all[0])
    temp_array = np.full(np.shape(Li_success_wo[0,:]), z)
    ages_array = np.append(ages_array, temp_array)


#plot tick and label parameters
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['xtick.labelsize'] = 8

plt.rcParams['xtick.major.size'] = 4
plt.rcParams['ytick.major.size'] = 4

plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['font.family'] = 'Palatino Linotype'



'''Heatmaps plot'''
fig1, ax1 = plt.subplots(nrows = 4, ncols = 2, sharex = 'col', figsize = (12,12))
fig1.subplots_adjust(hspace = 0.3, wspace = 0.4) #adjusts spacing between subplots


#plotting seawater Li
ax1[0,0].plot(age_all_plot, dLi_plot_real, 'indianred', ls = '--', label = 'dLi', alpha = 0.9)
ax1[0,0].plot(age_all_plot, upper_dLi_plot_real, 'indianred', alpha=0.9)
ax1[0,0].plot(age_all_plot, lower_dLi_plot_real, 'indianred', alpha=0.9)
h = ax1[0,0].hist2d(ages_array, Li_results_ravel, bins = (3000, 120), range = np.array([(0,3), (0,30)]), cmap=plt.cm.Reds_r, cmin = 1, rasterized=True)
ax1[0,0].xaxis.set_major_locator(MultipleLocator(0.5))
ax1[0,0].xaxis.set_major_formatter(ScalarFormatter())
ax1[0,0].xaxis.set_minor_locator(MultipleLocator(0.25))
ax1[0,0].set_ylabel('\u03B4$^\mathregular{7}$Li$_\mathregular{sw}$ (\u2030)', fontsize = 'large', labelpad = 5)
ax1[0,0].tick_params(axis = 'both', length = 4, width = 1, labelsize = 'large', direction = 'out', pad = 4, bottom = True, left = True)
ax1[0,0].set_facecolor('whitesmoke')
cb1 = fig1.colorbar(h[3], ax=ax1[0,0], aspect = 10, cax = fig1.add_axes([0.46,0.725,0.015,0.15]))
cb1.ax.tick_params(labelsize='large')
ax1[0,0].set_ylim(-5,30)

#plotting outgassing curve
t2_plot = np.linspace(3,0,50)
ax1[0,1].plot(t2_plot, degassing_all_mean, color = 'dimgrey')
ax1[0,1].fill_between(t2_plot,degassing_all_mean + degassing_all_std, degassing_all_mean - degassing_all_std, facecolor='dimgrey',alpha = 0.5)
ax1[0,1].set_ylabel('Outgassing (x modern)', fontsize = 'large', labelpad = 5)
ax1[0,1].tick_params(axis='y', labelsize = 'large')
ax1[0,1].yaxis.set_major_locator(MultipleLocator(1))
ax1[0,1].yaxis.set_major_formatter(ScalarFormatter())
ax1[0,1].xaxis.set_major_locator(MultipleLocator(0.5))
ax1[0,1].xaxis.set_major_formatter(ScalarFormatter())
ax1[0,1].xaxis.set_minor_locator(MultipleLocator(0.25))
ax1[0,1].set_facecolor('whitesmoke')
ax1[0,1].set_xlim(0,3)

#plotting dLi_riv
h2 = ax1[1,0].hist2d(ages_array, K_R_riv_results_ravel, bins = (3000, 120), range = np.array([(0,3), (0,30)]), cmap=plt.cm.Reds_r, cmin = 1, rasterized=True)
ax1[1,0].set_ylabel('\u03B4$^\mathregular{7}$Li$_\mathregular{riv}$ (\u2030)', fontsize = 'large', labelpad = 5)
ax1[1,0].tick_params(axis = 'both', length = 4, width = 1, labelsize = 'large', direction = 'out', pad = 4, bottom = True, left = True)
ax1[1,0].yaxis.set_major_locator(MultipleLocator(10))
ax1[1,0].yaxis.set_major_formatter(ScalarFormatter())
ax1[1,0].set_ylim(0,30)
ax1[1,0].set_facecolor('whitesmoke')
cb2 = fig1.colorbar(h2[3], ax=ax1[1,0], aspect = 10, cax = fig1.add_axes([0.46,0.52,0.015,0.15]))
cb2.ax.tick_params(labelsize='large')

#plotting dLi_lt
h3 = ax1[2,0].hist2d(ages_array, K_R_lt_results_ravel, bins = (3000, 60), range = np.array([(0,3), (10,25)]), cmap=plt.cm.Reds_r, cmin = 1, rasterized=True)
ax1[2,0].set_ylabel('\u0394$^\mathregular{7}$Li$_\mathregular{lowT}$ (\u2030)', fontsize = 'large', labelpad = 5)
ax1[2,0].tick_params(axis = 'both', length = 4, width = 1, labelsize = 'large', direction = 'out', pad = 4, bottom = True, left = True)
ax1[2,0].set_yticks([10, 17, 25])
ax1[2,0].set_ylim(10,25)
ax1[2,0].set_facecolor('whitesmoke')
cb3 = fig1.colorbar(h3[3], ax=ax1[2,0], aspect = 10, cax = fig1.add_axes([0.46,0.317,0.015,0.15]))
cb3.ax.tick_params(labelsize='large')

#plotting dLi_maac
h4 = ax1[3,0].hist2d(ages_array, K_R_maac_results_ravel, bins = (3000, 60), range = np.array([(0,3), (1,25)]), cmap=plt.cm.Reds_r, cmin = 1, rasterized=True)
ax1[3,0].set_ylabel('\u0394$^\mathregular{7}$Li$_\mathregular{maac}$ (\u2030)', fontsize = 'large', labelpad = 5)
ax1[3,0].tick_params(axis = 'both', length = 4, width = 1, labelsize = 'large', direction = 'out', pad = 4, bottom = True, left = True)
ax1[3,0].set_yticks([5., 15, 25])
ax1[3,0].set_ylim(1,25)
ax1[3,0].set_facecolor('whitesmoke')
cb4 = fig1.colorbar(h4[3], ax=ax1[3,0], aspect = 10, cax = fig1.add_axes([0.46,0.112,0.015,0.15]))
cb4.ax.tick_params(labelsize='large')

#plotting lt contr
h5 = ax1[2,1].hist2d(ages_array, lt_contr_success_ravel, bins = (3000, 50), range = np.array([(0,3), (0,1)]), cmap=plt.cm.Reds_r, cmin = 1, rasterized=True)
ax1[2,1].set_ylabel('$\mathit{f}$$_\mathregular{lowT}$', fontsize = 'large', labelpad = 5)
ax1[2,1].tick_params(axis = 'both', length = 4, width = 1, labelsize = 'large', direction = 'out', pad = 4, bottom = True, left = True)
ax1[2,1].yaxis.set_major_locator(MultipleLocator(0.5))
ax1[2,1].yaxis.set_major_formatter(ScalarFormatter())
ax1[2,1].yaxis.set_minor_locator(MultipleLocator(0.25))
ax1[2,1].set_facecolor('whitesmoke')
cb5 = fig1.colorbar(h5[3], ax=ax1[2,1], aspect = 10,  cax = fig1.add_axes([0.91,0.317,0.015,0.15]))
cb5.ax.tick_params(labelsize='large')


#plotting maac contr
h6 = ax1[3,1].hist2d(ages_array, maac_contr_success_ravel, bins = (3000, 50), range = np.array([(0,3), (0,1)]), cmap=plt.cm.Reds_r, cmin = 1, rasterized=True)
ax1[3,1].set_ylabel('$\mathit{f}$$_\mathregular{maac}$', fontsize = 'large', labelpad = 5)
ax1[3,1].tick_params(axis = 'both', length = 4, width = 1, labelsize = 'large', direction = 'out', pad = 4, bottom = True, left = True)
ax1[3,1].yaxis.set_major_locator(MultipleLocator(0.5))
ax1[3,1].yaxis.set_major_formatter(ScalarFormatter())
ax1[3,1].yaxis.set_minor_locator(MultipleLocator(0.25))
ax1[3,1].set_facecolor('whitesmoke')
cb6 = fig1.colorbar(h6[3], ax=ax1[3,1], aspect = 10, cax = fig1.add_axes([0.91,0.112,0.015,0.15]))
cb6.ax.tick_params(labelsize='large')


#plotting F_riv
h7 = ax1[1,1].hist2d(ages_array, K_riv_results_ravel, bins = (3000, 50), range = np.array([(0,3), (np.nanmin(K_riv_results_ravel)-1,np.nanmax(K_riv_results_ravel)+1)]), cmap=plt.cm.Reds_r, cmin = 1, rasterized=True)
offset = ax1[1,1].yaxis.get_major_formatter().get_offset()
ax1[1,1].set_ylabel('F$_\mathregular{riv}$ x 10$^\mathregular{10}$ (mol/yr)', fontsize = 'large', labelpad = 5)
ax1[1,1].tick_params(axis = 'both', length = 4, width = 1, labelsize = 'large', direction = 'out', pad = 4, bottom = True, left = True)
ax1[1,1].yaxis.offsetText.set_visible(False)
ax1[1,1].yaxis.set_minor_locator(MultipleLocator(0.5E10))
ax1[1,1].xaxis.set_major_locator(MultipleLocator(0.5))
ax1[1,1].xaxis.set_major_formatter(ScalarFormatter())
ax1[1,1].xaxis.set_minor_locator(MultipleLocator(0.25))
ax1[1,1].set_facecolor('whitesmoke')
cb7 = fig1.colorbar(h7[3], ax=ax1[1,1], aspect = 10, cax = fig1.add_axes([0.91,0.52,0.015,0.15]))
cb7.ax.tick_params(labelsize='large')