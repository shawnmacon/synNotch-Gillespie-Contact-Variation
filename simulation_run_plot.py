"""
synNotch Gillespie Simulations
Shawn Macon
09.09.2025
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})

# Multi-core processing for some more speed
import multiprocessing

# Save to Database Functions
#   -   create_schema()
#   -   insert_simulation_data(parameters, gfp_data, trial_id)
#   -   insert_simulation_data(parameters, gfp_data, trial_id)
#   -   get_simulation_results(synthesis_constant, decay_constant, LoP_average, LoP_deviation_value, trial_id=None)
import database_func

# CalTech Gillespie Algorithm, modified for synNotch
#   -   variation(LoP_input, deviation)
#   -   simple_propensity(propensities, population, t, k, gamma, LoP_initial, LoP_variation)
#   -   sample_discrete(probs)
#   -   gillespie_draw(propensity_func, propensities, population, t, args=())
#   -   gillespie_ssa(propensity_func, update, population_0, time_points, args=())
import gillespie_algorithm

# Supporting Functions
#   -   start_time() 
#   -   end_time(start)
#   -   plot_distribution(mu, sigma)
#   -   progress_bar(iter_count, iter_max, start_time, bar_len=50)
#   -   plot_distribution(mu, sigma)
#   -   truncated_normal_pdf(x, mu0, sigma, a=0.0, b=1.0)
#   -   truncated_moments(mu0, sigma, a=0.0, b=1.0)
#   -   predict_G_event_timeweighted(mu0, sigma, S, D, a=0.0, b=1.0)
import supp_func

# Timing Function
from tqdm import tqdm

# Tracking Simulation Progress
import sys
import time

# ODE Solution to GFP synthesis equation, for use as needed
from scipy.integrate import odeint

from scipy.stats import norm

def ODE_trajectory(variables,t,params):
    """
    Function call for plotting ODE trajectory of GFP concentration with time

    Parameters
    ----------
    variables : array that returns the current GFP concentration
    t : array of time values for simulation of ODE trajectory
    params : array of synthesis constant and decay constant values

    Returns
    -------
    None.

    """
    G = variables[0]
    
    k = params[0]
    
    gamma = params[1]
    
    lop = params[2]
    
    dGdt = k * lop - gamma * G
    return([dGdt])

start = supp_func.start_time()

# synNotch Simulation Variation Parameters
synthesis_constant = 30.
decay_constant = 2.

# Run 5000 for final results...
trials = 5000

LoP_average = [
    0.00, 0.08, 0.16, 0.24,
    0.28, 0.31, 0.34, 0.37, 0.40,
    0.43, 0.46, 0.49, 0.52,
    0.60, 0.68, 0.76, 0.84, 0.92, 1.00
]

LoP_deviation = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

time_points = np.arange(0, 50+0.5, 0.5)

redraw_int = [0.0]

# Multiprocessing Core Allocation
#cores = multiprocessing.cpu_count // 2

# Configure Database, as needed
#database_func.create_schema()


# Simulation Loop over Trials, LoP Deviation, Synthesis Constant, Decay Constant, plotting...
print('Iterating through trials:')
print('')
iter_max = len(LoP_deviation)*trials*len(redraw_int)*len(LoP_average)
iter_count = 0
start_time = time.time()



MSE_fig = 10000001
ABC_fig = 10000002


# NOTE: Need to add another nested loop -- iterate over redraw intervals to look at data that way too...
for LoP_deviation_val in LoP_deviation:
    mse_vals = []
    abc_vals = []
    for LoP_average_val in LoP_average:
        for redraw in redraw_int:
            trials_data = np.zeros([trials, len(time_points)])
            for trial in range(0, trials):
                
                # Printing a progress bar when running simulation
                iter_count += 1
                supp_func.progress_bar(iter_count, iter_max, start_time)
                
                
                #figure_number = int((100000 * synthesis_constant + 1000 * decay_constant) + 100 * LoP_deviation_val)
                
                # Predicted equilibrium GFP based on differential equation. dG/dt = 0, S/D*LoP = GFP_eq
                equilibrium_GFP = LoP_average_val*synthesis_constant/decay_constant
                
                # Set initial concentrations for Gillespie Algorithm
                initial_concentrations = np.array([0,0], dtype = float)
                
                #Run Gillespie Simulation Algorithm
                params = (synthesis_constant, decay_constant, LoP_average_val, LoP_deviation_val, redraw)
                dataNP = gillespie_algorithm.gillespie_ssa(gillespie_algorithm.simple_propensity, gillespie_algorithm.simple_update, initial_concentrations, time_points, args = params)
                data = dataNP[:,0].tolist()
                trials_data[trial] = data
                
                
                #Retrieve only used data, and group with time points
                gfp_data = list(zip(time_points, data))
                
                # Saving to database
                """
                #Convert numpy.float64 to float for database import
                syn = round(float(synthesis_constant),2)
                dec = round(float(decay_constant),2)
                try:
                    database_func.insert_simulation_data((syn, dec, round(LoP_average,2), round(LoP_deviation_value,2), round(redraw_int, 2)), gfp_data, trial)
                    #print("Success")
                except Exception as e:
                    print(f"Error inserting data for trial {trial}: {e}")
                """
                
                #Increase countcheck for clearing RAM occasionally...
                #countcheck+=1
                
                #Clear GFP_all from RAM every 1000 trials
                #if countcheck%1000 == 0:
                #    GFP_all.clear()
                
                
                
            #plt.figure(int((100000 * synthesis_constant + 1000 * decay_constant) + 100 * LoP_deviation_value), dpi = 600)
            
            #ODE Comparison Portion
            #params = (synthesis_constant, decay_constant, LoP_average_val, LoP_deviation_value)
            #y = odeint(ODE_trajectory,initial_concentrations[0],time_points,args=(params,))
            
            
            #Plot range of trials -- only 100 of all 10,000 -- not easy to visualize with that many ::: len(trials)
            #for i in range(30):
            #    if(i == 1):
            #        plt.plot(time_points, trials_data[i], color = 'darkgray', label = r'$G_i(t)$')
            #    else:
            #        plt.plot(time_points, trials_data[i], color = 'darkgray')
            #plt.xlabel(fr" time(s) ")
            #plt.ylabel("GFP Concentration")
            #plt.title(f'Averaged GFP Synthesis With Time ')
            # X-label: :: S = {round(synthesis_constant, 1)}, D = {round(decay_constant, 1)}, $\mu$ = {round(LoP_average, 1)}, $\sigma$ = {round(LoP_deviation_value, 1)}
            # title: Redraw Interval: {redraw_int}
            # Equilibrium GFP predictions
            #print(fr'Figure {figure_number}: S = {round(synthesis_constant, 1)}, D = {round(decay_constant, 1)}, $\mu$ = {round(LoP_average_val, 1)}, $\sigma$ = {round(LoP_deviation_value, 1)}, Redraw interval = {redraw_int} s')
            G_eq = synthesis_constant/decay_constant * LoP_average_val
            #G_eq1 = (synthesis_constant/decay_constant)*(2*LoP_average_val - np.sqrt(LoP_average_val**2 + LoP_deviation_value**2))
            
            #alph = (0 - LoP_average_val)/(LoP_deviation_value)
            #beta = (1 - LoP_average_val)/(LoP_deviation_value)
            #mu_T = LoP_average_val + ((norm.pdf(alph) - norm.pdf(beta))/(norm.cdf(beta) - norm.cdf(alph)))*LoP_deviation_value
            #sigma_T_sq = LoP_deviation_value**2 * (1 + (alph*norm.pdf(alph) - beta*norm.pdf(beta))/(norm.cdf(beta) - norm.cdf(alph)) - ((norm.pdf(alph) - norm.pdf(beta))/(norm.cdf(beta) - norm.cdf(alph)))**2 )
            
            #G_eq2 = (synthesis_constant/decay_constant)*(2*mu_T - np.sqrt(mu_T**2 + sigma_T_sq))
            
            #G_eq3 = supp_func.predict_G_event_timeweighted(LoP_average_val, LoP_deviation_value, synthesis_constant, decay_constant)
            
            #G_eq4 = supp_func.predict_G_timeblend_selftau(LoP_average_val, LoP_deviation_value, synthesis_constant, decay_constant, redraw_int, beta = 0.01) if redraw_int != 0.0 else None
            #labels = ['Simple', 'Perturbed', 'Truncated Perturbed', 'Event-Redraw', 'Time-Weighted']
            
            # Mean across trials
            GFP_mean = np.mean(trials_data, axis=0)
                    
            # Mean squared error
            mse_val = np.mean((GFP_mean - G_eq)**2)
            mse_vals.append(mse_val)
            
            # Area Between Curves Error
            abc_val = np.trapz(np.abs(GFP_mean - G_eq), time_points)
            abc_vals.append(abc_val)
                
                # Plotting All Results
                #plt.plot(time_points,y[:,0], color = 'yellow', label =r'$G_{ode}(t)$')
                #plt.axhline(G_eq, color = 'tab:red', linewidth = 1, label = r'$G_{\alpha, eq} = \frac{S}{D}\mathbb{E}[L]$')
                #plt.axhline(G_eq1, color = 'tab:orange', linewidth = 1, label = r'Perturbed')
                #plt.axhline(G_eq2, color = 'tab:green', linewidth = 1, label = r'Truncated Perturbed')
                #plt.axhline(G_eq3, color = 'tab:red', linewidth = 1, label = r'Event-Redraw')
                #plt.axhline(G_eq4, color = 'tab:purple', linewidth = 1, label = r'Time-Weighted') if redraw_int != 0.0 else None
                #plt.plot(time_points, GFP_mean, color = 'k', linewidth=1, label=r'$G_{avg}(t)$')
                #plt.xlim(0, 50)
                #plt.ylim(bottom = 0)
                #plt.legend(loc='upper right')
                #plt.show()
                
                
                #title = fr'Residuals vs Simulation for S = {round(synthesis_constant, 1)}, D = {round(decay_constant, 1)}, $\mu$ = {round(LoP_average, 1)}, $\sigma$ = {round(LoP_deviation_value, 1)}'
                #supp_func.plot_residuals([G_eq, G_eq1, G_eq2, G_eq3, G_eq4], GFP_mean, labels, title, time_points)
            #plt.figure(100000000 + i)
            #plt.hist(LoP_check, bins = len(range(int(1000 * min(LoP_check)), int(1000 * max(LoP_check)))), weights=np.ones(len(LoP_check)) / len(LoP_check))
            #plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
            #plt.xlabel("L/P Value")
            #plt.ylabel("Probability")
            #plt.title(f'LoP Mean {LoP_average}, Std Dev. = {LoP_deviation_value}')
    
    plt.figure(MSE_fig, dpi=1200)
    plt.plot(LoP_average, mse_vals, label = fr'$\sigma = {LoP_deviation_val}$')
    plt.xlim(0, 1)
    plt.ylim(bottom=0)
    plt.autoscale(axis='y')
    plt.xlabel(r'$\mu$')
    plt.ylabel('MSE')
    plt.title(r'Estimation Error for Varied Redraw Frequency')
    plt.legend()

    plt.figure(ABC_fig, dpi = 1200)
    plt.plot(LoP_average, abc_vals, label = fr'$\sigma = {LoP_deviation_val}$')
    plt.xlim(0, 1)
    plt.ylim(bottom=0)
    plt.autoscale(axis='y')
    plt.xlabel(r'$\mu$')
    plt.ylabel('ABC')
    plt.title(r'Estimation Error for Varied Redraw Frequency')
    plt.legend()



plt.show()
sys.stdout.write('')
sys.stdout.flush()


supp_func.end_time(start)