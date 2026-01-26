"""
CalTech Gillespie Algorithm synNotch Synthesis
Shawn Macon
09.09.2025
"""


import numpy as np
from scipy.stats import truncnorm

"""Simulation Algorithm"""

# Column 0 is change in m, column 1 is change in p
simple_update = np.array([[1, 0], #GFP Synthesis event 
                         [-1,0]], #GFP decay event
                         dtype=int)


"""Variation Function for Time-Domain LoP Variation"""
def variation(LoP_input, deviation):
    
    #LoP = np.random.normal(loc = LoP_input, scale = deviation)        #Note: loc = mean     scale = std dev.
    
    lower = 0.01
    upper = 1
    
    # Using scipy.stats.truncnorm to better approximate the truncated normal distribution...
    #a = (lower - LoP_input)/deviation
    #b = (upper - LoP_input)/deviation
    #LoP = truncnorm.rvs(a, b, loc = LoP_input, scale = deviation)
    
    while True:
        draw = np.random.normal(LoP_input, deviation)
        if lower <= draw <= upper:
            LoP = draw
            break
    
    return LoP


def simple_propensity(propensities, population, t, k, gamma, LoP_initial, LoP_variation, redraw_interval=None):
    """Updates an array of propensities given a set of parameters
    and an array of populations.
    """
    # Unpack population
    current_G = population[0]
    
    #print(t)
    
    if not hasattr(simple_propensity, "next_redraw"):
        simple_propensity.next_redraw = redraw_interval
        simple_propensity.LoP_draw = LoP_initial
        
    # Redraw LoP on fixed interval
    if redraw_interval is not None and t>=simple_propensity.next_redraw:
        simple_propensity.LoP_draw = variation(LoP_initial, LoP_variation)
        simple_propensity.next_redraw += redraw_interval
        #print(f"Redraw: {simple_propensity.LoP_draw}")
    
    if redraw_interval is None:
        simple_propensity.LoP_draw = LoP_initial
    # Redraw at every next reaction
    #simple_propensity.LoP_draw = variation(LoP_initial, LoP_variation)
    
    # No redraw of LoP
    #LoP_draw = LoP_initial
    
    # Clamp LoP to a range between 0.01 (effectively zero) and 1.0 -- Removed for truncnorm()!!
    #LoP_used = max(0.01, min(1.0, simple_propensity.LoP_draw))
    LoP_used = simple_propensity.LoP_draw
    
        #Appending LoP Variation to array of varied LoP values based on random draw, total count, and current number of draws to pick some LoP draws.
    #if (tick < 100):
    #    if(LoP_variation != 0):
    #        if(random.random() >= 0.5):
    #            if(len(LoP_check) < 1000000):
    #                LoP_check.append(LoP_used)
    #                #print(f"Added {LoP_used}")
            
    
    # Update propensities
    propensities[0] = k * LoP_used             # Make GFP
    propensities[1] = gamma * current_G        # Degrade GFP
    


def sample_discrete(probs):
    """Randomly sample an index with probability given by probs."""
    # Generate random number
    q = np.random.rand()
    
    # Find index
    i = 0
    p_sum = 0.0
    while p_sum < q:
        p_sum += probs[i]
        i += 1
    return i - 1



def gillespie_draw(propensity_func, propensities, population, t, args=()):
    """
    Draws a reaction and the time it took to do that reaction.
    
    Parameters
    ----------
    propensity_func : function
        Function with call signature propensity_func(population, t, *args)
        used for computing propensities. This function must return
        an array of propensities.
    population : ndarray
        Current population of particles
    t : float
        Value of the current time.
    args : tuple, default ()
        Arguments to be passed to `propensity_func`.
        
    Returns
    -------
    rxn : int
        Index of reaction that occured.
    time : float
        Time it took for the reaction to occur.
    """
    # Compute propensities
    propensity_func(propensities, population, t, *args)
    
    # Sum of propensities
    props_sum = propensities.sum()
    
    # Compute next time
    time = np.random.exponential(1.0 / props_sum)
    
    # Compute discrete probabilities of each reaction
    rxn_probs = propensities / props_sum
    
    # Draw reaction from this distribution
    rxn = sample_discrete(rxn_probs)
    
    return rxn, time



def gillespie_ssa(propensity_func, update, population_0, time_points, args=()):
    """
    Uses the Gillespie stochastic simulation algorithm to sample
    from probability distribution of particle counts over time.
    
    Parameters
    ----------
    propensity_func : function
        Function of the form f(params, t, population) that takes the current
        population of particle counts and return an array of propensities
        for each reaction.
    update : ndarray, shape (num_reactions, num_chemical_species)
        Entry i, j gives the change in particle counts of species j
        for chemical reaction i.
    population_0 : array_like, shape (num_chemical_species)
        Array of initial populations of all chemical species.
    time_points : array_like, shape (num_time_points,)
        Array of points in time for which to sample the probability
        distribution.
    args : tuple, default ()
        The set of parameters to be passed to propensity_func.        

    Returns
    -------
    sample : ndarray, shape (num_time_points, num_chemical_species)
        Entry i, j is the count of chemical species j at time
        time_points[i].
    """

    # Initialize output
    pop_out = np.empty((len(time_points), update.shape[1]), dtype=int)
    
    # Note: this is concatenating onto end identifiers for synthesis constant, decay, etc. -- not needed... can include as a check as needed...
    # args = [synthesis_constant, decay_constant, LoP_draw, LoP_deviation_value]
    
    #pop_out = np.empty((len(time_points) + 4, update.shape[1]), dtype=int)
    #pop_out[len(time_points)] = (int((100000 * args[0] + 1000 * args[1]) + 100 * args[3]))
    #pop_out[len(time_points) + 1] = 100 * args[0]
    #pop_out[len(time_points) + 2] = 100 * args[1]
    #pop_out[len(time_points) + 3] = 100 * args[3]
    #pop_out[len(time_points) + 3] = 100 * args[2]

    # Initialize and perform simulation
    i_time = 1
    i = 0
    t = time_points[0]
    population = population_0.copy()
    pop_out[0,:] = population
    propensities = np.zeros(update.shape[0])

    # unpack arguments
    k, gamma, LoP_initial, LoP_variation, redraw_interval = args

    # --- Reset LoP attributes for this run ---
    propensity_func.next_redraw = redraw_interval
    propensity_func.LoP_draw = LoP_initial
    # -----------------------------------------

    while i < len(time_points):
        while t < time_points[i_time]:
            event, dt = gillespie_draw(
                propensity_func,
                propensities,
                population,
                t,
                args=(k, gamma, LoP_initial, LoP_variation, redraw_interval)
            )
            population_previous = population.copy()
            population += update[event,:]
            t += dt

        i = np.searchsorted(time_points > t, True)
        pop_out[i_time:min(i,len(time_points))] = population_previous
        i_time = i

    return pop_out