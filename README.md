# synNotch-Gillespie-Contact-Variation
Outlines and codes for Gillespie simulations of synNotch synthesis with time-varied contact

## The Gillespie Algorithm
The Gillespie stochastic simulation algorithm (SSA) has become a valuable means of representing statistically correct system dynamics in stochastic systems with known reaction rates. Though originally designed for 
chemical reaction kinetics, for the concepts of event-based collision synthesis and decays, the Gillespie algorithm has seen widespread applications beyond this original case. Dan Gillespie originally presented the 
algorithm concept in 1977 for simple, accurate simulations of chemical and biochemical simulations with limited computational power. A primary further application of the Gillespie SSA that arose following the original
presentation and impelmentation developed in the area of biological systems analysis, where reaction dynamics in cells can be represented as a stochastic process and then tracked via the Gillespie SSA.

The key concept behind the SSA is the consideration of the variable of interest as a discretized item (Markov jump process) rather than one with continuous variation , such that it can be incremented by +1 or -1
based on the associated rates in the governing equation(s). When looking at a reaction chamber in infinitesimal time-slices, the approximate representation of overall dynamics can be described by the occurence of
reactions one at a time. A further separation and difference in the Gillespie SSA is in the consideration of time evolution. Generally, time-domain simulation of reaction dynamics with a governing differential
equation will be done by discretizing time steps, (dt = 0.1, 0.05, etc.) and iterating through these timesteps to evolve our system quantities. By considering the dynamics on the scale of individual event-type
reactions, time then needs to be discretized randomly and scaled by the probability of a reaction taking place. The Gillespie SSA does this by implementing a random waiting time until any reaction occurs (scaled
by the propensity/likelihood of some reaction occuring at a next time) and then again implements another random sample to decide what reaction occurs. In general then, the reaction time in the Gillespie SSA can 
be represented mathematically as the following.

$$
\tau = \frac{-ln(r_1)}{R} \text{ ,\hspace{1cm}   } r_1 \sim U(0,1)
$$

$$
R = \sum_j a_j(x)
$$

Following this, simulation time is updated, $t = t + \tau$, and a reaction is chosen based on a subsequent uniform random number draw, $r_2$. The reaction index that is allowed to occur is then chosen based on this 
random draw based on comparison with the sum of reaction propensities, $R$. This next random value draw is then compared against each of the individual reaction propensities to determine what the next reaction the
system sees will be. Effectively, we will be comparing the product of this integer draw $\mu$ against the reaction propensity sum to choose what would occur, with a randomness added to the weighting that normally occurs
based on the actual rate value at this time point

$$
\sum_{j=1}^\mu R_j \geq r_2 R_0, \hspace{1cm} r_2 \sim U(0,1)
$$

The repetiton of this time-reaction update routine then allows for incrementation through time to produce statistically correct individual trajectories for stochastic processes.


## synNotch Gillespie With Time Contact Variation
Here we consider the dynamics of GFP (Green Fluorescent Protein) synthesis as occurs due to the juxtacrine signalling pathway synNotch. The goal of doing this was to better understand the dynamics of synNotch
synthesis and how continuous-time variation in contact in complex cellular environments could affect overall GFP synthesis. In the system where this was initially investigated and imaged (developing fruit fly larvae),
one can reasonably expect a large amount of contact variation due to the complex dynamics of cell growth and replication in a developing organism. This then can be reasonably expected to introduce variation in the 
contact that any given cell has with its neighbors.

synNotch signaling is directly dependent upon contact between cells, where a "signal sending" cell which posesses the synNotch ligand protein within its cell membrane presents this to a "signal receiving" cell
by direct contact along the cell membrane of these two cells. This contact along the receiving cell, which contains synNotch receptor proteins, results in the binding, transcription, and expression of GFP in 
the receiving cell. The dynamics of this synthesis has been generally represented and investigated previously, with an equation present to describe the synthesis and decay of GFP, and we seek to expand upon 
this understanding here through the consideration of variation in the synthesis term due to contact variation.

$$
\frac{dG_\alpha}{dt} = \frac{l_{\alpha \beta}}{P_{\alpha}}S - D G_\alpha
$$

$$
L_{\alpha \beta} = \frac{l_{\alpha \beta}}{P_{\alpha}}
$$

Generally then in synNotch, it can be observed that GFP will be affected by some constant synthesis rate, scaled by the fraction of a cell's perimeter which is in contact with synNotch ligand cells, with a 
degradation of GFP being proportional to the current amount of GFP present. To affect this time-variation in contact length fraction, $L_{\alpha \beta}$, we can consider the chosen value 
of the overall contact fraction as a normally distributed quantity, with a mean and standard deviation. By then adding in a means of redrawing the contact length fraction from a normal distribution centered about
a chosen mean with a chosen deviation, and implement a separate frequency at which the contact length fraction will be varied. This then allows for contact variation to be investigated and understood, with a means
of comparing different redraw frequencies and understanding how that will affect synNotch GFP output, as represented in the code contained here which has been adapted from that originally created by Justin Bois 
and Michael Elowitz, Caltech.





# References
Exact stochastic simulation of coupled chemical reactions. Daniel T. Gillespie, The Journal of Physical Chemistry 1977 81 (25), 2340-2361 DOI: 10.1021/j100540a008, https://pubs.acs.org/doi/10.1021/j100540a008

Contact area and tissue growth dynamics shape synthetic juxtacrine signaling patterns. Dawson, Jonathan E. et al. Biophysical Journal, Volume 124, Issue 1, 93 - 106 DOI: 10.1016/j.bpj.2024.11.007, https://www.cell.com/biophysj/fulltext/S0006-3495(24)00716-1

Stochastic simulation of biological circuits. Justin Bois and Michael Elowitz. http://www.be150.caltech.edu/2019/handouts/12_stochastic_simulation_all_code.html
