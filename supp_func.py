"""
Supporting Functions
Shawn Macon
09.09.2025
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import truncnorm


# Simulation Time Elapsed Functions
import timeit

def start_time():
    # Start a timer for evaluating total simulation time
    start = timeit.default_timer()
    return start


def end_time(start):
    # Print total time elapsed for simulations, function of start time
    end = timeit.default_timer()
    print("")
    print("")
    if(end-start > 86400):
        days = int((end-start) / 86400)
        hours = int(((end-start) - (86400 * days)) / 3600)
        minutes = int(((end-start) - (86400 * days) - (3600 * hours)) / 60)
        seconds = ((end-start) - (86400 * days) - (3600 * hours) - (60 * minutes))
        print(f"Time Elapsed: {days} days {hours} hours {minutes} minutes {round(seconds, 3)} seconds ({round(end-start, 3)} seconds)")
    elif(end-start > 3600):
        hours = int((end-start) / 3600)
        minutes = int(((end-start) - (3600 * hours)) / 60)
        seconds = ((end-start) - (3600 * hours) - (60 * minutes))
        print(f"Time Elapsed: {hours} hours {minutes} minutes {round(seconds, 3)} seconds ({round(end-start, 3)} seconds)")
    elif(end-start > 60):
        minutes = int((end-start) / 60)
        seconds = ((end-start) - (60 * minutes))
        print(f"Time Elapsed: {minutes} minutes {round(seconds, 3)} seconds ({round(end-start, 3)} seconds)")
    else:
        print(f"Time Elapsed: {round(end-start, 3)} seconds")



# Progress Bar of Simulation Function

import sys
import time

def progress_bar(iter_count, iter_max, start_time, bar_len=50):
    """
    Prints a progress bar with ETA.

    Parameters
    ----------
    iter_count : int
        Current iteration number.
    iter_max : int
        Total number of iterations.
    start_time : float
        Time when the loop started (time.time()).
    bar_len : int
        Length of the progress bar.
    """
    percent = (iter_count / iter_max) * 100
    elapsed = time.time() - start_time
    
    if iter_count > 0:
        time_per_iter = elapsed / iter_count
        remaining_iters = iter_max - iter_count
        eta = remaining_iters * time_per_iter
    else:
        eta = 0

    eta_str = time.strftime("%H:%M:%S", time.gmtime(eta))
    filled_len = int(bar_len * percent / 100)
    bar = '#' * filled_len + ' ' * (bar_len - filled_len)

    progress_string = (
        f"\rProgress: {iter_count}/{iter_max} [{percent:.2f}%] "
        f"[{bar}] ETA: {eta_str}"
    )
    sys.stdout.write(progress_string)
    sys.stdout.flush()


# Plotting the PDF of a Truncated Normal Distribution of mu, sigma within bounds of 0.01 - 1.0
def plot_distribution(mu, sigma):
    # Define the truncated normal: lower bound 0, upper bound 1
    a, b = (0.01 - mu) / sigma, (1 - mu) / sigma
    truncated = truncnorm(a, b, loc=mu, scale=sigma)
    
    # Draw samples
    vals = truncated.rvs(1000000)
    
    def draw_value(lower, upper, mu, sigma):
        while True:
            draw = np.random.normal(mu, sigma)
            if lower <= draw <= upper:
                return draw
            
    vals_custom = np.array([draw_value(0.01, 1, mu, sigma) for _ in range(1000000)])
    
    
    
    # Plot histogram
    plt.figure(dpi=600)
    plt.hist(vals, bins=100, density=True, alpha=0.7, color='darkred')
    #plt.hist(vals_custom, bins=100, density=True, alpha=0.7, color='darkgreen')
    plt.title(f'Contact Length Fraction (mu={mu}, sigma={sigma})')
    plt.xlabel(r'$\frac{l_{\alpha \beta}}{P_\alpha}$')  
    plt.ylabel('Density')
    plt.show()



# Functions for calculating the equilibrium GFP reduction due to event-time considerations in Gillespie

from scipy.stats import norm
from scipy.integrate import quad
from scipy.optimize import brentq

def truncated_normal_pdf(x, mu0, sigma, a=0.01, b=1.0):
    """
    PDF of a Normal(mu0, sigma^2) truncated to [a,b].
    """
    alpha = (a - mu0) / sigma
    beta = (b - mu0) / sigma
    Z = norm.cdf(beta) - norm.cdf(alpha)
    return np.where((x >= a) & (x <= b),
                    norm.pdf((x - mu0) / sigma) / (sigma * Z),
                    0.0)

def truncated_moments(mu0, sigma, a=0.01, b=1.0):
    """
    Mean and variance of a truncated Normal(mu0, sigma^2) on [a,b].
    """
    alpha = (a - mu0) / sigma
    beta = (b - mu0) / sigma
    Z = norm.cdf(beta) - norm.cdf(alpha)
    mu_T = mu0 + (norm.pdf(alpha) - norm.pdf(beta)) / Z * sigma
    sigma_T2 = sigma**2 * (
        1 + (alpha * norm.pdf(alpha) - beta * norm.pdf(beta)) / Z
        - ((norm.pdf(alpha) - norm.pdf(beta)) / Z) ** 2
    )
    return mu_T, sigma_T2

def predict_G_event_timeweighted(mu0, sigma, S, D, a=0.01, b=1.0):
    """
    Solve the time-weighted self-consistent equation for equilibrium GFP:

        D*G = S * [∫ L p(L)/(L S + D G) dL] / [∫ p(L)/(L S + D G) dL]

    where p(L) is the truncated Normal(mu0, sigma^2) on [a,b].
    """

    def p_trunc(L):
        return truncated_normal_pdf(L, mu0, sigma, a, b)

    def f(G):
        if G <= 0:
            G = 1e-12
        def integrand_A(L):
            return (L / (L*S + D*G)) * p_trunc(L)
        def integrand_B(L):
            return (1.0 / (L*S + D*G)) * p_trunc(L)
        A = quad(integrand_A, a, b, epsabs=1e-9, epsrel=1e-8, limit=200)[0]
        B = quad(integrand_B, a, b, epsabs=1e-9, epsrel=1e-8, limit=200)[0]
        return S*A/B - D*G

    # Initial guess from time-coupled ODE
    mu_T, sigma_T2 = truncated_moments(mu0, sigma, a, b)
    G_guess = (S / D) * mu_T
    lower, upper = 1e-8, max(G_guess * 10, 1.0)

    # Expand bracket until root is bracketed
    f_low, f_up = f(lower), f(upper)
    attempts = 0
    while f_low * f_up > 0 and attempts < 30:
        upper *= 2
        f_up = f(upper)
        attempts += 1

    if f_low * f_up > 0:
        # fallback: dense grid search
        grid = np.logspace(-8, 3, 200)
        fvals = [f(g) for g in grid]
        sign_changes = [(grid[i], grid[i+1])
                        for i in range(len(grid)-1)
                        if fvals[i]*fvals[i+1] < 0]
        if not sign_changes:
            raise RuntimeError("Unable to bracket root for G*. Try different parameters.")
        lower, upper = sign_changes[0]

    G_root = brentq(f, lower, upper, xtol=1e-12, rtol=1e-10, maxiter=200)
    return G_root

def predict_G_timeblend_selftau(mu0, sigma, S, D, delta_t,
                                a=0.01, b=1.0, beta=1.0, gamma=None):
    """
    Self-consistent time-blend with tau(G) = ∫ p(L)/(L*S + D*G) dL.
    Optional scaling beta and optional stretched exponent gamma (if gamma!=None).
    """
    def p_trunc(L):
        return truncated_normal_pdf(L, mu0, sigma, a, b)

    def tau_of_G(G):
        # tau(G) = ∫ p(L)/(L*S + D*G) dL
        def integrand(L):
            return p_trunc(L) / (L*S + D*G)
        val = quad(integrand, a, b, epsabs=1e-9, epsrel=1e-8, limit=300)[0]
        return max(val, 1e-12)

    def alpha_of_G(G):
        tauG = tau_of_G(G)
        x = beta * (delta_t / tauG)
        if gamma is None:
            return 1.0 - np.exp(-x)
        else:
            return 1.0 - np.exp(-(x**gamma))

    def w(L, G):
        aLG = L*S + D*G
        inv = 1.0 / max(aLG, 1e-12)
        alpha = alpha_of_G(G)
        return (1.0 - alpha)*inv + alpha*1.0

    def integrals_N_D(G):
        def numint(L):
            return L * p_trunc(L) * w(L, G)
        def denint(L):
            return p_trunc(L) * w(L, G)
        N = quad(numint, a, b, epsabs=1e-9, epsrel=1e-8, limit=300)[0]
        Dn= quad(denint, a, b, epsabs=1e-9, epsrel=1e-8, limit=300)[0]
        return N, Dn

    def residual(G):
        if G <= 0:
            G = 1e-12
        N, Dn = integrals_N_D(G)
        if Dn == 0 or not np.isfinite(N) or not np.isfinite(Dn):
            return 1e6
        return S * N / Dn - D * G

    # initial guess
    mu_T, _ = truncated_moments(mu0, sigma, a, b)
    G_guess = max((S / D) * mu_T, 1e-8)

    # bracket & solve (robustly, as before)
    lower, upper = 1e-8, max(G_guess * 10, 1.0)
    f_low, f_up = residual(lower), residual(upper)
    attempts = 0
    while (np.isnan(f_low) or np.isnan(f_up) or f_low * f_up > 0) and attempts < 60:
        upper *= 2; f_up = residual(upper); attempts += 1
    if np.isnan(f_low) or np.isnan(f_up) or f_low * f_up > 0:
        # fallback grid search
        grid = np.logspace(-12, 6, 400)
        fvals = [residual(g) for g in grid]
        sign_changes = [(grid[i], grid[i+1]) for i in range(len(grid)-1)
                        if np.isfinite(fvals[i]) and np.isfinite(fvals[i+1]) and fvals[i]*fvals[i+1] < 0]
        if not sign_changes:
            return G_guess
        lower, upper = sign_changes[0]

    G_root = brentq(residual, lower, upper, xtol=1e-12, rtol=1e-10, maxiter=300)
    return G_root


# Residuals Plotting for Result Accuracy Comparison
def plot_residuals(predictions: np.ndarray, GFP_mean: np.ndarray, 
                               labels: list, title: list, time_points: np.ndarray = None):
    """
    Plot residuals of equilibrium predictions vs simulation mean as a line over time,
    with optional shaded standard error, and print overall difference between each
    prediction and the mean across GFP_mean next to the labels.
    
    Parameters
    ----------
    predictions : np.ndarray
        1D array of scalar predictions (e.g., [G_eq1, G_eq2, G_eq3]).
    GFP_mean : np.ndarray
        1D array of mean GFP values over time (average across simulation runs).
    labels : list of str
        Names corresponding to each prediction for plotting.
    title : list of str
        1 element list of string title for plot.
    time_points : np.ndarray, optional
        1D array of time points corresponding to GFP_mean. If None, uses indices.
    """
    n_timepoints = len(GFP_mean)

    if len(predictions) != len(labels):
        raise ValueError("Number of predictions must match number of labels.")

    if time_points is None:
        time_points = np.arange(n_timepoints)
    elif len(time_points) != n_timepoints:
        raise ValueError("time_points must have the same length as GFP_mean.")

    # Compute overall mean across the simulation
    GFP_overall_mean = np.mean(GFP_mean)

    # Plot residuals for each prediction
    fig, ax = plt.subplots(figsize=(8, 5), dpi = 600)

    for pred, label in zip(predictions, labels):
        residual = pred - GFP_mean  # residual as a function of time
        # Compute difference to overall mean
        delta = np.abs(pred - GFP_overall_mean)
        label_with_delta = f"{label} (Δ={delta:.3f})"
        ax.plot(time_points, residual, label=label_with_delta)

    ax.axhline(0, color='black', linestyle='--', linewidth=1)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Residual (Prediction - Simulation)")
    ax.set_title(title)
    
    ax.legend()
    plt.show()

    # Optionally return the differences for programmatic use
    overall_differences = predictions - GFP_overall_mean
    return overall_differences


#plot_distribution(0.5, 0.1)