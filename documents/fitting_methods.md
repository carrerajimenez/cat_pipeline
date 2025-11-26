# Comparison of Fitting Algorithms in `lmfit`

This document outlines the different optimization algorithms available in `lmfit` that can be used within the CaT Pipeline. These methods are selected via the `fitter_method` parameter in `config.yaml`.

The choice of algorithm depends on the quality of your data (SNR), the accuracy of your initial guesses, and the computational time available.

## 1. Local Optimizers (Drop-in Replacements)

These methods search for a minimum starting from the initial parameter guesses. They are generally fast but can get stuck in "local minima" if the initial guess is poor.

### Levenberg-Marquardt (`leastsq`)
* **Type:** Gradient-based.
* **Default:** This is the default method in `lmfit` and `scipy`.
* **Pros:**
    * **Fast:** Extremely efficient for well-behaved functions like spectral lines.
    * **Native Errors:** Automatically calculates the covariance matrix and standard errors.
* **Cons:**
    * **Local Minima:** If the initial guess is off (e.g., wrong line center), it may fit the noise instead of the line.
    * **Stability:** Can be sensitive to numerical instability if parameters are not scaled well.
* **Best For:** Standard analysis where initial guesses are reasonably close.

### Nelder-Mead (`nelder`)
* **Type:** Direct Search (Simplex).
* **Pros:**
    * **Robust to Noise:** Does not calculate gradients (derivatives), so it is less likely to be confused by jagged, noisy data.
    * **Simple:** Works well on "rough" functions.
* **Cons:**
    * **Slower:** Generally slower to converge than Levenberg-Marquardt.
    * **No Native Errors:** Does *not* calculate a covariance matrix. (Note: The CaT pipeline handles this by using a continuum perturbation method to estimate errors).
* **Best For:** Noisy spectra where gradient-based methods fail.

### Least Squares / Trust Region Reflective (`least_squares`)
* **Type:** Gradient-based (Modern version of `leastsq`).
* **Pros:**
    * **Bounds:** Handles parameter boundaries (`min`/`max`) more robustly than the older `leastsq`.
* **Cons:**
    * Behavior is very similar to `leastsq`.
* **Best For:** A modern alternative to `leastsq`.

---

## 2. Global Optimizers

These methods search the *entire* parameter space within the defined bounds. They are computationally expensive but are the best choice if you do not know where the line is.

### Differential Evolution (`differential_evolution`)
* **Type:** Stochastic / Genetic Algorithm.
* **Pros:**
    * **Global Minimum:** Excellent at finding the true best fit even with terrible initial guesses.
    * **Robust:** Very unlikely to get stuck in local noise.
* **Cons:**
    * **Slow:** Can take 10x to 100x longer than local methods.
    * **Bounds Required:** **Crucial:** Every parameter must have finite `min` and `max` bounds, or it will crash.
    * **No Native Errors:** Does not calculate uncertainties (The pipeline handles this via perturbation).
* **Best For:** Cases where the line center is uncertain, or where `leastsq` consistently fits the wrong feature.

### Dual Annealing (`dual_annealing`)
* **Type:** Stochastic / Simulated Annealing.
* **Pros:**
    * Similar to Differential Evolution but often converges faster for certain types of problems.
* **Cons:**
    * Also requires finite bounds.
    * Slower than local methods.
* **Best For:** Difficult optimization landscapes.

---

## 3. Bayesian / MCMC

### MCMC (`emcee`)
* **Type:** Markov Chain Monte Carlo (Sampling).
* **Pros:**
    * **The "Gold Standard" for Errors:** It doesn't just find the best value; it maps the entire probability distribution of your parameters.
    * **Correlations:** Reveals if parameters are degenerate (e.g., if amplitude and sigma are correlated).
* **Cons:**
    * **Slowest:** Requires thousands of function evaluations.
    * **Requires Burn-in:** You must run it long enough for the chains to stabilize.
* **Best For:** Final publication-quality results where robust uncertainty estimation is the priority.

---

## Summary Table

| Method Name in Config | Type | Speed | Robustness | Native Errors? | Best Use Case |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `leastsq` | Local | ⚡ Fast | Medium | ✅ Yes | Default, high SNR data. |
| `nelder` | Local | 🚀 Med | High | ❌ No | Noisy data, rough profiles. |
| `differential_evolution`| Global | 🐢 Slow | 🛡️ Very High| ❌ No | Unknown initial parameters. |
| `emcee` | Sampling| 🐌 Slowest| 🛡️ High | ✅ Yes (Dist) | Robust Error analysis. |

---

## Not Recommended

* `brute`: Brute force grid search. Too slow for 5+ parameters.
* `ampgo`, `shgo`: Require specific complex dependencies.
* `newton`, `cg`, `bfgs`: Gradient methods often less stable for non-linear curve fitting than `leastsq`.