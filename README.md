<h2 align="center">Tempered Expectation-Maximization algorithm</h2>

<h5 align="center">**Luca Brusa** &middot; **Francesco Bartolucci** &middot; **Fulvia Pennoni**</h5>

<br>

<h4>Overview</h4>

This repository contains the <tt>R</tt> functions to perform maximum likelihood estimation of the parameters of two main classes of discrete latent variable models: latent class (LC) and hidden Markov (HM) models. Both the standard Expectation-Maximization algorithm and different versions of the proposed tempered EM (T-EM) algorithm are implemented

---

<h4>Description of the content</h4>

The repository contains:

- the functions used to perform the maximum likelihood estimation of model parameters, implementing both EM and T-EM algorithms (in the subfolder *Functions*); in particular:
  - the functions <tt>est_LC.R</tt>, <tt>est_HM.R</tt> and <tt>est_HM_cont.R</tt> contains the code to run different versions of the EM algorithm on LC model, HM model with categorical responses and HM model with continuous responses respectively;
  - the remaining functions serve as auxiliary functions;
- some datasets used to test the performance of the proposed algorithm (in the subfolder *Dataset*); in particular:
  - the files <tt>simulated_sample_LC.RData</tt>, <tt>simulated_sample_HM_cat.RData</tt> and <tt>simulated_sample_HM_cont.RData</tt> are three simulated datasets, drawn from an LC model, an HM model with categorical responses and an HM model with continuous responses respectively;
  - the file <tt>anxiety_depression.RData</tt> cross-sectional contains data about measurement of anxiety and depression in oncological patients (see also the <tt>MultiLCIRT</tt> package);
  - the file <tt>criminal.RData</tt> contains longitudinal data about crimes committed by a cohort of subjects (see also the <tt>LMest</tt> package);
- a short list of examples, useful to inspect the behavior of the proposed T-EM algorithm and the relative functions, contained in the file <tt>main.R</tt>.

---

<h4>Usage</h4>

1. **Latent Class model** (function `est_LC`)

    The general formulation of the function is the following:

    ```r
    est_LC(S, yv, k, sv = NULL, tol = 1e-08, maxit = 1e+06,
           algorithm = 0, profile_pars = list(alpha = NULL, beta = NULL, ro = NULL, T0 = NULL))
    ```
    The arguments are:
    - `S`: matrix of all response sequences observed at least once in the sample and listed row-by-row
    - `yv`: vector of the frequencies of every response configuration in `S`
    - `k`: number of latent classes
    - `sv`: list of initial model parameters (piv, Piv, Phi, Psi); default is `NULL`, in which case, they are randomly drawn
    - `tol`: tolerance level for checking convergence of the algorithm as relative difference between two consecutive log-likelihood values
    - `maxit`: maximum number of iterations
    - `algorithm`: possible algorithms: 0 for standard EM algorithm (default), 1, 2, 3 for tempered EM algorithm with different tempering profiles
    - `profile_pars`: tempering constants required for the selected tempering profile
    
    <br>
  
2. **Hidden Markov model** (function `est_LM` for categorical response variables and function `est_LM_cont` for continuous response variables)
    
    The general formulation of the functions is the following:
    
    ```r
    est_HM(data, index, k, modBasic = 0, sv = NULL, tol = 1e-08, maxit = 1e+06,
           algorithm = 0, profile_pars = list(alpha = NULL, beta = NULL, ro = NULL, T0 = NULL))
    est_HM_cont(data, index, k, modBasic = 0, sv = NULL, tol = 1e-08, maxit = 1e+06,
           algorithm = 0, profile_pars = list(alpha = NULL, beta = NULL, ro = NULL, T0 = NULL))
    ```
    The arguments are:
    - `data`: dataframe in long format with id and time columns
    - `index`: character vector with two elements, the first indicating the name of the unit identifier, and the second the time occasions
    - `k`: number of latent states
    - `modBasic`: model on the transition probabilities (0 for time-heterogeneity, 1 for time-homogeneity, from 2 to (TT-1) partial time-homogeneity of a certain order)
    - `sv`: list of initial model parameters (piv, Pi, Phi if categorical responses are considered and piv, Pi, Mu, Si if continuous responses are considered); default is `NULL`, in which case, they are randomly drawn
    - `tol`: tolerance level for checking convergence of the algorithm as relative difference between two consecutive log-likelihood values
    - `maxit`: maximum number of iterations
    - `algorithm`: possible algorithms: 0 for standard EM algorithm (default), 1, 2, 3 for tempered EM algorithm with different tempering profiles
    - `profile_pars`: tempering parameters required for the selected tempering profile

---

<h4>Examples</h4>

  ```r
  load("./Dataset/simulated_sample_LC.RData")
  
  std_est <- est_LC(S = sample$S, yv = sample$yv, k = 4)
  mon_est <- est_LC(S = sample$S, yv = sample$yv, k = 4
                    algorithm = 2, profile_pars = list(alpha = 15, beta = 1))
  osc_est <- est_LC(S = sample$S, yv = sample$yv, k = 4
                    algorithm = 3, profile_pars = list(alpha = 0.3, beta = 100, ro = 5, T0 = 10))
  ```

In the example above, `std_est` denotes the results obtained with the EM algorithm, `mon_est` the results obtained with the T-EM algorithm with monotonic profile, and `osc_est` the results obtained with the T-EM algorithm with oscillating tempering profile. In particular:

- `lk` is the value of the maximized log-likelihood function
- `lkv` is the sequence of log-likelihood values at each step of the algorithm
- `it` is the number of iteration of the algorithm
- `piv`, `Piv`, `Phi` and `Psi` are the estimated parameters of the model
- `k` is the number of latent classes
- `np` is the number of parameters
- `aic` and `bic` are the values of AIC and BIC associated to the estimated model
- `algorithm` and `profile_pars` are the type of tempering profile and the corresponding tempering constants used

---

