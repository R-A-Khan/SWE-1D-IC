# SWE-1D-IC
Reconstruction of missing initial data for the 1-D shallow water equations using data assimilation of observed data, with applications to tsunami modelling. Full description of model derivation and evaluation can be found at:

Kevlahan, N.KR., Khan, R. & Protas, B. On the convergence of data assimilation for the one-dimensional shallow water equations with sparse observations. Adv Comput Math 45, 3195–3216 (2019). https://doi.org/10.1007/s10444-019-09733-6

## Usage
Adjust model parameters in 
* __run_data_assimil.m__  
or
* __run_data_assimil_CG_tikhonov_multnorm.m__ (if using Conjugate Gradient descent and Tikhonov Regularisation). 
 
Utilises the data assimilation algorithm in __data_assimil.m__ or __data_assimil_CG_tikhonov_multnorm.m__.

Verification of numerical solvers used in DA scheme done using a kappa convergence test in 
* __Kappa_test_numerical.m__

Visualize results using 
* __plots.m__ 

Built using Matlab 2018a, including Optimization toolbox.
