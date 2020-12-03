# stat_shock_patch
Implementation of the shock-simulating patch dynamics of Maclean, Bunder, Kevrekidis and Roberts.

Top-level file is script_patchSol.m. Running it will:
* initialise and compute a solution by patch dynamics
* get a trusted solution at the space and time points computed by the patch simulation.

The trusted solution is computed by a numerical quadrature via the Cole-Hopf transformation. It may produce warnings, or even be unstable, for particularly small values of t and/or *very* near to a shock with small diffusion. However it is generally stable on long times at distances only slightly displaced from the shock.
