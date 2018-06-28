# PRED_PREY_IMSP_IMSP1
Simulating Predator-Prey interactions in 1D and 2D 

PRED_PREY_IMSP is a collection of simple MATLAB routines using the finite difference/finite element method for simulating the dynamics of predator-prey interactions modelled by a nonlinear reaction-diffusion system by means of Implicit-Symplectic Schemes. 

IMSP1 is a collection of MATLAB routines using the finite element / difference method for the dynamics of predator-prey interactions in 1 and 2 spatial dimension and time with the IMSP first order scheme (part of PRED_PREY_IMSP).

The MATLAB code is mostly self explanatory, with the names of variables and parameters corresponding to the symbols used in the finite element / difference methods described in the papers referenced below. Copies of the MATLAB codes are freely available via the link below.

## Reference

    Garvie M.R. , "Finite difference schemes for reaction-diffusion equations modeling predator-prey interactions in MATLAB," Bulletin of Mathematical Biology (2007) 69:931-956
    Garvie M.R., Burkardt J., Morgan J. "Simple Finite Element Methods for Approximating Predator-Prey Dynamics in Two Dimensions using MATLAB," Bulletin of Mathematical Biology (2015) 77:548-578
    Diele, F., Garvie M.R., Trenchea, C. " Numerical analysis of a first-order in time implicit-symplectic scheme for predator-prey systems," submitted 2016

## Download codes for IMSP1

Files you may copy include:

    fd1dKin2_IMSP1.m    MATLAB code for IMSP first order scheme applied to Kinetics (ii) in 1D.
    fe2d_n_fast_IMSP1.m    MATLAB code for IMSP first order scheme applied to Kinetics (i) in 2D.
    fe2d_n_fast_IMSP1_test.m    MATLAB code for running IMSP first order scheme solving the two dimensional example in Section 4.2. of Reference [3].

PRED_PREY_IMSP is distributed under the GNU GPL; see the License and Copyright notice for more information. 
