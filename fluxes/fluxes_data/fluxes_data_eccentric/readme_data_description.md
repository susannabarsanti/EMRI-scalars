*eccentric_scalar.m*
Code for the scalar eccentric scalar flux computation. 

eccentric_tensorial.m
Code for the gravitational eccentric scalar flux computation. 

Interpolator.nb
It imports the grav_data.wdx and scalar_total_level12.wdx and interpolate the fluxes by making use of the Maarteen’s interpolator in ChebInt.m.

ChebInt.m
Package used for the fluxes interpolation. Given to me by Lorenzo, written by Maarten I guess. I never changed it. 

EmriEcc09.m
Code for the EMRI evolution, used for the trajectories. It needs EdotinfIntofaueSC_2_total.wdx and ChebInt.m to obtain the function of the fluxes. The trajectories were computed by using the fluxes EdotinfIntofaueSC_2_total.wdx, which was obtained with Interpolator.nb and considering scalar_level2.wdx.
I computed some trajectories with the correct fluxes (EdotinfIntofaueSC_new.wdx from scalar_total_level12.wdx) and they are slightly different (not exported). 

evol.nb
Notebook for the trajectories computation, exported as evolution_d0.1_a0.9_rp6_ra11_level2included.dat.

scalar_total_level12.wdx
Scalar fluxes for the grid points (computed by Niels&me and used by Lorenzo. 13x13x17=2873 points). 
- Hai fatto il check con il toolkit per il caso circolare e visto che quello eccentrico si avvicina con continuità a quello circolare. 
- The flux is total, in the sense that it considers the sum over the negative m too (for the circular case). 
- To use this scalar flux with our normalization, we need to divide them by a factor 4. 
- See the description in the folder EMRI_MCMC/DATA

scalar_level2.wdx
Scalar fluxes computed by me only. There are some differences with the Niels’ one. To check: the one I used to complete scalar_total_level12.wdx were different as well? Are they bad? 
