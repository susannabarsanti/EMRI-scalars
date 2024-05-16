# Scalar_EMRIs 
Mathematica codes for EMRIs with scalar fields. 

The repository includes: 

1) Codes for the computations of massless scalar and gravitational fluxes for:
   
   i) equatorial circular orbits;
   
   ii) equatorial eccentric orbits;
   
   iii) circular inclined orbits.

3) Code for the computation of massive scalar field fluxes;

4) Codes for the EMRI adiabatic orbital evolution, GW template and fisher analysis.  
   

## Table of Contents  
[Fluxes](#Fluxes)  
[EMRIEvolution](#EMRIEvolution)

### Fluxes

#### Fluxes codes 

 
##### Circular
Codes for the scalar and tensorial fluxes computation that make use of the Black Hole Perturbation Toolkit https://bhptoolkit.org/mathematica-install-dev.html. Point particle in circular orbit around a Kerr black hole.  

##### Eccentric
Codes for the scalar and tensorial fluxes computation. Point particle in circular orbit around a Kerr black hole. If you make use of this codes, please cite https://arxiv.org/abs/2203.05003. 

Moreover, note that: 
- The convergence criteria have been implemented by following Viktor Skoupy: https://arxiv.org/abs/2201.07044. 
- The functions for the homogeneous solutions have been implemented by Gabriel Piovano, i.e. boundary  conditions for the radial Teukolsky equation in horizon penetrating, hyperboloidal slicing coordinates. If you make use of these boundary conditions, please acknowledge https://arxiv.org/pdf/2105.07083.  

Finally:	this notebook makes use  of the BHPToolkit https://bhptoolkit.org/mathematica-install-dev.html

(*IMPORTANT: FindRoot errors in the paclet version of the toolkit, use the GitHub one.*)

##### Inclined

##### Massive scalars 
Codes for the massive scalar fluxes computation. Point particle in circular orbit around a Kerr black hole.

If you make use of this code, please cite https://arxiv.org/abs/2212.03888.

Important Notes: 
- Toolkit vs massless : factor $\pi/4$

- THINGS TO IMPLEMENT: A FUNCTION FO R2 THAT VARIES WITH THE MASS OF THE SCALAR

#### Fluxes data

##### Circular
Some circular scalar and gravitational fluxes for circular orbit and different primary spin. 

##### Eccentric
For the paper: https://arxiv.org/abs/2203.05003
The "Scalar" folder contains the scalar fluxes for spin $ a=0.2M $ and a=0.9M.
The "Tensorial" folder contains the gravitational fluxes for spin a=0.2M and a=0.9M.

For the MCMC ANALYSIS:
"gravitational_eccentric_fluxes.dat"
Gravitational eccentric fluxes for the grid points (13x13x17=2873 points). The columns read $(a,p,e,\dot{E}_{inf}, \dot{E}_h)$
"scalar_total_level12.wdx"
Scalar eccentric fluxes for the grid points (13x13x17=2873 points).  The columns read $(a, p, e , \dot{E}_\textnormal{inf},  \dot{E}_\textnormal{h},  \dot{L}_\textnormal{inf},  \dot{L}_\textnormal{h})$.
- To use this scalar flux with our normalization, we need to divide them by a factor 4. 


### EMRI Evolution&Fishers

#### Circular
AK template with fully relativistic inspiral 

#### Equatorial eccentric
AK template with fully relativistic inspiral

#### Inclined circular
NK template with fully relativistic inspiral

#### Massive
AK template with fully relativistic circular inspiral 
