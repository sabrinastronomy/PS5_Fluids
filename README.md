Solution to integrating numerical 1D fluids equations as described in Problem Set 5 of the Astrophysical Fluids class taught by Prof. Eve Lee.

Answer to conceptual question in #5:
As we increase the amplitude of the perturbation, we see a shock. The shock manifests as a splitting of the gaussian perturbation with two separate peaks running to the edges of the boundaries, and because of reflective boundary conditions, they run back to the center. 

The width of the shock is set by viscosity. Though we don't explicitly set the viscosity, the second order errors to our 1D hydro solver, namely d^2 f/ dx^2, are what causes the change in viscosity and thus width of shock. 

Instructions for use: 

`python ps_5_3.py` shows a real time solution to the the advection equation using the FTCS and Lax-Friedrich methods.

`python ps_5_4.py` shows a real time solution to the the advection-diffusion equation using the Lax-Friedrich methods.

`python ps_5_5.py` shows a 1D hydro solver simulation using the donor cell advection scheme.


Package Dependencies:

matplotlib==3.3.2  
numpy==1.17.0 

