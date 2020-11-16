Solution to integrating numerical 1D fluids equations as described in Problem Set 5 of the Astrophysical Fluids class taught by Prof. Eve Lee.

**Answer to conceptual questions in #5:**
As we increase the amplitude of the perturbation, we see that the shock occurs more quickly. If we set the amplitude too high, then the integration becomes unstable and causes values of $\rho$ to go to infinity. This creates an unstable numerical solution.

The shock manifests as a splitting of the gaussian perturbation with two separate peaks running to the edges of the boundaries, and because of reflective boundary conditions, they run back to the center. 

The width of the shock is set by viscosity. In an ideal numerical simulation of this kind with no error, the shock width would be infinitesimally small. Though we don't explicitly set the viscosity, the second order errors to our 1D hydro solver from are what causes the viscosity and thus width of shock. These second order errors come about because we cannot set dx and dt to be 0. The numerical velocity is proportional to our dt and dx values since the viscosity term (if we had added it) would be \frac{dx^2}{2 dt}

You can balance viscosity with advection to find the width of the shock. Then, you can find that the final width of the shock = mean free path / mach number. So a higher viscosity, higher mean free path, and more friction within the fluid, corresponds to a larger width of the shock.

Instructions for use: 

`python ps_5_3.py` shows a real time solution to the the advection equation using the FTCS and Lax-Friedrich methods.

`python ps_5_4.py` shows a real time solution to the the advection-diffusion equation using the Lax-Friedrich methods.

`python ps_5_5.py` shows a 1D hydro solver simulation using the donor cell advection scheme. Also used reference: https://www.ita.uni-heidelberg.de/~dullemond/lectures/num_fluid_2011/Chapter_5.pdf


Package Dependencies:

matplotlib==3.3.2  
numpy==1.17.0 

