Solution to integrating numerical 1D fluids equations as described in Problem Set 5 of the Astrophysical Fluids class taught by Prof. Eve Lee.

Answer to conceptual question in #5:
As we increase the amplitude of the perturbation, we see a shock. The shock manifests as a splitting of the gaussian perturbation with two separate peaks running to the edges of the boundaries, and because of reflective boundary conditions, they run back to the center. 

The width of the shock is set by viscosity. Though we don't explicitly set the viscosity, the second order errors to our 1D hydro solver, namely d^2 f/ dx^2, are what causes the change in viscosity and thus width of shock. 

You can balance viscosity with advection to find the width of the shock, where viscosity was v del^2 u. Then, you can find that the final width of the shock = mean free path / mach number. So a higher viscosity, higher mean free path, and more friction within the fluid, corresponds to a larger width of the shock.

Although we don’t have an explicit viscosity term as we did in the advection-diffusion equation, we have an increase in entropy in our simulations.  Entropy matters because in hydrodynamics without viscosity considerations, a shock is the only place where there will be an increase in the entropy of the gas flow. Since we must keep all quantities conserved from the Euler and continuity equations we’re integrating, we still have a shock that has a width partially set by the numerical viscosity error (alongside the bulk motion velocity we set). 

Instructions for use: 

`python ps_5_3.py` shows a real time solution to the the advection equation using the FTCS and Lax-Friedrich methods.

`python ps_5_4.py` shows a real time solution to the the advection-diffusion equation using the Lax-Friedrich methods.

`python ps_5_5.py` shows a 1D hydro solver simulation using the donor cell advection scheme. Also used reference: https://www.ita.uni-heidelberg.de/~dullemond/lectures/num_fluid_2011/Chapter_5.pdf


Package Dependencies:

matplotlib==3.3.2  
numpy==1.17.0 

