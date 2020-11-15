"""
Physics 643 Fluids
@author: Sabrina Berger
Problem Set #5
Due November 16th, 2020
"""

# importing packages
import numpy as np
import matplotlib.pyplot as plt

# Problem 4 Function
def advection(dt, dx, ngrid, nsteps, u):
    """
    This function solves the advection equation using the FTCS and Lax-Friedrich methods.
    :param dt: time spacing
    :param dx: position spacing
    :param ngrid: size of 1D grid
    :param nsteps: a proxy for time
    :param u: velocity
    """

    # Initializing position and f arrays
    x_arr = np.arange(ngrid) * dx
    f_n_ftcs = np.copy(x_arr) / ngrid
    f_n_lf = np.copy(x_arr) / ngrid

    # Prefactor in (8) and (11) of numerical handout
    prefactor = (u * dt) / (2 * dx)

    # Initializing interactive plots, two subplots on one panel
    plt.ion()
    fig, axes = plt.subplots(1, 2)

    axes[0].set_title("FTCS (unstable)")
    axes[1].set_title("Lax-Friedrich (stable)")

    for axe in axes:
        axe.set_xlim([0, ngrid])
        axe.set_ylim([0, 1.5])

    # Plotting original to keep on plot throughout numerical integration
    axes[0].plot(x_arr, f_n_ftcs, 'm-')
    axes[1].plot(x_arr, f_n_lf, 'm-')

    # Initiating plotting objects to use animation within for loop
    ftcs_plot, = axes[0].plot(x_arr, f_n_ftcs, 'ko', markersize=2)
    fl_plot, = axes[1].plot(x_arr, f_n_lf, 'ko', markersize=2)
    fig.canvas.draw()
    plt.pause(0.1) # Pausing plotting before beginning integration


    for i in range(nsteps):
        # FTCS (forward time centered space) method: (8) in handout
        f_n_ftcs[1:ngrid - 1] = f_n_ftcs[1:ngrid - 1] - prefactor * (f_n_ftcs[2:] - f_n_ftcs[:ngrid - 2])
        # Lax-Friedrich method: (11) in handout
        f_n_lf[1:ngrid - 1] = 0.5*(f_n_lf[2:] + f_n_lf[:ngrid - 2]) - prefactor * (f_n_lf[2:] - f_n_lf[:ngrid - 2])

        # Updating plots with current iteration
        ftcs_plot.set_ydata(f_n_ftcs)
        fl_plot.set_ydata(f_n_lf)
        fig.canvas.draw()
        plt.pause(0.0001) # Pausing plot to see evolution before moving on to next iteration

    plt.close('all') # Closing plotting instance
    return


# Running Problem 3

# Setting the grid paramters
ngrid = 40  # grid size
nsteps = 4000  # time steps

# Setting up grid spacing
dt, dx = 1, 1  # grid spacing
u = -0.1  # velocity
# Note: satisfies Courant condition because dt < dx/u

# Calling advection function to begin integration
advection(dt, dx, ngrid, nsteps, u)

