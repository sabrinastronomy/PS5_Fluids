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

def implicit_diffusion_term(f_n, diff_factor):
    """
    This function calculates the diffusion term using the implicit scheme outlined in section 2.5.2 of handout.
    :param f_n:
    :param beta:
    :return:
    """
    beta = diff_factor*dt/(dx**2)
    # Initializing tri-diagonal matrix with elements -beta, (1+2beta), -beta according to handout
    a_mat = np.eye(ngrid)*(1 + 2*beta) + np.eye(ngrid, k=1)*(-beta) + np.eye(ngrid, k=-1)*(-beta)

    # Imposing no-slip boundary condition at left boundary
    a_mat[0][0] = 1
    a_mat[0][1] = 0

    # Imposing stress-free boundary condition at right boundary
    a_mat[-1][-1] = 1 + beta

    # Solving for f_(n+1) by solving linear algebra equation:
    # \vec{f_(n+1)} = A_inv \vec{f_(n)}
    f_n_1 = np.linalg.solve(a_mat, f_n)
    return f_n_1

def advection_diffusion(dt, dx, ngrid, nsteps, u, d_arr = [0.1, 1]):
    """
    This function solves the advection-diffusion equation using the Lax-Friedrich method.
    :param dt: time spacing
    :param dx: position spacing
    :param ngrid: size of 1D grid
    :param nsteps: a proxy for time
    :param u: velocity
    :param d_arr: different diffusion factors
    """

    # Initializing position and f arrays
    x_arr = np.arange(ngrid) * dx
    f_n_lf = np.copy(x_arr) / ngrid
    f_n_lf_2 = np.copy(x_arr) / ngrid

    # Prefactor in (11) of numerical handout for Lax-Friedrich method
    prefactor = (u * dt) / (2 * dx)

    # Initializing interactive plots, two subplots on one panel
    plt.ion()
    fig, axes = plt.subplots(1, 2)

    axes[0].set_title("D = {}".format(d_arr[0]))
    axes[1].set_title("D = {}".format(d_arr[1]))


    # Plotting original to keep on plot throughout numerical integration
    axes[0].plot(x_arr, f_n_lf, 'm-')
    axes[1].plot(x_arr, f_n_lf_2, 'm-')

    # Initiating plotting objects to use animation within for loop
    lf_plot, = axes[0].plot(x_arr, f_n_lf, 'ko', markersize=2)
    lf_plot2, = axes[1].plot(x_arr, f_n_lf_2, 'ko', markersize=2)
    fig.canvas.draw()
    plt.pause(0.1) # Pausing plotting before beginning integration

    for i in range(nsteps):
        # Using operator splitting
        # Updating f with diffusion term using the implicit method
        f_n_lf[1:ngrid - 1] = implicit_diffusion_term(f_n_lf, d_arr[0])[1:ngrid - 1]
        f_n_lf_2[1:ngrid - 1] = implicit_diffusion_term(f_n_lf_2, d_arr[1])[1:ngrid - 1]

        # Updating f with advection term using Lax-Friedrich method
        f_n_lf[1:ngrid - 1] = 0.5*(f_n_lf[2:] + f_n_lf[:ngrid - 2]) - prefactor * (f_n_lf[2:] - f_n_lf[:ngrid - 2])
        f_n_lf_2[1:ngrid - 1] = 0.5*(f_n_lf_2[2:] + f_n_lf_2[:ngrid - 2]) - prefactor * (f_n_lf_2[2:] - f_n_lf_2[:ngrid - 2])

        # Updating plots with current iteration
        lf_plot.set_ydata(f_n_lf)
        lf_plot2.set_ydata(f_n_lf_2)
        fig.canvas.draw()
        plt.pause(0.0001) # Pausing plot to see evolution before moving on to next iteration
    plt.close('all') # Closing plotting instance
    return


# Running Problem 4

# Setting up the grid
ngrid = 40  # grid size
nsteps = 4000  # time steps

# Setting up grid spacing
dt, dx = 1, 1  # grid spacing
u = -0.1  # velocity
# Note: satisfies Courant condition because dt < dx/u

# Calling advection-diffusion function to begin integration
advection_diffusion(dt, dx, ngrid, nsteps, u)

