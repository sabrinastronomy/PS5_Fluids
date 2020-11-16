"""
Physics 643 Fluids
@author: Sabrina Berger
Problem Set #5
Due November 12th, 2020
"""

# importing packages
import numpy as np
import matplotlib.pyplot as plt

# Problem 5

def find_j(f, u, ngrid):
    """
    Function to find j array as described in section 3.1.1 in handout
    :param f: current f array
    :param u: velocity
    :param ngrid: grid size
    :return: j array
    """
    # create j array
    j = np.zeros(ngrid-1)
    # indices where u is greater than 0
    indices_gt = np.where(u > 0)[0]
    # assigning u_j*f_j value to j indices where u > 0 (
    j[indices_gt] = u[indices_gt] * f[indices_gt]

    # indices where u is less than 0
    indices_lt = np.where(u < 0)[0]
    # assigning u_j*f_(j+1) value to j indices where u > 0
    j[indices_lt] = u[indices_lt] * f[indices_lt+1]
    return j


def update_f(f, j, ngrid, prefactor):
    """
    Function to update f with no source term.
    :param f: current f array
    :param j: current j array
    :param ngrid: grid size
    :param prefactor: prefactor which is dt/dx in our case
    :return: Equation (25) in handout
    """
    return f[1:ngrid -1] - prefactor * (j[1:ngrid - 1] - j[:ngrid - 2])


def update_f_wsource(f_2, f_1, ngrid, prefactor_source):
    """
    Function to update f_2 with a source term.
    :param f_2: rho*u
    :param f_1: rho
    :param ngrid: grid size
    :param prefactor_source: prefactor which is dt/dx * (speed of sound)^2
    :return: Equation (28) in handout which now considers the source term
    """
    return f_2[1:ngrid-1] - (prefactor_source) * (f_1[2:] - f_1[:ngrid-2])

def gaussian(x_arr, pos, stdev, gauss_height=0.2, offset=5):
    """
    Returns a Gaussian function evaluation at the x values inputted.
    :param x_arr: input x values to calculate gaussian at
    :param pos: mean of gaussian (\mu)
    :param stdev: standard deviation of gaussian (\sigma)
    :param gauss_height: setting height of gaussian to scale it
    :param offset: offsetting Gaussian by some value to avoid values being too small.
    :return:
    """
    return gauss_height*(np.exp(-((x_arr-pos)/stdev)**2) + offset)

def hydro_solver(dt, dx, ngrid, nsteps, c_s=343):
    """
    :param dt: time spacing
    :param dx: position spacing
    :param ngrid: size of 1D grid
    :param nsteps: a proxy for time
    :param c_s: speed of sound, default 343 m/s
    """
    # Prefactor of the integration without source term
    prefactor = dt/dx

    # u is an interface value (left side of leftmost cell to right side of rightmost cell)
    u_arr = np.zeros(ngrid-1)

    # Setting resolution of grid
    x_arr = np.arange(ngrid) * dx

    # Setting the initial f_2 values to 0
    # rho = f_1
    # rho u = f_2
    f_2 = np.zeros(ngrid)

    # Using the mean of the x array to set the gaussian mean
    gauss_pos = np.mean(x_arr)

    # Setting the standard deviation of the gaussian
    gauss_stdev = 20

    # Setting the initial values of f_1 = rho to a gaussian perturbation
    f_1 = gaussian(x_arr, gauss_pos, gauss_stdev)

    # Setting f_2 = rho*u values to their initial values based on f_1 and u_arr
    f_2[:-1] = u_arr*f_1[:-1]

    # Prefactor of the integration with source term
    prefactor_source = (prefactor)*c_s**2

    # Initiating plotting objects to use animation within for loop
    plt.ion()
    fig = plt.figure()
    plt.xlabel("Position")
    plt.ylabel("Density")
    plt.title("1D Hydrosolver with Gaussian Perturbation")
    plotter, = plt.plot(f_1, c='k')
    fig.canvas.draw()
    plt.pause(0.1) # Pausing plotting before beginning integration

    for i in range(nsteps):
        # Computing the velocities at the interfaces as shown in (26) of the handout
        # where the indices with 1/2 in them are taken to be the left side
        u_arr = 0.5 * ((f_2[0:ngrid-1]/f_1[0:ngrid-1]) + (f_2[1:ngrid]/f_1[1:ngrid]))

        # Determing the j arrays to be used in updating f_1 and f_2
        j_1 = find_j(f_1, u_arr, ngrid)
        j_2 = find_j(f_2, u_arr, ngrid)

        # Updating f_1 and f_2 without the source term
        f_1[1:ngrid-1] = update_f(f_1, j_1, ngrid, prefactor)
        f_2[1:ngrid-1] = update_f(f_2, j_2, ngrid, prefactor)

        # Updating f_1 and f_2 with the source term
        f_2[1:ngrid-1] = update_f_wsource(f_2, f_1, ngrid, prefactor_source)

        # Setting reflective boundary conditions as shown in 3.1.3 of handout
        f_1[0] = f_1[0] - (prefactor)*j_1[0]
        f_1[-1] = f_1[-1] + (prefactor)*j_1[-1]

        f_2[0] = f_2[0] - (prefactor)*j_2[0]
        f_2[-1] = f_2[-1] + (prefactor)*j_2[-1]
        if i % 500 == 0:
            # Updating animation plot instance with new data
            plotter.set_ydata(f_1)
            # plt.ylim(np.min(f_1), np.max(f_1))
            fig.canvas.draw()
            plt.pause(0.00001)
    plt.close()

# Running Problem 5
# Setting up grid
dt, dx = 0.0001, 0.1
ngrid = 10000 # grid size

nsteps = int(1e5) # time steps

# Calling 1D hydro solver with input parameters above
hydro_solver(dt, dx, ngrid, nsteps)
