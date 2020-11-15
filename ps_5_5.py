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

def find_j(f, u, ngrid): # Check this
    j = np.zeros(ngrid-1)
    indices_gt = np.where(u > 0)[0]
    j[indices_gt] = u[indices_gt] * f[indices_gt]
    indices_lt = np.where(u < 0)[0]
    j[indices_lt] = u[indices_lt] * f[indices_lt+1]
    return j

def update_f(f, j, ngrid, prefactor):
    return f[1:ngrid - 1] - prefactor * (j[1:ngrid - 1] - j[:ngrid - 2])

def update_f_wsource(f_2, f_1, ngrid, prefactor_source):
    return f_2[1:ngrid-1] - (prefactor_source) * (f_1[2:] - f_1[:ngrid-2])

def gaussian(x_arr, pos, stdev, offset=0):
    return np.exp(-((x_arr-pos)/stdev)**2) + offset

def hydro_solver(dt, dx, ngrid, nsteps, c_s=343):
    prefactor = dt/dx

    # u is an interface value (left side of leftmost cell to right side of rightmost cell)
    u_arr = np.arange(ngrid-1)

    # Setting resolution of grid
    x_arr = np.arange(ngrid) * dx

    # rho = f_1
    # rho u = f_2
    f_2 = np.zeros(ngrid)

    gauss_pos = np.mean(x_arr)
    gauss_stdev = x_arr[-1]/20
    f_1 = 0.1*gaussian(x_arr, gauss_pos, gauss_stdev)
    f_2[:-1] = u_arr*f_1[:-1]
    prefactor_source = (prefactor)*c_s**2


    # plt.ion()
    # fig = plt.figure()
    # plotter, = plt.plot(f_1)
    # fig.canvas.draw()
    # plt.pause(0.1)

    for i in range(100):
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

        print(max(f_2))
        print(max(f_1))

        f_1[0] = f_1[0] - (dt/dx)*j_1[0]
        f_1[-1] = f_1[-1] + (dt/dx)*j_1[-1]

        f_2[0] = f_2[0] - (dt/dx)*j_2[0]
        f_2[-1] = f_2[-1] + (dt/dx)*j_2[-1]

        f_1[:-100] = 0
        f_2[:-100] = 0

        # if i % 1 == 0:
        #     plotter.set_ydata(f_1)
        #     plt.ylim(-2*np.max(f_1), 2*np.max(f_1))
        #     fig.canvas.draw()
        #     plt.pause(1)
    plt.close()

dt = 0.00001
dx = 0.02
ngrid = 5000
nsteps = 120000
# Running Problem 5
hydro_solver(dt, dx, ngrid, nsteps)
