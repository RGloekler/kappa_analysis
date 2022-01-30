# Ryan Gloekler
# ECSyD - Dr. Mahdi Nikdast
# Last updated: 1/3/22

# approximation method adapted from:
# https://www.engr.colostate.edu/~mnikdast/files/papers/Mahdi_J22.pdf
import os, sys
import numpy as np
import pandas as pd
import scipy.special as sp
import matplotlib.pyplot as plt

# define pi
pi = np.pi

# parse the input file and obtain coupling parameters
def get_parameters():
    data = open('kappa_input.txt', 'r')
    input_list, param_list = [line for line in data if line], [] # non-blank lines

    for val in input_list:
        if not val == '\n': # skip the empty lines
            val = float(val.split(':')[1]) # get each numerical value
            param_list.append(val)
    return param_list

# different curvature functions for kappa calculation
def b_approx(gamma, R, w):
    return np.sqrt(2*pi*gamma*(R + w/2))

def b_exact(gamma, R, w):
    x = gamma*(R + w/2)
    return (pi * x * np.exp(-x)) * (sp.i1(x) + sp.modstruve(-1,x))

# plot the computed kappa values vs. simulated data
def plot_sweep(min, max, parameters, Ae, Ao, Ge, Go):
    width_string = str(int(parameters[3]))

    # get the generated kappa values
    kappas = []
    gaps   = []
    for gap in range(min, max + 50, 50):
        kappas.append(compute_kappa(parameters[0], gap, parameters[2] * 1000,
        parameters[3], Ae, Ao, Ge, Go))
        gaps.append(gap)

    t_vals = [np.sqrt(1-kappa**2) for kappa in kappas]

    # handle the simulated kappa/through-port values
    kappa_data = open('data/Si_Drop_Results_' + width_string + '.txt', 'r')
    through_data = open('data/Si_Through_Results_' + width_string + '.txt', 'r')

    kappas_sim, through_sim = [], []
    gaps_sim = [x for x in range(50, 450 + 50, 50)]

    # fill lists from local files
    for line in kappa_data: kappas_sim.append(float(line))
    for line in through_data: through_sim.append(float(line))

    fig, ax = plt.subplots()

    # plot the generated kappa/through values
    k1 = ax.plot(gaps, kappas, label='kappa')
    t1 = ax.plot(gaps, t_vals, label='through')

    # plot the simulated kappa/through data
    k2 = ax.plot(gaps_sim, kappas_sim, '--', label='kappa_sim')
    t2 = ax.plot(gaps_sim, through_sim, '--', label='through_sim')

    # compute and print the average error between FDTD and model...
    err_val = compute_error(kappas[1:], kappas_sim[1:])
    print('\nKappa Percent Error: ', err_val)
    err_val = compute_error(t_vals[1:], through_sim[1:])
    print('Through Port Percent Error: ', err_val)

    # plot details
    plt.title('Kappa vs. Waveguide Gap [nm] - W = ' + width_string + 'nm', fontweight='bold')
    plt.ylabel('Kappa', fontweight='bold')
    plt.xlabel('Waveguide Gap [nm]', fontweight='bold')

    plt.grid(color='k', alpha = 0.1, linestyle='-', linewidth=2)
    plt.yticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1])
    plt.xticks(gaps)
    ax.set_xlim(50, 450)
    ax.legend(loc='upper right')
    plt.show()

# solve for the numerical approximation of kappa
def compute_kappa(wavelength, d, R, w, Ae, Ao, Ge, Go):
    return (pi / wavelength) * (Ae/Ge * np.exp(-Ge*d) * b_approx(Ge, R, w)
    + Ao/Go * np.exp(-Go*d) * b_approx(Go, R, w))

def compute_error(actual, pred):
    actual, pred = np.array(actual), np.array(pred)
    return np.mean(np.abs((pred - actual) / actual)) * 100

# main function for computing the kappa for a given ring structure
def main():
    # define the variables from the input file
    parameters = get_parameters()
    wavelength, d, R, w = parameters[0], parameters[1], parameters[2] * 1000, parameters[3]

    csv_data = pd.read_csv('data/kappa_coefficients.csv')
    coeffs = list(csv_data[str(int(parameters[3])) + 'nm'])
    Ae, Ao, Ge, Go = coeffs

    kappa   = compute_kappa(wavelength, d, R, w, Ae, Ao, Ge, Go)
    through = np.sqrt(1-kappa**2)
    print("Kappa (" + str(int(d)) + "nm gap): ", kappa)
    print("Through-Port ("+ str(int(d)) + "nm gap): ", through)

    # handle optional program argument
    if len(sys.argv) >= 2:
        if sys.argv[1] == '-plotR10': plot_sweep(50, 450, parameters, Ae, Ao, Ge, Go)

if __name__ == "__main__":
    main()
