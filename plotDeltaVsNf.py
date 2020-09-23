# -*- encoding: utf-8 -*-

## @package TBTKview
#  @file plotLDOS.py
#  @brief Plot local density of states
#
#  @author Kristofer Björnson

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes
import matplotlib.cm
import scipy.ndimage.filters
import mpl_toolkits.mplot3d
import sys
from scipy.signal import find_peaks
import os
import scipy.interpolate as interpolate

def findpeakdiscrete(input_array, peak_nr):
    threshold = 0.1
    counter = 0
    length = len(input_array)
    for i in range(int(length/2), length):
        if input_array[i] > threshold:
            counter +=1
        if counter == peak_nr:
            return i
    raise Exception("No peak found")

#if(len(sys.argv) != 2):
#	print( "Error, the following parameters are needed: .hdf5-file")
#	exit(1)
    
filenamebase = "vz_"

#results_dir = ""
#extension = "_diag_size21.hdf5"
#figure_extension =  "_diag"

results_dir = "" # "kebnekaiseResults/"
extension_sc = "_diag_size21.hdf5"
extension_normal = "_diag_size21_normalState.hdf5" 
figure_extension =  "_diagDelta0_1_normalState"

sigma = 0.2
peak_height = 7.5

delta_bulk = 0.12

delta_at_imp = []

vz_points = []

text_pos = [1,1]

interesting_vz_values = [0, 0.5, 1.0, 1.5, 1.75, 2.0]#[0, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
#interesting_vz_values = np.arange(0.0, 10.1, 0.0625)

delta_fig_x, delta_ax_x = plt.subplots()
delta_fig_xy, delta_ax_xy = plt.subplots()
ldos_normal_fig, ldos_normal_ax = plt.subplots()


for vz in np.arange(0.0, 10.1, 0.0625):
    vz_format = "{:.6f}".format(vz)
    fName_sc = filenamebase + vz_format + extension_sc
    fName_normal = filenamebase + vz_format + extension_normal
    if os.path.exists(results_dir + fName_normal) and os.path.exists(results_dir + fName_sc) and \
    vz in interesting_vz_values:
        
        
        vz_points.append(vz)

        filename_sc = results_dir + fName_sc
        filename_normal = results_dir + fName_normal
        
        file_sc = h5py.File(filename_sc, 'r');
        file_normal = h5py.File(filename_normal, 'r');
        dataset_sc = file_sc['LDOS']
        dataset_normal = file_normal['LDOS']
        
        dataset_charge_spin0 = file_normal['chargeDensity_spin_0Real']
        dataset_charge_spin1 = file_normal['chargeDensity_spin_1Real']
        dataset_charge_spin2 = file_normal['chargeDensity_spin_2Real']
        dataset_charge_spin3 = file_normal['chargeDensity_spin_3Real']
        
        vz_format += figure_extension

        
        
        data_dimensions = dataset_sc.shape
        y_plane = int(np.floor(data_dimensions[1]/2))
        physical_dimensions = len(data_dimensions) - 1 #Last dimensions are for energy.
        energy_resolution = data_dimensions[physical_dimensions];
        limits = dataset_sc.attrs['UpLowLimits']
        
        
        datasetDeltaReal = file_sc['deltaReal0']
        datasetDeltaImag = file_sc['deltaImag0']
        delta = np.array(datasetDeltaReal)
        
        
        size_x = data_dimensions[0]
        size_y = data_dimensions[1]
        
        sigma_discrete_units = sigma*energy_resolution/(limits[0] - limits[1])
        
        x = np.arange(0, data_dimensions[0], 1)
        y = np.arange(limits[1], limits[0], (limits[0] - limits[1])/energy_resolution)
        energies = y
        X, Y = np.meshgrid(x, x)

        energy_range = delta_bulk *1.0
        middle_energies = np.logical_and(energies>=-energy_range, energies<=energy_range)
        
        dos_fermi = np.zeros((size_x, size_x))
        
        for x_pos in x:
            for y_pos in x:
                dos_fermi[x_pos, y_pos] = np.sum(dataset_normal[x_pos, y_pos, :][middle_energies])/np.sum(middle_energies)
                
        dos_fermi_normal_fig, dos_fermi_normal_ax = plt.subplots()
        
        cs = dos_fermi_normal_ax.contourf(X.transpose(), Y.transpose(), dos_fermi) #, cmap=matplotlib.cm.coolwarm)
        dos_fermi_normal_ax.contour(cs, colors='k', linewidths=0.2)

        
        cbar = dos_fermi_normal_fig.colorbar(cs)
        cbar.ax.set_ylabel(r'$N_0$')
        
        dos_fermi_normal_ax.set_xlabel(r'$x$')
        dos_fermi_normal_ax.set_ylabel(r'$y$')
        dos_fermi_normal_ax.text(text_pos[0], text_pos[1], r'$V_Z =$ ' + str(vz))
        dos_fermi_normal_fig.savefig("figures/dos" + vz_format + ".png")
        
        plt.figure()
        dos_xy = []
        for i in range(len(x)):
            dos_xy.append(dos_fermi[i,i])
        plt.plot(x, dos_xy)
        plt.savefig("figures/dos_xy_vz_" + vz_format +  ".png")
        plt.close()
        
        plt.figure()
        plt.plot(x, dos_fermi[:,y_plane])
        plt.savefig("figures/dos_x_vz_" + vz_format +  ".png")
        plt.close()
        
        ratio = delta/dos_fermi
        
        ratio_fig, ratio_ax = plt.subplots()
        
        cs = ratio_ax.contourf(X.transpose(), Y.transpose(), ratio) #, cmap=matplotlib.cm.coolwarm)
        ratio_ax.contour(cs, colors='k', linewidths=0.2)

        
        cbar = ratio_fig.colorbar(cs)
        cbar.ax.set_ylabel(r'$\frac{\Delta}{N_0}$')
        
        ratio_ax.set_xlabel(r'$x$')
        ratio_ax.set_ylabel(r'$y$')
        ratio_ax.text(text_pos[0], text_pos[1], r'$V_Z =$ ' + str(vz))
        ratio_fig.savefig("figures/ratio" + vz_format + ".png")
        

        
        charge_density = np.array(dataset_charge_spin0) + np.array(dataset_charge_spin1)
                
                
        charge_fermi_normal_fig, charge_fermi_normal_ax = plt.subplots()
        
        cs = charge_fermi_normal_ax.contourf(X.transpose(), Y.transpose(), charge_density) #, cmap=matplotlib.cm.coolwarm)
        charge_fermi_normal_ax.contour(cs, colors='k', linewidths=0.2)

        
        cbar = charge_fermi_normal_fig.colorbar(cs)
        cbar.ax.set_ylabel(r'$\rho_{el}$')
        
        charge_fermi_normal_ax.set_xlabel(r'$x$')
        charge_fermi_normal_ax.set_ylabel(r'$y$')
        charge_fermi_normal_ax.text(text_pos[0], text_pos[1], r'$V_Z =$ ' + str(vz))
        charge_fermi_normal_fig.savefig("figures/charge" + vz_format + ".png")
        
        charge_fermi_normal_xy_fig, charge_fermi_xy_normal_ax = plt.subplots()
        charge_xy = []
        delta_xy = []
        for i in range(len(x)):
            charge_xy.append(charge_density[i,i])
            delta_xy.append(delta[i,i])
        charge_xy = np.array(charge_xy)
        delta_xy = np.array(delta_xy)
        charge_fermi_xy_normal_ax.plot(x, charge_xy, label=r'$\rho_{el}$')
        charge_fermi_xy_normal_ax.set_xlabel(r'$xy$')
        charge_fermi_xy_normal_ax.set_ylabel(r'$\rho_{el}$')
        delta_xy_ax = charge_fermi_xy_normal_ax.twinx()
        delta_xy_ax.plot(x, delta_xy, color="r", label=r'$\Delta$')
        delta_xy_ax.set_ylabel(r'$\Delta$')
        ratio_xy_ax = charge_fermi_xy_normal_ax.twinx()
        ratio_xy_ax.plot(x, delta_xy/charge_xy, color="k", label=r'$\frac{\Delta}{\rho_{el}}$')
        ratio_xy_ax.set_ylabel(r'$\frac{\Delta}{\rho_{el}}$')
        prod_xy_ax = charge_fermi_xy_normal_ax.twinx()
        prod_xy_ax.plot(x, delta_xy*charge_xy, color="g", label=r'$\Delta \rho_{el}$')
        prod_xy_ax.set_ylabel(r'$\Delta \rho_{el}$')
        charge_fermi_normal_xy_fig.legend()
        charge_fermi_normal_xy_fig.savefig("figures/charge_xy_vz_" + vz_format +  ".png")
        
        charge_fermi_normal_x_fig, charge_fermi_x_normal_ax = plt.subplots()
        charge_fermi_x_normal_ax.plot(x, charge_density[y_plane,:])
        charge_fermi_x_normal_ax.set_xlabel(r'$x$')
        charge_fermi_x_normal_ax.set_ylabel(r'$\rho_{el}$')
        charge_fermi_normal_x_fig.savefig("figures/charge_x_vz_" + vz_format +  ".png")
        
        ratio = delta/dos_fermi
        
        ratio_fig, ratio_ax = plt.subplots()
        
        cs = ratio_ax.contourf(X.transpose(), Y.transpose(), ratio) #, cmap=matplotlib.cm.coolwarm)
        ratio_ax.contour(cs, colors='k', linewidths=0.2)

        
        cbar = ratio_fig.colorbar(cs)
        cbar.ax.set_ylabel(r'$\frac{\Delta}{N_0}$')
        
        ratio_ax.set_xlabel(r'$x$')
        ratio_ax.set_ylabel(r'$y$')
        ratio_ax.text(text_pos[0], text_pos[1], r'$V_Z =$ ' + str(vz))
        ratio_fig.savefig("figures/ratio" + vz_format + ".png")
        
        product = delta*dos_fermi
        
        product_fig, product_ax = plt.subplots()
        
        cs = product_ax.contourf(X.transpose(), Y.transpose(), product) #, cmap=matplotlib.cm.coolwarm)
        product_ax.contour(cs, colors='k', linewidths=0.2)
        


        
        cbar = product_fig.colorbar(cs)
        cbar.ax.set_ylabel(r'$\frac{\Delta}{N_0}$')
        
        product_ax.set_xlabel(r'$x$')
        product_ax.set_ylabel(r'$y$')
        product_ax.text(text_pos[0], text_pos[1], r'$V_Z =$ ' + str(vz))
        product_fig.savefig("figures/product" + vz_format + ".png")
        
        
        ratio = delta/charge_density
        
        ratio_fig, ratio_ax = plt.subplots()
        
        cs = ratio_ax.contourf(X.transpose(), Y.transpose(), ratio) #, cmap=matplotlib.cm.coolwarm)
        ratio_ax.contour(cs, colors='k', linewidths=0.2)

        
        cbar = ratio_fig.colorbar(cs)
        cbar.ax.set_ylabel(r'$\frac{\Delta}{\rho_{el}}$')
        
        ratio_ax.set_xlabel(r'$x$')
        ratio_ax.set_ylabel(r'$y$')
        ratio_ax.text(text_pos[0], text_pos[1], r'$V_Z =$ ' + str(vz))
        ratio_fig.savefig("figures/ratio_chargedensity" + vz_format + ".png")
        
        product = delta*charge_density
        
        product_fig, product_ax = plt.subplots()
        
        cs = product_ax.contourf(X.transpose(), Y.transpose(), product) #, cmap=matplotlib.cm.coolwarm)
        product_ax.contour(cs, colors='k', linewidths=0.2)
        


        
        cbar = product_fig.colorbar(cs)
        cbar.ax.set_ylabel(r'$\Delta \rho_{el}$')
        
        product_ax.set_xlabel(r'$x$')
        product_ax.set_ylabel(r'$y$')
        product_ax.text(text_pos[0], text_pos[1], r'$V_Z =$ ' + str(vz))
        product_fig.savefig("figures/product_charge" + vz_format + ".png")
        
        
        

        delta_ax_x.plot(x, delta[:,y_plane], label=r'$V_Z = $' + str(vz) )

        
        delta_xy = []
        for i in range(len(x)):
            delta_xy.append(delta[i,i])
        delta_ax_xy.plot(x, delta_xy, label=r'$V_Z = $' + str(vz))
        
        
        Z1 = dataset_normal[y_plane, y_plane, :]
        signal = Z1[: int(data_dimensions[2]/2)]


        


        # 300 represents number of points to make between T.min and T.max
        xnew = np.linspace(y.min(), y.max(), 300000)  
        
        spline = interpolate.splrep(y, Z1, k=1)
        ldos_new = interpolate.splev(xnew, spline)
        ldos_new = scipy.ndimage.filters.gaussian_filter1d(ldos_new, sigma_discrete_units)
        ldos_normal_ax.plot(xnew, ldos_new, label=r'$V_Z = $' + str(vz))
#        plt.plot(y[peaks[-1]], Z1[peaks[-1]], 'x')
#        
        lim = np.max(Z1)*1.2
#        plt.vlines([Eg[y_plane], -Eg[y_plane]], 0,  lim)



        

delta_ax_x.set_xlabel(r'$x$')
delta_ax_x.set_ylabel(r'$\Delta$')
delta_ax_x.set_xlim([0,30])
delta_ax_x.legend(loc='lower left')
delta_fig_x.savefig("figures/delta_profile_compare_x.png")

delta_ax_xy.set_xlabel(r'$xy$')
delta_ax_xy.set_ylabel(r'$\Delta$')
delta_ax_xy.legend(loc='lower left')
delta_ax_xy.set_xlim([0,30])
delta_fig_xy.savefig("figures/delta_profile_compare_xy.png")

ldos_normal_ax.set_xlabel(r'$E$')
ldos_normal_ax.set_ylabel(r'LDOS')
ldos_normal_ax.legend(loc='upper left')
ldos_normal_ax.set_xlim([-2*delta_bulk, 2*delta_bulk])
ldos_normal_ax.set_ylim([-0.5, 10])
ldos_normal_fig.savefig("figures/LDOS_compare_normal.png")

