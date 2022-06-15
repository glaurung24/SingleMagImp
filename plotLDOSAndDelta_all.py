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
extension = "_diag_size21.hdf5" #"_chebychev_GPUsize20.hdf5"
figure_extension =  "_diagDelta0_1"

sigma = 0.0001 #0.002
peak_height = 7.5

delta_at_imp = []

vz_points = []

text_pos = [1,1]



for vz in np.arange(0.0, 10.1, 0.0025):
    vz_format = "{:.6f}".format(vz)
    fName = filenamebase + vz_format + extension
    if os.path.exists(results_dir + fName):
        
        vz_points.append(vz)

        filename = results_dir + fName
        
        file = h5py.File(filename, 'r');
        dataset = file['LDOS']
        
        vz_format += figure_extension
        
        
        eigenvals = file["EigenValues"][:]
        
        positive_energies = np.sort(abs(eigenvals))
        
        Eg = positive_energies[2]
        
        
        data_dimensions = dataset.shape
        y_plane = int(np.floor(data_dimensions[1]/2))
        physical_dimensions = len(data_dimensions) - 1 #Last dimensions are for energy.
        energy_resolution = data_dimensions[physical_dimensions];
        limits = dataset.attrs['UpLowLimits']
        
        
        datasetDeltaReal = file['deltaReal0']
        datasetDeltaImag = file['deltaImag0']
        delta = np.array(datasetDeltaReal)
        
        
        size_x = data_dimensions[0]
        size_y = data_dimensions[1]
        
        x = np.arange(0, data_dimensions[0], 1)
        y = np.arange(limits[1], limits[0], (limits[0] - limits[1])/energy_resolution)
        energies = y
        X, Y = np.meshgrid(x, y)
        
        fig = matplotlib.pyplot.figure()
        Z = dataset[:,y_plane,:]
        sigma_discrete_units = sigma*energy_resolution/(limits[0] - limits[1])
        for xp in range(0, size_x):
        	Z[xp,:] = scipy.ndimage.filters.gaussian_filter1d(Z[xp,:], sigma_discrete_units)
        
        #Color map figure
        ax = fig.gca()
        im = ax.pcolormesh(X.transpose(), Y.transpose(), Z, cmap=matplotlib.cm.coolwarm)
        
        plt.ylim([-0.2, 0.2])
        fig.colorbar(im)
        fig.savefig('figures/LDOS_vz_x_' + vz_format +  '.png')
        
        
        fig = matplotlib.pyplot.figure()
        Z = np.zeros_like(Z)
        for i in range(size_x):
            Z[i, :] = dataset[i,i,:]
        sigma_discrete_units = sigma*energy_resolution/(limits[0] - limits[1])
        for xp in range(0, size_x):
        	Z[xp,:] = scipy.ndimage.filters.gaussian_filter1d(Z[xp,:], sigma_discrete_units)
        
        #Color map figure
        ax = fig.gca()
        im = ax.pcolormesh(X.transpose(), Y.transpose(), Z, cmap=matplotlib.cm.coolwarm)
        
        plt.ylim([-0.2, 0.2])
        fig.colorbar(im)
        fig.savefig('figures/LDOS_vz_xy_' + vz_format +  '.png')
        
        
        # Find peaks in the LDOS (Eg):
        
        
        Eg = np.zeros_like(delta)
        ratio = np.zeros_like(delta)
        
        if vz != 0:
            
            signal = dataset[y_plane, y_plane, :]
            ysr_peak = findpeakdiscrete(signal, 1)
#            print(ysr_peak)
        else:
            ysr_peak = 0
        
        for i in range(data_dimensions[0]):
            for j in range(data_dimensions[1]):
                signal = dataset[i, j, : ]
                peak_found = False
#                print(signal)
                for peak_counter in range(1,10):
                    peaks2 = findpeakdiscrete(signal, peak_counter)
                    if peaks2 != ysr_peak:
                        peak_found = True
                        if(peak_counter > 2):
                            print('oi')
                        break
                if not peak_found:
                    raise Exception("peak not found")
#                signal = scipy.ndimage.filters.gaussian_filter1d(signal, sigma_discrete_units)
#                peaks, _ = find_peaks(signal) #, height=peak_height)
#                print(energies[peaks[-2]], energies[peaks2])
                Eg[i,j] = energies[peaks2]#energies[peaks[-2]]
                ratio[i,j] = Eg[i,j]/delta[i,j]
        
        sigma = 0.002
        sigma_discrete_units = sigma*energy_resolution/(limits[0] - limits[1])
        
        
        
        
        Z1 = dataset[y_plane, y_plane, :]
        signal = Z1[: int(data_dimensions[2]/2)]
        
        
        Z2 = dataset[0, 0, :]
        plt.figure()
#        Z1 = scipy.ndimage.filters.gaussian_filter1d(Z1, sigma_discrete_units)
#        peaks, _ = find_peaks(signal, height=peak_height)
#        Z2 = scipy.ndimage.filters.gaussian_filter1d(Z2, sigma_discrete_units)
        plt.plot(y, Z2, '--', label="Bulk")
        plt.plot(y, Z1, label=r'$x_0$')
#        plt.plot(y[peaks[-1]], Z1[peaks[-1]], 'x')
#        
        lim = np.max(Z1)*1.2
#        plt.vlines([Eg[y_plane], -Eg[y_plane]], 0,  lim)
        plt.xlim([-0.2, 0.2])
        plt.legend()
        plt.xlabel(r'$E$')
        plt.ylabel(r'LDOS')
        plt.savefig('figures/LDOS_middle_vz_' + vz_format +  '.png')
        plt.ylim([-0.01*lim, 0.05*lim ])
        plt.savefig('figures/LDOS_middle_zoom_vz_' + vz_format +  '.png')
        plt.close()
        
        fig, ax = plt.subplots()
        x = np.arange(0, data_dimensions[0], 1)
        y = np.arange(0, data_dimensions[1], 1)
        X, Y = np.meshgrid(x, y)
        Z = delta
        cs = ax.contourf(X.transpose(), Y.transpose(), Z) #, cmap=matplotlib.cm.coolwarm)
        ax.contour(cs, colors='k', linewidths=0.2)
        
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        
        cbar = fig.colorbar(cs)
        cbar.ax.set_ylabel(r'$\Delta_{i}$')
        
        ax.text(text_pos[0], text_pos[1], r'$V_Z =$ ' + str(vz))
        
        fig.savefig("figures/delta_vz_" + vz_format +  ".png")
        plt.close()
        
        plt.figure()
        plt.plot(x, delta[:,y_plane])
        plt.savefig("figures/delta_profile_x_vz_" + vz_format +  ".png")
        plt.close()
        
        plt.figure()
        delta_xy = []
        for i in range(len(x)):
            delta_xy.append(delta[i,i])
        plt.plot(x, delta_xy)
        plt.savefig("figures/delta_profile_xy_vz_" + vz_format +  ".png")
        plt.close()
        

        delta_at_imp.append(delta[y_plane, y_plane])
        

        
        fig, ax = plt.subplots()
        x = np.arange(0, data_dimensions[0], 1)
        y = np.arange(0, data_dimensions[1], 1)
        X, Y = np.meshgrid(x, y)
        Z = Eg
        cs = ax.contourf(X.transpose(), Y.transpose(), Z) #, cmap=matplotlib.cm.coolwarm)
        ax.contour(cs, colors='k', linewidths=0.2)
        
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        
        cbar = fig.colorbar(cs)
        cbar.ax.set_ylabel(r'$E_{g}$')
        
        ax.text(text_pos[0], text_pos[1], r'$V_Z =$ ' + str(vz))
        
        fig.savefig("figures/Eg_vz_" + vz_format +  ".png")
        plt.close()
        
        
        
        plt.figure()
        Eg_xy = []
        for i in range(len(x)):
            Eg_xy.append(Eg[i,i])
        plt.plot(x, Eg_xy)
        plt.savefig("figures/Eg_xy_vz_" + vz_format +  ".png")
        plt.close()
        
        plt.figure()
        plt.plot(x, Eg[:,y_plane])
        plt.savefig("figures/Eg_x_vz_" + vz_format +  ".png")
        plt.close()
        
        
        fig, ax = plt.subplots()
        x = np.arange(0, data_dimensions[0], 1)
        y = np.arange(0, data_dimensions[1], 1)
        X, Y = np.meshgrid(x, y)
        Z = ratio
        cs = ax.contourf(X.transpose(), Y.transpose(), Z) #, cmap=matplotlib.cm.coolwarm)
        ax.contour(cs, colors='k', linewidths=0.2)
        
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        
        cbar = fig.colorbar(cs)
        cbar.ax.set_ylabel(r'$\frac{E_{g}}{\Delta}$')
        ax.text(text_pos[0], text_pos[1], r'$V_Z =$ ' + str(vz))
        
        fig.savefig("figures/ration_vz_" + vz_format +  ".png")
        plt.close()
                
        plt.figure()
        ratio_xy = []
        for i in range(len(x)):
            ratio_xy.append(Eg[i,i])
        plt.plot(x, ratio_xy)
        plt.savefig("figures/ratio_xy_vz_" + vz_format +  ".png")
        plt.close()
        
        plt.figure()
        plt.plot(x, ratio[:,y_plane])
        plt.savefig("figures/ratio_x_vz_" + vz_format +  ".png")
        plt.close()
        
        
        
plt.figure()
plt.plot(vz_points, delta_at_imp)
plt.savefig("figures/deltaAtImp.png")
plt.close()

        
         
