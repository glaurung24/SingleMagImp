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

#if(len(sys.argv) != 2):
#	print( "Error, the following parameters are needed: .hdf5-file")
#	exit(1)
    
results_dir = ""

filenamebase = "vz_"
extension = "_diag_size21.hdf5"
sigma = 0.002
peak_height = 7.5

for vz in np.arange(0.1, 7.6, 0.1):
    vz_format = "{:.6f}".format(vz)
    fName = filenamebase + vz_format + extension
    if os.path.exists(results_dir + fName):

        filename = results_dir + fName
        
        file = h5py.File(filename, 'r');
        dataset = file['LDOS']
        
        
        
        data_dimensions = dataset.shape
        y_plane = int(np.floor(data_dimensions[1]/2))
        physical_dimensions = len(data_dimensions) - 1 #Last dimensions are for energy.
        energy_resolution = data_dimensions[physical_dimensions];
        limits = dataset.attrs['UpLowLimits']
        
        
        datasetDeltaReal = file['deltaReal0']
        datasetDeltaImag = file['deltaImag0']
        delta = abs(np.array(datasetDeltaReal) + 1j*np.array(datasetDeltaImag))
        
        
        size_x = data_dimensions[0]
        size_y = data_dimensions[1]
        
        x = np.arange(0, data_dimensions[0], 1)
        y = np.arange(limits[1], limits[0], (limits[0] - limits[1])/energy_resolution)
        X, Y = np.meshgrid(x, y)
        
        fig = matplotlib.pyplot.figure()
        Z = dataset[:,y_plane,:]
        sigma_discrete_units = sigma*energy_resolution/(limits[0] - limits[1])
        for xp in range(0, size_x):
        	Z[xp,:] = scipy.ndimage.filters.gaussian_filter1d(Z[xp,:], sigma_discrete_units)
        
        #Color map figure
        ax = fig.gca()
        im = ax.pcolormesh(X.transpose(), Y.transpose(), Z, cmap=matplotlib.cm.coolwarm)
        
        plt.ylim([-1, 1])
        fig.colorbar(im)
        fig.savefig('figures/LDOS_vz_' + vz_format +  '.png')
        
        sigma = 0.001
        sigma_discrete_units = sigma*energy_resolution/(limits[0] - limits[1])
        
        
        
        
        Z1 = dataset[y_plane, y_plane, :]
        signal = Z1[: int(data_dimensions[2]/2)]
        
        
        Z2 = dataset[0, 0, :]
        plt.figure()
        Z1 = scipy.ndimage.filters.gaussian_filter1d(Z1, sigma_discrete_units)
        peaks, _ = find_peaks(signal, height=peak_height)
        Z2 = scipy.ndimage.filters.gaussian_filter1d(Z2, sigma_discrete_units)
        plt.plot(y, Z1)
        plt.plot(y[peaks[-1]], Z1[peaks[-1]], 'x')
        plt.plot(y, Z2, '--')
        plt.xlim([-1, 1])
        plt.savefig('figures/LDOS_middle_vz_' + vz_format +  '.png')
        plt.close()
        
        plt.figure()
        x = np.arange(0, data_dimensions[0], 1)
        y = np.arange(0, data_dimensions[1], 1)
        X, Y = np.meshgrid(x, y)
        Z = delta
        plt.pcolormesh(X.transpose(), Y.transpose(), Z, cmap=matplotlib.cm.coolwarm)
        plt.colorbar()
        plt.savefig("figures/delta_vz_" + vz_format +  ".png")
        plt.close()
        
        plt.figure()
        plt.plot(x, delta[:,y_plane])
        plt.savefig("figures/delta_profile_vz_" + vz_format +  ".png")
        plt.close()
        
        
        # Find peaks in the LDOS (Eg):
        
        
        Eg = np.zeros_like(delta)
        ratio = np.zeros_like(delta)
        for i in range(data_dimensions[0]):
            for j in range(data_dimensions[1]):
                signal = dataset[i, j, : int(np.floor(data_dimensions[2]/2))]
                signal = scipy.ndimage.filters.gaussian_filter1d(signal, sigma_discrete_units)
                peaks, _ = find_peaks(signal, height=peak_height)
                Eg[i,j] = signal[peaks[-1]]
                ratio[i,j] = Eg[i,j]/delta[i,j]
        
        
        plt.figure()
        x = np.arange(0, data_dimensions[0], 1)
        y = np.arange(0, data_dimensions[1], 1)
        X, Y = np.meshgrid(x, y)
        Z = Eg
        plt.pcolormesh(X.transpose(), Y.transpose(), Z, cmap=matplotlib.cm.coolwarm)
        plt.colorbar()
        plt.savefig("figures/Eg_vz_" + vz_format +  ".png")
        plt.close()
        
        
        plt.figure()
        x = np.arange(0, data_dimensions[0], 1)
        y = np.arange(0, data_dimensions[1], 1)
        X, Y = np.meshgrid(x, y)
        Z = ratio
        plt.pcolormesh(X.transpose(), Y.transpose(), Z, cmap=matplotlib.cm.coolwarm)
        plt.colorbar()
        plt.savefig("figures/ration_vz_" + vz_format +  ".png")
        plt.close()
        
        
        
        
        
        
         
