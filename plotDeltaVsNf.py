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
extension = "_diag_size21_normalState.hdf5" #"_chebychev_GPUsize20.hdf5"
figure_extension =  "_diagDelta0_1_normalState"

sigma = 0.2
peak_height = 7.5

delta_bulk = 0.12

delta_at_imp = []

vz_points = []

text_pos = [1,1]

interesting_vz_values = [0, 1.0, 1.75]

delta_fig_x, delta_ax_x = plt.subplots()
delta_fig_xy, delta_ax_xy = plt.subplots()
ldos_fig_xy, ldos_ax_xy = plt.subplots()

for vz in np.arange(0.0, 10.1, 0.0025):
    vz_format = "{:.6f}".format(vz)
    fName = filenamebase + vz_format + extension
    if os.path.exists(results_dir + fName) and vz in interesting_vz_values:
        
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
        
        sigma_discrete_units = sigma*energy_resolution/(limits[0] - limits[1])
        
        x = np.arange(0, data_dimensions[0], 1)
        y = np.arange(limits[1], limits[0], (limits[0] - limits[1])/energy_resolution)
        energies = y
        X, Y = np.meshgrid(x, y)
        

        delta_ax_x.plot(x, delta[:,y_plane], label=r'$V_Z = $' + str(vz) )

        
        delta_xy = []
        for i in range(len(x)):
            delta_xy.append(delta[i,i])
        delta_ax_xy.plot(x, delta_xy, label=r'$V_Z = $' + str(vz))
        
        
        Z1 = dataset[y_plane, y_plane, :]
        signal = Z1[: int(data_dimensions[2]/2)]


        


        # 300 represents number of points to make between T.min and T.max
        xnew = np.linspace(y.min(), y.max(), 300000)  
        
        spline = interpolate.splrep(y, Z1, k=1)
        ldos_new = interpolate.splev(xnew, spline)
        ldos_new = scipy.ndimage.filters.gaussian_filter1d(ldos_new, sigma_discrete_units)
        ldos_ax_xy.plot(xnew/delta_bulk, ldos_new, label=r'$V_Z = $' + str(vz))
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

ldos_ax_xy.set_xlabel(r'$E\,/\,\Delta$')
ldos_ax_xy.set_ylabel(r'LDOS')
ldos_ax_xy.legend(loc='upper left')
ldos_ax_xy.set_xlim([-2, 2])
ldos_fig_xy.savefig("figures/LDOS_compare.png")

