"""
    pepico_library

    author:: Suddhasattwa Mandal
    contact:: m.suddhasattwa@gmail.com
    copyright:: Suddhasattwa Mandal 2021
"""
import os
import math
import h5py
import subprocess
import numpy as np
from scipy import ndimage
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from tkinter import messagebox
from tkinter import simpledialog


class pepico:

    def __init__(self, file_name, file_path, tof_address, vmi_address, channel_address):

        self.file_name = file_name
        self.file_path = file_path
        self.full_file_path = os.path.join(self.file_path, self.file_name)

        try:
            self.data_file = h5py.File(self.full_file_path, 'r')
        except Exception as e:
            print(e)
        else:
            self.channel_address = channel_address
            self.channel = self.data_file[self.channel_address]
            self.channel = np.array(self.channel[0])

            self.tof_address = tof_address
            self.tof = self.data_file[self.tof_address]
            self.tof = np.array(self.tof[0])

            self.vmi_address = vmi_address
            self.vmi = self.data_file[self.vmi_address]
            self.vmi = np.array(self.vmi)

            self.data_file.close()

    def get_tof(self, channel_requested):

        tof_out = []
        for index, channel_id in enumerate(self.channel):
            if channel_id == channel_requested:
                tof_out.append(self.tof[index])
        tof_out = np.array(tof_out)

        return tof_out

    def get_vmi(self, channel_requested):

        vmi_out = []
        for index, channel_id in enumerate(self.channel):
            if channel_id == channel_requested:
                vmi_out.append([self.vmi[0, index], self.vmi[1, index]])
        vmi_out = np.array(vmi_out)
        vmi_out = vmi_out.transpose()

        return vmi_out

    @staticmethod
    def hist_1d(data, start, end, nbin):

        bins = np.linspace(start, end, nbin)
        hist, bins = np.histogram(data, bins)
        axis = (bins[:-1] + bins[1:]) / 2

        return axis, hist

    @staticmethod
    def hist_2d(xy_data, x_start, x_end, x_nbin, y_start, y_end, y_nbin):

        x_data = xy_data[0, :]
        y_data = xy_data[1, :]
        x_bins = np.linspace(x_start, x_end, x_nbin)
        y_bins = np.linspace(y_start, y_end, y_nbin)
        xy_hist, x_bins, y_bins = np.histogram2d(x_data, y_data, bins=(x_bins, y_bins))
        x_axis = (x_bins[:-1] + x_bins[1:]) / 2
        y_axis = (y_bins[:-1] + y_bins[1:]) / 2
        xy_hist = xy_hist.transpose()

        return xy_hist, x_axis, y_axis

    @staticmethod
    def get_tof_correlated_vmi(tof_data, vmi_data, tof_lim):

        vmi_out = []
        for index, tof in enumerate(tof_data):
            if tof >= tof_lim[0] and tof <= tof_lim[1]:
                vmi_out.append([vmi_data[0, index], vmi_data[1, index]])
        vmi_out = np.array(vmi_out)
        vmi_out = vmi_out.transpose()

        return vmi_out

    @staticmethod
    def invert_image(image_name, centre, angle, radius, crop):

        imagename, extn = os.path.splitext(image_name)
        image_in = np.loadtxt(f'./ana_dat/{image_name}.dat', delimiter='\t')

        cx, cy = centre
        image_translated = np.zeros(image_in.shape)
        nrow, ncol = image_in.shape
        for irow, row in enumerate(image_in):
            for icol, num in enumerate(row):
                ir = irow + round(nrow / 2) - cy
                ic = icol + round(ncol / 2) - cx
                if (ir >= 0 and ir < nrow) and (ic >= 0 and ic < ncol):
                    image_translated[ir, ic] = num

        image_rotated = ndimage.rotate(image_translated, angle)

        if crop:
            mask = np.zeros(image_rotated.shape)
            nrow, ncol = image_rotated.shape
            crow = round(nrow / 2)
            ccol = round(ncol / 2)
            for irow, row in enumerate(mask):
                for icol, num in enumerate(row):
                    iradius = ((irow - crow) ** 2 + (icol - ccol) ** 2) ** 0.5
                    if iradius <= radius:
                        mask[irow, icol] = 1
            image_out = mask * image_rotated
        else:
            image_out = image_rotated

        np.savetxt(f'./ana_dat/Proc_{imagename}.dat', image_out, delimiter='\t')

        nrow, ncol = image_out.shape
        ix = round(ncol/2)
        iz = round(nrow/2)
        command = f'F2QC.elf ./ana_dat/Proc_{imagename}.dat -IX{ix} -IZ{iz} -M0'
        subprocess.run(command, shell=True)

        command = 'Meveler2.elf DefaultQ.dat'
        subprocess.run(command, shell=True)

        im1 = np.loadtxt('MEXmap.dat', delimiter=',')
        im2 = np.fliplr(im1)
        im3 = np.flipud(im2)
        im4 = np.fliplr(im3)
        mev1 = np.append(im4, im1, axis=0)
        mev2 = np.append(im3, im2, axis=0)
        mev_image = np.append(mev2, mev1, axis=1)

        with open('MEXdis.dat', 'r') as f:
            data = []
            for line in f:
                line = line.replace('D', 'E')
                data.append(line.split())
        mexdis = np.array(data, dtype=float)


        np.savetxt(f'./ana_dat/MEV_{image_name}.dat', mev_image, delimiter='\t')
        np.savetxt(f'./ana_dat/MEXdis_{image_name}.dat', mexdis, delimiter='\t')

        os.remove('DefaultQ.dat')
        os.remove('MEXini.dat')
        os.remove('MEXdis.dat')
        os.remove('MEXmap.dat')
        os.remove('MEXres.dat')
        os.remove('MEXsim.dat')

        return image_out, mev_image, mexdis




    @staticmethod
    def calibrate_energy(xdata, ydata, calibration_factor):

        cali_data = np.zeros([len(xdata), 2])
        for i in range(len(xdata)):
            cali_data[i, 0] = calibration_factor[0]+calibration_factor[1]*xdata[i]+calibration_factor[2]*(xdata[i]**2)+calibration_factor[3]*(xdata[i]**3)+calibration_factor[4]*(xdata[i]**4)+calibration_factor[5]*(xdata[i]**5)
            cali_data[i, 1] = ydata[i]/(calibration_factor[1]+2*calibration_factor[2]*xdata[i]+3*calibration_factor[3]*(xdata[i]**2)+4*calibration_factor[4]*(xdata[i]**3)+5*calibration_factor[5]*(xdata[i]**4))

        return cali_data

    @staticmethod
    def calibrate_tof(xdata, ydata):

        def tellme(s):
            plt.title(s, fontsize=16)
            plt.draw()

        def Gaussian(x, A, sigma, mu):
            y = A*(1/(math.sqrt(2*math.pi)*sigma))*np.exp(-0.5*((x-mu)/sigma)**2)
            return y

        def line(x, a, b):
            y = a*x + b
            return y

        fig1, ax1 = plt.subplots()
        ax1.plot(xdata, ydata)
        npeak = 0
        ToF_peak = []
        while True:
            tellme(f'Peak no.{str(npeak+1)}: Press any key to continue.')
            keypressed = plt.waitforbuttonpress()
            if keypressed:
                tellme('Select both side of the peak')
                pts = plt.ginput(2)
                (x0, y0), (x1, y1) = pts
                xmin, xmax = sorted([x0, x1])
                # gaussian fitting
                xdata_fit = xdata[(xdata >= xmin) & (xdata <= xmax)]
                ydata_fit = ydata[(xdata >= xmin) & (xdata <= xmax)]
                mm = np.amax(ydata_fit)
                ii = np.argmax(ydata_fit)
                ff = abs(ydata_fit - mm/2)
                m1 = np.amin(ff)
                lim1 = np.argmin(ff)
                ff[lim1] = ff[lim1]+np.amax(ff)
                m2 = np.amin(ff)
                lim2 = np.argmin(ff)
                sigma0 = abs(xdata_fit[lim1]-xdata_fit[lim2])/(2*math.sqrt(2*math.log(2)))
                guess = [0.7*mm, sigma0, xdata_fit[ii]]
                parameters, covariance = curve_fit(Gaussian, xdata_fit, ydata_fit, p0=guess)
                A = parameters[0]
                sigma = parameters[1]
                mu = parameters[2]
                fit_y = Gaussian(xdata_fit, A, sigma, mu)
                ax1.plot(xdata, ydata, label='data')
                ax1.plot(xdata_fit, fit_y, label='fit')
                npeak = npeak + 1
                ToF_peak.append(mu)
                answer = messagebox.askyesno("Question", "Do you want to add mass peak?")
                if not answer:
                    tellme('All Done!')
                    break

        mass_sqrt_peak = []
        for i in range(npeak):
            m = simpledialog.askfloat("Question", f"Eneter the correseponding masss for peak no.{str(i+1)}")
            mass_sqrt_peak.append(math.sqrt(m))


        # fitting
        mass_sqrt = np.array(mass_sqrt_peak)
        ToF = np.array(ToF_peak)
        a0 = (ToF[-1]-ToF[0])/(mass_sqrt[-1]-mass_sqrt[0])
        b0 = ToF[0]-a0*mass_sqrt[0]
        guess = [a0, b0]
        print(guess)
        parameters, covariance = curve_fit(line, mass_sqrt, ToF, p0=guess)
        a = parameters[0]
        b = parameters[1]
        print(parameters)
        fit_y = line(mass_sqrt, a, b)
        fig2, ax2 = plt.subplots()
        ax2.plot(mass_sqrt, ToF, label='data')
        ax2.plot(mass_sqrt, fit_y, label='fit')

        mass_data = ((xdata-b)/a)**2
        fig3, ax3 = plt.subplots()
        ax3.plot(mass_data, ydata)
        data_to_save = np.array([mass_data, ydata])
        data_to_save = data_to_save.transpose()
        np.savetxt(f'mass_cali_coeff.dat', parameters, delimiter=' ')


        plt.show()

    @staticmethod
    def calibrate_mass(tof, intensity, mass_calibration):
        mass_spec = np.zeros([2, tof.size])
        for i, val in enumerate(tof):
            mass_spec[0, i] = ((val-mass_calibration[1])/mass_calibration[0])**2
            mass_spec[1, i] = intensity[i]*(0.5*mass_calibration[0]**2)/(val-mass_calibration[1])
        return mass_spec[0,:], mass_spec[1,:]

#     @staticmethod
#     def hist_1d_rebining(xy_data, x_start, x_end, nbin):

#         x_data = xy_data[0, :]
#         y_data = xy_data[1, :]
#         x_bins = np.linspace(x_start, x_end, x_nbin)

#         xy_hist, x_bins, y_bins = np.histogram2d(x_data, y_data, bins=(x_bins, y_bins))
#         x_axis = (x_bins[:-1] + x_bins[1:]) / 2
#         y_axis = (y_bins[:-1] + y_bins[1:]) / 2
#         xy_hist = xy_hist.transpose()

#         x_data = np.linspace(0, 10, 24)
#         y_data = np.random.randint(0,10,24)
#         x_start = 0
#         x_end = 10
#         x_nbin = 10
#         dx = (x_end-x_start)/x_nbin
#         x_bin_edges = np.linspace(x_start, x_end, x_nbin+1)
#         x_axis = (x_bin_edges[:-1] + x_bin_edges[1:]) / 2
#         y_axis = np.zeros_like(x_axis)
#         for i in range(x_nbin):
#             for j, val in enumerate(y_data):
#                 if((x_data[j]>=x_axis[i]-dx/2)&(x_data[j]<x_axis[i]+dx/2)):
#                     y_axis[i]=y_axis[i]+val

#         return axis, hist
