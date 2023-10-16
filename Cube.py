import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils import CircularAperture, aperture_photometry, ApertureStats

c_sp = 299792458        # Speed of light (m/s)

class Cube:
    def __init__(self, file_name):

        if not isinstance(file_name, str):
            raise TypeError("The input file name is invalid. It must be a character string")
        else:
            self.file_name = file_name

            hdul = fits.open(self.file_name)
            primary_hdu = hdul[0]
            sci_hdu = hdul[1]
            err_hdu = hdul[2]
            dq_hdu = hdul[3]
            varPoisson_hdu = hdul[4]
            varRnoise_hdu = hdul[5]
            asdf_hdu = hdul[6]

            self.primary_header = primary_hdu.header
            self.data_header = sci_hdu.header
            self.data = sci_hdu.data
            self.size = self.data.shape
            self.px_area = float(self.data_header['PIXAR_SR'])      # steradian


    def wvs(self, units='um'):

        all_units = ['um', 'A', 'nm']

        if (not isinstance(units, str)) or (units not in all_units):
            raise Exception("The input units are invalid. It must be one of the following units: um, A or nm.")
        else:
            head = self.data_header

            ref_pix = float(head['CRPIX3'])
            lambda_ref = float(head['CRVAL3'])
            inc_lambda = float(head['CDELT3'])
            nw = int(head['NAXIS3'])

            wgrid = lambda_ref + (np.arange(nw) - ref_pix + 1) * inc_lambda     # µm

            if units == all_units[0]: return wgrid
            elif units == all_units[1]: return wgrid * 1e4
            elif units == all_units[2]: return wgrid * 1e3


    def extract_spec(self, radius, position, err=None, units='Jy'):

        all_units = all_units = ['Jy', 'erg s-1 cm-2 micron-1', 'erg s-1 cm-2 Hz-1']

        if not isinstance(radius, int) and radius > 0:
            raise TypeError("The input radius is invalid. It must be an positive integer.")
        elif (np.size(position) != 2):
            raise Exception("The input position is invalid. The parameter must be a list of two elements and contain the spatial position of the aperture.")
        elif (not isinstance(position[0], int)) or (not isinstance(position[1], int)) or position[0] < 0 or position[1] < 0:
            raise TypeError("The input position is invalid. Values must be positive integers.")
        elif (units not in all_units) or (not isinstance(units, str)):
            raise Exception("The input units are not invalid. They must be one of the following units: Jy, erg s-1 cm-2 micron-1, erg s-1 cm-2 Hz-1")
        else:

            spec_len = self.size[0]
            spec = np.zeros(spec_len)

            for i in range(spec_len):

                ch_map = self.data[i,:,:]           # MJy/sr
                nan_map = np.isnan(ch_map)
                ch_map[nan_map] = 0
                ch_map *= 1e6 * self.px_area        # Jy

                apert = CircularAperture((position[0], position[1]), r=radius)
                apertstats = ApertureStats(ch_map, apert)

                spec[i] = apertstats.sum


            if units == all_units[0]:
                return spec
            else:
                spec *= 1e-23
                if units == all_units[2] :
                    return spec
                elif units == all_units[1]:
                    spec *= (c_sp * 1e6 / (self.wvs(units='um') ** 2))
                    return spec









