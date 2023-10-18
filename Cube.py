import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from photutils import CircularAperture, aperture_photometry, ApertureStats
from scipy import ndimage

from Spec import Spec

c_sp = 299792458        # Speed of light (m/s)


"""
.info (Taille, Nom, Surf px, coords centre, instrum, wvs min/max/step, ID, filtre, units)
.copy(Cube)
.reduce(taille spatiale ou spectrale)
.plot ?
.mean_spectrum
"""

class Cube:
    def __init__(self, file_name):

        if not isinstance(file_name, str):
            raise TypeError("The input file name is invalid. It must be a character string")
        else:
            self.file_name = file_name                              # Data cube name

            hdul = fits.open(self.file_name)
            primary_hdu = hdul[0]
            sci_hdu = hdul[1]
            err_hdu = hdul[2]
            dq_hdu = hdul[3]
            varPoisson_hdu = hdul[4]
            varRnoise_hdu = hdul[5]
            asdf_hdu = hdul[6]

            self.primary_header = primary_hdu.header                # Primary header
            self.data_header = sci_hdu.header                       # Data header
            self.data = sci_hdu.data                                # Data cube
            self.errs = err_hdu.data                                # Data cube errors
            self.size = self.data.shape                             # Data cube size
            self.px_area = float(self.data_header['PIXAR_SR'])      # Pixel area (steradian)
            self.instrument = self.primary_header['INSTRUME']       # Instrument
    def wvs(self, units='um'):
        """Returns the wavelength grid of the data cube

        Parameters
        ----------
        units : str, optional
            The character string specifies the units of the wavelengths.

        Returns
        ----------
        list
            The wavelength grid
        """

        all_units = ['um', 'A', 'nm']

        if (not isinstance(units, str)) or (units not in all_units):
            raise Exception("The input units are invalid. They must be one of the following units: um, A or nm.")
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
    def extract_spec(self, radius, position, err=False, units='Jy'):
        """Extracts a summed spectrum in a circular aperture

        Parameters
        ----------
        radius : int
            Radius in pixel of the integration aperture
        position : list
            Position in pixels of the aperture center. It must contains two values: the horizontale and
            verticale coordinates respectively.
        err : bool, optional
            If True, return the errors of each spectral flux value.
        units : str, optional
            The character string specifies the units of the output spectrum.
        Returns
        ----------
        list
            If err is False, the routine returns flux values of the summed spectrum
            If err is True, the routine returns 2 sub-lists. The first containing flux values of the summed spectrum
            and the second containing erros associated with flux values.
        """

        all_units = all_units = ['Jy', 'erg s-1 cm-2 micron-1', 'erg s-1 cm-2 Hz-1']

        if not isinstance(radius, int) and radius > 0:
            raise TypeError("The input radius is invalid. It must be an positive integer.")
        elif (np.size(position) != 2):
            raise Exception("The input position is invalid. The parameter must be a list of two elements and contain the spatial position of the aperture.")
        elif (not isinstance(position[0], int)) or (not isinstance(position[1], int)) or position[0] < 0 or position[1] < 0:
            raise TypeError("The input position is invalid. Values must be positive integers.")
        elif (units not in all_units) or (not isinstance(units, str)):
            raise Exception("The input units are invalid. They must be one of the following units: Jy, erg s-1 cm-2 micron-1, erg s-1 cm-2 Hz-1")
        else:

            spec_len = self.size[0]
            spec = np.zeros(spec_len)

            for i in range(spec_len):

                ch_map = self.data[i,:,:]           # MJy/sr
                nan_idxs = np.isnan(ch_map)
                ch_map[nan_idxs] = 0
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
    def get_world_coords(self, coords):
        """Returns the coordinates in degrees (R.A., Dec.) of one or more pixel positions in the data cube.

        Parameters
        ----------
        coords : list
            Coordinates in pixels to be converted into degrees. It can contain two elements (corresponding to the
            position of a single point) or two sub-lists containing the horizontal and vertical positions of several
            points respectively.

        Returns
        ----------
        list
             If the coordinates of a single point have been given, the list contains two elements being the R.A., Dec.
             coordinates converted into degrees. If the coordinates are those of several points, the list contains two
             sub-lists containing respectively the R.A., Dec. positions of the different points.
        """

        if not isinstance(coords, list):
            raise TypeError('The input coordinates are invalid. They must be a list of two elements or a list of sublist as follow: [[x1, x2, ...], [y1, y2, ...]]')
        else:

            sci_header_mod = self.data_header.copy()
            sci_header_mod["NAXIS"] = 2
            sci_header_mod["WCSAXES"] = 2
            for keyword in ["CTYPE3", "CRVAL3", "CDELT3", "CRPIX3", "CUNIT3", "PC1_3", "PC2_3", "PC3_1", "PC3_2", "PC3_3"]:
                del sci_header_mod[keyword]

            wcs_sci = WCS(sci_header_mod)
            coords_proj = wcs_sci.pixel_to_world_values(coords[0], coords[1])

            return coords_proj
    def get_px_coords(self, coords):
        """Returns the coordinates in pixels (x,y) of one or more pixel positions in the data cube.

        Parameters
        ----------
        coords : list
            Coordinates in degrees (R.A., Dec.) to be converted into pixel coordinates. It can contain two elements
            (corresponding to the position of a single point) or two sub-lists containing the R.A. and Dec. positions of
            several points respectively.

        Returns
        ----------
        list
             If the coordinates of a single point have been given, the list contains two elements being the (x,y)
             coordinates converted into pixel coordinates. If the coordinates are those of several points, the list
             contains two sub-lists containing respectively the x and y positions of the different points.
        """

        if not isinstance(coords, list):
            raise TypeError('The input coordinates are invalid. They must be a list of two elements or a list of sublist as follow: [[x1, x2, ...], [y1, y2, ...]]')
        else:

            sci_header_mod = self.data_header.copy()
            sci_header_mod["NAXIS"] = 2
            sci_header_mod["WCSAXES"] = 2
            for keyword in ["CTYPE3", "CRVAL3", "CDELT3", "CRPIX3", "CUNIT3", "PC1_3", "PC2_3", "PC3_1", "PC3_2", "PC3_3"]:
                del sci_header_mod[keyword]

            wcs_sci = WCS(sci_header_mod)
            coords_proj = wcs_sci.world_to_pixel_values(coords[0], coords[1])

            return coords_proj
    def rotate(self, angle):

        new_cube = np.full(self.size, np.nan)
        cube_flag = np.full(self.size, 1)

        for k in range(self.size[0]):
            ch_map = self.data[k,:,:]
            max_ch_map = np.nanmax(ch_map)
            ch_map += max_ch_map

            nan_idxs = np.isnan(ch_map)
            ch_map_masked = np.copy(ch_map)
            ch_map_masked[nan_idxs] = -1e6

            ch_map_rotated = ndimage.rotate(ch_map_masked, angle, reshape=False)

            ch_map_rotated[np.where(ch_map_rotated <= 0)] = np.nan
            ch_map_rotated -= max_ch_map

            new_cube[k,:,:] = ch_map_rotated

        return new_cube

