"""
Cube class, for manipulating JWST spectro-imaging data.

This class is used to create and manage Cube objects. 
JWST spectro-imaging data consists of a series of images, 
each associated with a specific wavelength. 
The array of values is then in the form of a 3D list. 
The first dimension is a spectral dimension, defaulting to 
micron wavelengths. 
The wavelength range is specific to the instrument's grism 
and/or filter. The other two dimensions are spatial dimensions, 
forming images. 

The Cube object has two headers containing all data information. 
The structure and information of the headers are identical to 
those of the files output by the reduction pipeline. 
The first header is called 'primary' and contains all general 
information about the observations (PI, instrument, 
date, time and duration of observations, configuration, etc.).
The second header provides more information about 
the data, such as 3D array size, 3-axis sampling and units. 
A summary of the information can be displayed using the .info() method.

The values (in surface brightness if units are the default) of 
the 3D array are stored in the .data attribute. 
The uncertainties at each pixel of the cube are also stored in 
an .errs attribute, an array of the same size as the data. 

When creating a Cube object, you must provide the file name in .fits format. 

Parameters
----------
file_name : str
    The name of the file in .fits format. For JWST spectro-imaging data, 
    the default name contains the suffix “_s3d”.

Attributes
----------
    primary_header : 'astropy.io.fits.Header'
        The FITS primary header, using astropy.io tools.
    data_header : 'astropy.io.fits.Header'
        The FITS header associated with the data, using astropy.io tools.
    data : array_like 
        Data stored as a cube (3D array). The first dimension is the 
        spectral dimension, the other two dimensions are the spatial dimensions. 
    errs : array_like
        Uncertainties associated with 'science' data stored in the .data attribute.
    size : array_like
        The number of points in each dimension. The first value gives 
        the number of spectral pixels, the second the number of x-axis pixels 
        and the third the number of y-axis pixels. 
    px_area : float
        Area of a spatial pixel in images. The value is given in steradian.
    units : str
        The unit of values stored in the .data table. Default values are 
        surface brightness in MJy/sr. 

"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.io import fits
from astropy.wcs import WCS

try:
    from photutils.aperture.circle import CircularAperture
    from photutils.aperture import aperture_photometry, ApertureStats
except ImportError:
    from photutils import CircularAperture, aperture_photometry, ApertureStats
    
from scipy import ndimage
from tqdm import tqdm
import warnings

from JWSToolKit.Spec import Spec

c_sp = 299792458        # Speed of light (m/s)


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
            #dq_hdu = hdul[3]
            #varPoisson_hdu = hdul[4]
            #varRnoise_hdu = hdul[5]
            #asdf_hdu = hdul[6]

            self.primary_header = primary_hdu.header                # Primary header
            self.data_header = sci_hdu.header                       # Data header
            self.data = sci_hdu.data                                # Data cube
            self.errs = err_hdu.data                                # Data cube errors
            self.size = self.data.shape                             # Data cube size
            self.px_area = float(self.data_header['PIXAR_SR'])      # Pixel area (steradian)
            self.units = self.data_header['BUNIT']



    def get_wvs(self, units='um'):
        """Returns the wavelength grid of the data cube

        Parameters
        ----------
        units : str, optional
            The character string specifying the units of the wavelengths.

        Returns
        ----------
        array_like
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

    def info(self):
        """Prints information stored in headers associated with the data cube. 
        """

        dither_bool = False
        if self.primary_header['NUMDTHPT'] > 1:
            dither_bool = True


        print()
        print('__________ DATA CUBE INFORMATION __________')
        print('Data file name:' + self.file_name)
        print('Program PI: ' + self.primary_header['PI_NAME'] + ', for the project: ' + self.primary_header['TITLE'])
        print('Program ID: ' + self.primary_header['PROGRAM'])
        print('Target: ' + self.primary_header['TARGNAME'])
        print('Telescope: ' + self.primary_header['TELESCOP'] + ' \ Instrument: ' + self.primary_header['INSTRUME'])
        print('Configuration: ' + self.primary_header['GRATING'] + ' + ' + self.primary_header['FILTER'])
        print('Number of integrations, groups and frames: ' + str(self.primary_header['NINTS']) + ', ' + str(self.primary_header['NGROUPS']) + ', ' + str(self.primary_header['NFRAMES']))
        print('Dither strategy: ' + str(dither_bool))

        if dither_bool:
            print('Dither patern type: ' + self.primary_header['PATTTYPE'])

        print()
        print('Date and time of observations: ' + self.primary_header['DATE-OBS'] + ' | ' + self.primary_header['TIME-OBS'])
        print('Target position in the sky: RA(J2000) = ' + str(self.primary_header['TARG_RA']) + ' , Dec(J2000) = ' + str(self.primary_header['TARG_DEC']))
        print('Effecive Exposure Time: ' + str(self.primary_header['EFFEXPTM']) + ' s')
        print('Total Exposure Time (with overheads): ' + str(self.primary_header['DURATION']) + ' s')

        print()

        dim_data = self.data_header['NAXIS']
        data_type = 'None'
        data_shape = []

        for i in range(dim_data):
            data_shape.append(self.data_header['NAXIS{}'.format(int(i+1))])

        if dim_data == 3:
            data_type = 'Data Cube'
            print('Data type and shape: ' + data_type + ' | ' + str(data_shape[0]) + ', ' + str(data_shape[1]) + ', ' + str(data_shape[2]) + ' (x, y, wvs)')

            pixel_unit = self.data_header['CUNIT1']

            if pixel_unit == 'deg':
                x_px_size_deg = self.data_header['CDELT1']
                y_px_size_deg = self.data_header['CDELT2']

                print('Spatial pixel sizes in ' + pixel_unit + ' (dx, dy): ' + str(x_px_size_deg) + ', ' + str(y_px_size_deg))
                print('Spatial pixel sizes in arcsec (dx, dy): ' + str(round(x_px_size_deg * 3600, 4)) + ', ' + str(round(y_px_size_deg * 3600, 4)))


            dwvs_unit = self.data_header['CUNIT3']

            if dwvs_unit == 'um':
                print('Spectral pixel size (µm): ' + str(round(self.data_header['CDELT3'], 6)))

            print('Spectral range of observations (µm): ' + str(self.data_header['WAVSTART'] * 1e6) + ' - ' + str(self.data_header['WAVEND'] * 1e6))

            print('Units of spectral pixel values: ' + self.data_header['BUNIT'])

        print()

    def extract_spec_circ_aperture(self, radius, position, err=False, units='Jy'):
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
        array_like
            If err is False, the routine returns flux values of the summed spectrum
            If err is True, the routine returns 2 sub-lists. The first containing flux values of the summed spectrum
            and the second containing erros associated with flux values.
        """

        all_units = all_units = ['Jy', 'erg s-1 cm-2 um-1', 'erg s-1 cm-2 Hz-1']

        if not isinstance(radius, int) and radius > 0:
            raise TypeError("The input radius is invalid. It must be an positive integer.")
        elif (np.size(position) != 2):
            raise Exception("The input position is invalid. The parameter must be a list of two elements and contain the spatial position of the aperture.")
        elif (not isinstance(position[0], int)) or (not isinstance(position[1], int)) or position[0] < 0 or position[1] < 0:
            raise TypeError("The input position is invalid. Values must be positive integers.")
        elif (units not in all_units) or (not isinstance(units, str)):
            raise Exception("The input units are invalid. They must be one of the following units: Jy, erg s-1 cm-2 um-1, erg s-1 cm-2 Hz-1")
        else:

            spec_len = self.size[0]
            spec = np.zeros(spec_len)

            print('__________ Spectrum extraction __________')
            for i in tqdm(range(spec_len)):

                ch_map = self.data[i,:,:]           # MJy/sr
                nan_idxs = np.isnan(ch_map)
                ch_map[nan_idxs] = 0
                ch_map *= 1e6 * self.px_area        # Jy

                apert = CircularAperture((position[0], position[1]), r=radius)
                apertstats = ApertureStats(ch_map, apert)

                spec[i] = apertstats.sum

            print()

            if units == all_units[0]:
                return spec
            else:
                spec *= 1e-23
                if units == all_units[2] :
                    return spec
                elif units == all_units[1]:
                    spec *= (c_sp * 1e6 / (self.get_wvs(units='um') ** 2))
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
        array_like
             If the coordinates of a single point have been given, the list contains two elements being the R.A., Dec.
             coordinates converted into degrees. If the coordinates are those of several points, the list contains two
             sub-lists containing respectively the R.A., Dec. positions of the different points.
        """

        warnings.filterwarnings("ignore")

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
        array_like
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

    def line_emission_map(self, wv_line, continuum_range=2000, line_width=400, continuum_degree=1, map_units='MJy um sr-1', control_plot=False):
        """Builds the integrated emission map of a line at a given wavelength

        Parameters
        ----------
        wv_line : float
            Wavelength in vacuum and at rest of the emission line, given in the same unit as the x-axis of the spectra.
        continuum_range : float, optional
            Spectral half-interval used to adjust the spectrum continuum, given in km/s. The interval is centered on the wavelength of the line.
        line_width : float, optional
            Spectral width of the emission line, given in km/s.
        continuum_degree : int, optional
            Polynomial order used to fit the continuum around the line.
        map_units : str, optional
            Map pixel units.
        control_plot : bool, optional
            If True, show the integrated emission map.

        Returns
        ----------
        array_like
            The integrated emission map, with the same dimensions as the spatial dimensions of the initial data cube. 
        """

        all_map_units = ['MJy um sr-1', 'erg s-1 cm-2 sr-1']
        wvs_units = 'um'

        if not isinstance(map_units, str) or map_units not in all_map_units:
            raise TypeError('The input units are invalid. They must be one of the following units: MJy um sr-1, erg s-1 cm-2 sr-1')
        else:

            integrated_map = np.full((self.size[1], self.size[2]), np.nan)
            err_integrated_map = np.copy(integrated_map)

            wvs_cube = self.get_wvs(units=wvs_units)


            print('__________ Construction of the integrated emission map of the line at {} µm __________'.format(wv_line))
            for i in tqdm(range(self.size[1])):
                for j in range(self.size[2]):

                    full_spectrum_array = self.data[:,i,j]
                    full_spectrum = Spec(wvs_cube, full_spectrum_array, units=self.units)

                    if map_units == 'erg s-1 cm-2 sr-1':                        # MJy/sr -> erg s-1 cm-2 sr-1 µm-1
                        full_spectrum.convert(units='erg s-1 cm-2 um-1 sr-1')


                    sliced_spectrum = full_spectrum.cut(-continuum_range, continuum_range, units='vel', wv_ref=wv_line)


                    if not True in np.isnan(sliced_spectrum.values):

                        sliced_spectrum_baseline_subtracted = sliced_spectrum.sub_baseline(wv_line=wv_line, mask_rv=line_width/2, deg=continuum_degree, control_plot=False)

                        flux_line, err_flux_line = sliced_spectrum_baseline_subtracted.line_integ(wv_line=wv_line, profile=None, line_width=line_width, control_plot=False)
                        integrated_map[i,j] = flux_line
                        err_integrated_map[i,j] = err_flux_line


            if control_plot:

                std_map = np.nanstd(integrated_map)
                mediane_map = np.nanmedian(integrated_map)

                fig, ax = plt.subplots(figsize=(6,5))

                im = ax.imshow(integrated_map, origin='lower', cmap='magma', vmin=std_map)
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                cbar = fig.colorbar(im, cax=cax)
                cbar.set_label("Integrated emission (" + map_units  +')' , fontsize=12)

                ax.set_title('Integrated line emission map @ {:.3f} µm'.format(wv_line))

                fig.tight_layout()
                plt.show()
                plt.close()



            print()

            return integrated_map




"""
- get_band_image (à renommer)
    Crée une image à partir d'un nom de filtre d'un instrument d'imagerie (Type MIRI ou NIRCam)
- getters pour les différents sub-list du fichier (dq, ...)
- convolve 
    Convolution Gaussienne en spécifiant une FWHM à une longueur d'onde donnée
- rotate 
    Tourne tout le cube de données à partir d'un angle 
"""




"""
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
"""










# ________ FIN _________
