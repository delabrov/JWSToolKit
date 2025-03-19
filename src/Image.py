"""
Image class, for manipulating JWST image data.

This class is used to create and manipulate Image objects. 
Many of the JWST instruments provide data in image form (2D arrays). 
Reduced images in the JWST pipeline have the default suffix “_i2d”. 
This class can also be used to manipulate images generated from a 
data cube. 

The Image object has two headers, both containing information about 
the data. The structure and information of the headers are identical to 
those of the files output by the reduction pipeline. 
The first header is called 'primary' and contains all general 
information about the observations (PI, instrument, 
date, time and duration of observations, configuration, etc.).
The second header provides more information about 
the data, such as 2D array size, 2-axis sampling and units. 
A summary of the information can be displayed using the .info() method.

The values (in surface brightness if units are the default) of 
the 2D array are stored in the .data attribute. 
The uncertainties at each pixel of the image are also stored in 
an .errs attribute, an array of the same size as the data. 

When creating a Cube object, you must provide the file name in .fits format. 

Parameters
----------
file_name : str
    The name of the file in .fits format. For JWST imaging, 
    the default name contains the suffix “_i2d”.

Attributes
----------
    primary_header : 'astropy.io.fits.Header'
        The FITS primary header, using astropy.io tools.
    data_header : 'astropy.io.fits.Header'
        The FITS header associated with the data, using astropy.io tools.
    data : array_like 
        Data stored as an image (2D array). The two dimensions are the 
        spatial dimensions. 
    errs : array_like
        Uncertainties associated with 'science' data stored in the .data attribute.
    size : array_like
        The number of points in each dimension.
    px_area : float
        Area of the spatial pixels. The value is given in steradian.
    units : str
        The unit of values stored in the .data table. Default values are 
        surface brightness in MJy/sr. 
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.colors as colors
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import rotate, gaussian_filter
from scipy.special import voigt_profile
from scipy.signal import fftconvolve
from tqdm import tqdm
import warnings


class Image:

    def __init__(self, file_name):

        if not isinstance(file_name, str):
            raise TypeError("The input file name is invalid. It must be a character string")
        else:

            self.file_name = file_name                              # Image name
            self.primary_header, self.data_header, self.data, self.errs = self._load_fits(file_name)

            self.size = self.data.shape                             # Image size
            self.px_size = float(self.data_header['CDELT1']) * 3600 # Spatial pixel size (arcsec)
            self.px_area = float(self.data_header['PIXAR_SR'])      # Pixel area (steradian)
            self.units = self.data_header['BUNIT']                  # Values unit

    
    @classmethod
    def from_file_extension(cls, primary_header, data_header, data, errs=None):
        """Builds a 'Image' object from file headers and data.
        
        Parameters
        -----------
        primary_header : astropy.io.fits.header.Header
            The JWST image primary header, extract with astropy.io.
        data_header : astropy.io.fits.header.Header
            The science header for JWST images, extract with astropy.io.
        data : array_like
            Values from the image, stored in a 2D array.
        errs : array_like, optional
            Error data associated with the data array, stored in a 2D array

        Returns
        ---------
        Image object
            A Image object.
        """

        obj = cls.__new__(cls) 
        obj.file_name = None  
        obj.primary_header = primary_header
        obj.data_header = data_header
        obj.data = data
        obj.errs = errs
        obj.size = obj.data.shape                                       # Array shape
        obj.px_size = float(obj.data_header['CDELT1']) * 3600           # Spatial pixel size (arcsec)
        obj.px_area = float(obj.data_header['PIXAR_SR'])                # Pixel area (steradian)
        obj.units = obj.data_header['BUNIT']                            # Values unit

        return obj


    def _load_fits(self, file_name: str):
        """Returns file headers and data in .fits format

        Parameters
        -----------
        file_name : str
            The name of the file in .fits format.

        Returns
        ---------
        list 
            The primary header, data header, data and file errors.
        """

        hdul = fits.open(self.file_name)
        primary_hdu = hdul[0]
        sci_hdu = hdul[1]
        err_hdu = hdul[2]

        primary_header  = primary_hdu.header
        data_header     = sci_hdu.header
        data            = sci_hdu.data
        errs            = err_hdu.data

        return primary_header, data_header, data, errs

    def info(self):
        """Prints information stored in headers associated with the image. 
        """

        dither_bool = False
        if self.primary_header['NUMDTHPT'] > 1:
            dither_bool = True


        print()
        print('__________ IMAGE INFORMATION __________')
        if self.file_name != None:
            print('Data file name:' + self.file_name)
        else:
            print('No file name or unknown file.')
        print('Program PI: ' + self.primary_header['PI_NAME'] + ', for the project: ' + self.primary_header['TITLE'])
        print('Program ID: ' + self.primary_header['PROGRAM'])
        print('Target: ' + self.primary_header['TARGNAME'])
        print('Telescope: ' + self.primary_header['TELESCOP'] + ' \\ Instrument: ' + self.primary_header['INSTRUME'])
        print('Configuration:')
        print('     Detector: ' + self.primary_header['DETECTOR'])
        if self.primary_header['INSTRUME'] == 'NIRCAM':
            print('     Channel: ' + self.primary_header['CHANNEL'])            
        print('     Filter: ' + self.primary_header['FILTER'])
        if self.primary_header['INSTRUME'] == 'NIRCAM':
            print('     Pupil: ' + self.primary_header['PUPIL'])                
        
        print('Number of integrations, groups and frames: ' + str(self.primary_header['NINTS']) + ', ' + str(self.primary_header['NGROUPS']) + ', ' + str(self.primary_header['NFRAMES']))
        print('Dither strategy: ' + str(dither_bool))

        if dither_bool:
            print('Dither patern type: ' + self.primary_header['PATTTYPE'])
            if self.primary_header['INSTRUME'] == 'NIRCAM':
                print('Primary dither points: ' + str(self.primary_header['PRIDTYPE']) + ' \\ # points: ' + str(self.primary_header['PRIDTPTS']))    
            print('Total points in pattern: ' + str(self.primary_header['NUMDTHPT']))

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

        if dim_data == 2:
            data_type = 'Image'
            print('Data type and shape: ' + data_type + ' | ' + str(data_shape[0]) + ', ' + str(data_shape[1]) + ' (x, y)')

            pixel_unit = self.data_header['CUNIT1']

            if pixel_unit == 'deg':
                x_px_size_deg = self.data_header['CDELT1']
                y_px_size_deg = self.data_header['CDELT2']

                print('Spatial pixel sizes in ' + pixel_unit + ' (dx, dy): ' + str(x_px_size_deg) + ', ' + str(y_px_size_deg))
                print('Spatial pixel sizes in arcsec (dx, dy): ' + str(round(x_px_size_deg * 3600, 4)) + ', ' + str(round(y_px_size_deg * 3600, 4)))

            print('Unit of pixel values: ' + self.data_header['BUNIT'])

        print()

    def plot(self, scale: str = 'lin', 
            use_wcs: bool = False, 
            lims: list[float] = None, 
            abs_transform: bool = False, 
            save: bool = False, 
            colorbar: bool = False,
            origin_arcsec: list[float] = None, 
            draw_compass: bool = False):

        """Display the image via matplotlib 

        Parameters
        -----------
        scale : str, optional
            Transformation for normalizing image values, i.e. colorbar scaling. Accepted 
            transformations are : 'lin', 'log', 'sqrt', 'asinh'.
        use_wcs : bool, optional
            If True, the figure axes are given in RA Dec world coordinates using the wcs 
            of the observations.
        lims : list, optional
            The minimum and maximum values to be displayed on the image (may depend 
            on the type of normalization of the figure). Limits must respect the form [x,y].
        abs_transform : bool, optional 
            If True, displays the absolute value of the image.
        save : bool, optional
            If True, saves the figure in png format.
        colorbar : bool, optional 
            If True, displays the color scale on figure.
        origin_arcsec : list, optional 
            The pixel position of the axis origin converted to arcsec. It must be given 
            in the form [x0, y0].
        Returns
        -----------
        """

        warnings.filterwarnings("ignore")

        all_scales = ['lin', 'log', 'asinh', 'sqrt']
        cmap = 'inferno'
        img = self.data

        vmin = np.nanmin(img)
        vmax = np.nanmax(img)
        
        if abs_transform:
            img = abs(self.data)

        if lims != None:
            vmin = lims[0]
            vmax = lims[1]


        if scale == all_scales[0]:
            normalization = colors.Normalize(vmin=vmin, vmax=vmax)
        elif scale == all_scales[1]:
            normalization = colors.LogNorm(vmin=vmin, vmax=vmax)
        elif scale == all_scales[2]:
            normalization = colors.AsinhNorm(vmin=vmin, vmax=vmax)
        elif scale == all_scales[3]:
            normalization = colors.PowerNorm(gamma=0.5, vmin=vmin, vmax=vmax)
        else:
            print("The normalization mode given as a parameter is invalid; those allowed are: 'lin', 'log', 'asinh', 'sqrt'")


        wcs = WCS(self.data_header)

        if use_wcs:
            fig, ax = plt.subplots(subplot_kw={'projection': wcs})
        else:
            fig, ax = plt.subplots()


        if origin_arcsec != None:
            
            x0, y0 = origin_arcsec
            x_axis = (np.arange(self.size[1]) - x0) * self.px_size
            y_axis = (np.arange(self.size[0]) - y0) * self.px_size

            img_mpl = ax.pcolormesh(x_axis, y_axis, img, cmap=cmap, norm=normalization)

        else:
            img_mpl = ax.imshow(img, cmap=cmap, origin='lower', norm=normalization)


        if colorbar:
            fig.colorbar(img_mpl, pad=0.05, label='Pixel values (' + self.units + ')')


        if use_wcs:
            if origin_arcsec != None:
                raise Exception("You cannot specify an arcsec axis origin and display the RA Dec coordinates of the image. ")
            else:
                ax.grid(color='grey', ls='--')
                ax.set_ylabel('Right Ascension (RA J2000)')
                ax.set_xlabel('Declination (Dec J2000)')

        if origin_arcsec != None:
            ax.set_xlabel('$\Delta$X (arcsec)')
            ax.set_ylabel('$\Delta$Y (arcsec)')

        if draw_compass:
            
            NE_convention = - np.radians(90)   # degree
            l_arrow = 0.1               # % Axis full size
            pad = 0.06                  # % Axis full size

            xc, yc = 0.20, 0.80         # % Axis full size
            
            compass_color = 'lightgrey'
            
            # Image orientation
            wcs = WCS(self.data_header)
            wcs_matrix = wcs.wcs.pc
            angle = np.degrees(np.arctan2(wcs_matrix[1,0], wcs_matrix[0, 0]))
            angle_compass = np.arctan2(wcs_matrix[1,0], wcs_matrix[0, 0]) + NE_convention

        
            xN, yN = np.cos(angle_compass) * l_arrow, np.sin(angle_compass) * l_arrow
            xE, yE = -np.sin(angle_compass) * l_arrow, np.cos(angle_compass) * l_arrow

            x_N, y_N = np.cos(angle_compass) * (l_arrow + pad), np.sin(angle_compass) * (l_arrow + pad)
            x_E, y_E = -np.sin(angle_compass) * (l_arrow + pad), np.cos(angle_compass) * (l_arrow + pad)

            ax.arrow(x=xc, y=yc, dx=xN, dy=yN, color=compass_color, transform=ax.transAxes, head_width=0.015)
            ax.arrow(x=xc, y=yc, dx=xE, dy=yE, color=compass_color, transform=ax.transAxes, head_width=0.015)
            #ax.annotate("", xytext=(xc, yc), xy=(xN, yN), arrowprops=dict(arrowstyle="->"), color=compass_color)
            #ax.annotate("", xytext=(xc, yc), xy=(xE, yE), arrowprops=dict(arrowstyle="->"), color=compass_color)
            ax.text(xc + x_N, yc + y_N, 'N', color=compass_color, ha='center', va='center', fontsize=15, transform=ax.transAxes)
            ax.text(xc + x_E, yc + y_E, 'E', color=compass_color, ha='center', va='center', fontsize=15, transform=ax.transAxes)

        fig.tight_layout()
        if save:
            fig.savefig('image.png', dpi=300)
        #plt.show()

    def save_as_dat(self, filename: str = None):
        """Saves the image as a .dat file

        Parameters
        ----------
        filename : str, optional
            Output file name.
        
        Returns 
        ----------
        """

        data_filename_dat = 'jwst_image.dat'
        file_comments = '##################################################################################\n' + 'JWST image | Unit: ' + self.units + ' | Px Size: ' + str(self.px_size) + ' (arcsec)\n' + '##################################################################################' 

        if filename != None:
            data_filename_dat = filename + '.dat'

        np.savetxt(data_filename_dat, self.data, fmt='%f', delimiter='  ', header=file_comments)

        print()
        print("__________ Image saved successfully __________")
        print()
        
    def save_as_fits(self, filename: str = None):
        """Saves the image as a .fits file

        Parameters
        ----------
        filename : str, optional
            Output file name.
        
        Returns 
        ----------
        """

        new_primary_header = self.primary_header.copy()
        new_primary_header['COMMENT'] = 'Edited with JWSToolKit'

        primary_hdu = fits.PrimaryHDU(header=self.primary_header)
        science_hdu = fits.ImageHDU(data=self.data, header=self.data_header)

        hdul = fits.HDUList([primary_hdu, science_hdu])

        if filename != None:
            hdul.writeto(filename + '.fits', overwrite=True)
        else:
            hdul.writeto('new_jwst_image.fits', overwrite=True)

        print()
        print("__________ Image saved successfully __________")
        print()
   
    def get_px_coords(self, coords: list):
        """Returns the coordinates in pixels (x,y) of one or more pixel positions in the image.

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
            wcs_sci = WCS(self.data_header)
            coords_proj = wcs_sci.world_to_pixel_values(coords[0], coords[1])

            return coords_proj

    def get_world_coords(self, coords: list):
        """Returns the coordinates in degrees (R.A., Dec.) of one or more pixel positions in the image.

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
            wcs_sci = WCS(self.data_header)
            coords_proj = wcs_sci.pixel_to_world_values(coords[0], coords[1])

            return coords_proj

    def crop(self, width: int, height: int, center: list[float] = None):
        """Cut out a portion of the image based on width and height

        Parameters
        -----------
        width : int
            Width of final image, in pixel.
        height = int
            Height of final image, in pixel.
        center : list, optional
            The central position of the final image in the reference frame of the 
            initial image, in pixels. It must be in the form [x,y].
        Returns
        --------
        Image object
            An image object with modified header considering cropping parameters. 
        """

        warnings.filterwarnings("ignore")

        data_cropped = np.copy(self.data)

        cx, cy = self.size[1] // 2 , self.size[0] // 2

        if center != None:
            cx, cy = int(center[0]), int(center[1])

        data_cropped = data_cropped[cy - height//2 : cy + height//2 , cx - width//2 : cx + width//2]
        data_cropped_size = np.shape(data_cropped)


        wcs = WCS(self.data_header)
        new_wcs = wcs.deepcopy()
        wcs_matrix = wcs.wcs.pc 

        x_refpx, y_refpx = float(self.data_header['CRPIX1']), float(self.data_header['CRPIX2'])     # In the initial image

        x_refpx_new, y_refpx_new = width // 2, height // 2                                          # In the cropped image
        x_refdeg_new, y_refdeg_new = wcs.pixel_to_world_values(cx, cy)

        new_wcs.wcs.crpix = [x_refpx_new+1, y_refpx_new+1]
        new_wcs.wcs.crval = [x_refdeg_new, y_refdeg_new]

        new_data_header = self.data_header.copy()
        new_data_header.update(new_wcs.to_header())
        new_data_header['NAXIS1'] = data_cropped_size[1]
        new_data_header['NAXIS2'] = data_cropped_size[0]

        cropped_image = Image.from_file_extension(self.primary_header, new_data_header, data_cropped)

        return cropped_image

    def rotate(self, angle: float, control_plot: bool = False):
        """Rotates the image by modifying the WCS of the file headers.

        Parameters
        -----------
        angle : float
            Angle of rotation to be applied to data. The angle follows 
            the counter-clockwise convention.
        control_plot : float, optional
            If True, show the image before and after rotation.

        Returns
        ---------
        Image object
            Image rotated, with headers updated.
        """

        wcs = WCS(self.data_header)

        # Rotation matrix definition
        angle_radian = np.radians(angle)

        # Counter-clockwise rotation 
        rotation_matrix = np.array([[np.cos(angle_radian),  np.sin(angle_radian)], 
                                    [-np.sin(angle_radian), np.cos(angle_radian)]])         

        wcs_rotated = wcs.deepcopy()
        wcs_rotated.wcs.pc = np.dot(rotation_matrix, wcs.wcs.pc)

        # Update header with new WCS information
        data_header_rotated = self.data_header.copy()
        data_header_rotated.update(wcs_rotated.to_header())

        # Rotate image without changing pixel size
        rotated_image = rotate(self.data, angle, reshape=False, order=1, mode='nearest')


        if control_plot:

            fig, axs = plt.subplots(1,2)

            axs[0].imshow(abs(self.data), cmap='inferno', origin='lower', norm=colors.LogNorm())
            axs[1].imshow(abs(rotated_image), cmap='inferno', origin='lower', norm=colors.LogNorm())

            axs[0].set_title('Before rotation')
            axs[1].set_title('After rotation: ${\\theta} = $' + '{}'.format(angle) + '$^\degree$')

            fig.tight_layout()
            #fig.savefig('check_rotation.png', dpi=300)
            plt.show()

        return Image.from_file_extension(self.primary_header, data_header_rotated, rotated_image)

    def convolve(self, fwhm: float, psf: str = 'gaus', control_plot: bool = False):
        """Convolves the image with a convolution kernel.

        Parameters
        -----------
        fwhm : float
            Full-width at half-maximum of convolution profile. Whatever the convolution 
            kernel, the 2D profile is symmetrical. The value must be given in arcsec.
        psf : str, optional
            Convolution kernel. Possible choices are: 'gaussian', 'voigt', 'lorentz'.
        control_plot: bool, optional    
            If True, shows the image before and after convolution.
        Returns
        -----------
        Image object
            The initial image convoluted by a PSF profile.
        """

        all_psf = ['gaussian', 'voigt', 'lorentz']

        fwhm_px = fwhm / self.px_size                           # arcsec into pixel conversion

        convolved_image = np.full(self.size, np.nan)
        
        if psf == 'gaussian':

            sigma_px = fwhm_px / (2 * np.sqrt(2 * np.log(2)))   # Gaussian width
            convolved_image = gaussian_filter(self.data, sigma=sigma_px)

        elif psf == 'voigt':                                    # Gaussian profile convolved with Lorentz profile

            x_px = fwhm_px / 1.63759                            # Voigt FWHM such as FWHM_G = FWHM_L
            sigma_gaus_px = x_px / (2 * np.sqrt(2 * np.log(2))) # Gaussian width
            gamma_lorentz_px = x_px / 2                         # Lorentz width

            # Kernel parameters
            kernel_radius = int(np.ceil(3 * max(sigma_gaus_px, gamma_lorentz_px)))
            kernel_size = 2 * kernel_radius + 1

            x_values = np.linspace(-kernel_radius, kernel_radius, kernel_size)

            # 1D Voigt profile
            voigt_profile_1d = voigt_profile(x_values, sigma_gaus_px, gamma_lorentz_px)
            # 2D Voigt profile
            voigt_kernel_2d = np.outer(voigt_profile_1d, voigt_profile_1d)
            # Profile Normalization
            voigt_kernel_2d /= np.sum(voigt_kernel_2d)
            # 2D Convolution
            convolved_image = fftconvolve(self.data, voigt_kernel_2d, mode='same')

        elif psf == 'lorentz':

            gamma_lorentz_px = fwhm_px / 2                      # Lorentz width

            # Kernel parameters
            kernel_radius = int(np.ceil(3*gamma_lorentz_px))
            kernel_size = 2 * kernel_radius + 1

            x_values = np.linspace(-kernel_radius, kernel_radius, kernel_size)

            # 1D lorentz profile | Formula : L(x) = (1/pi) * (gamma / (x^2 + γ^2))
            lorentz_1d = (1 / np.pi) * (gamma_lorentz_px / (x_values**2 + gamma_lorentz_px**2))

            # 2D Lorentz profile
            lorentz_kernel_2d = np.outer(lorentz_1d, lorentz_1d)
            # Profile normalization
            lorentz_kernel_2d /= np.sum(lorentz_kernel_2d)
            # 2D Convolution
            convolved_image = fftconvolve(self.data, lorentz_kernel_2d, mode='same')


        if control_plot:

            fig, axs = plt.subplots(1,2)

            axs[0].imshow(self.data, origin='lower', cmap='inferno')
            axs[1].imshow(convolved_image, origin='lower', cmap='inferno')

            ny_image, nx_image = convolved_image.shape
            ellipse_center = (nx_image * 0.1, ny_image * 0.9)

            psf_ellipse = Ellipse(xy=ellipse_center, width=fwhm_px, height=fwhm_px, 
                    edgecolor='white', facecolor='none', lw=1)

            axs[1].add_patch(psf_ellipse)
            axs[1].annotate('FWHM PSF', xy=ellipse_center,
                 xytext=(ellipse_center[0] + fwhm_px, ellipse_center[1]),
                 arrowprops=dict(facecolor='red', arrowstyle='->'),
                 color='white', verticalalignment='center')

            axs[0].set_title('Originale image')
            axs[1].set_title('Convolved image')

            fig.tight_layout()
            plt.show()


        return Image.from_file_extension(self.primary_header, self.data_header, convolved_image)







    