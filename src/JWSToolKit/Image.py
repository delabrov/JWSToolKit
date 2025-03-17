import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import rotate
import matplotlib.colors as colors
from tqdm import tqdm


class Image:

    def __init__(self, file_name):

        if not isinstance(file_name, str):
            raise TypeError("The input file name is invalid. It must be a character string")
        else:

            self.file_name = file_name                              # Image name
            self.primary_header, self.data_header, self.data, self.errs = self._load_fits(file_name)

            self.size = self.data.shape                             # Data cube size
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


    def _load_fits(self, file_name):

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
            lims=None, 
            abs_transform: bool = False, 
            save: bool = False, 
            colorbar: bool = False,
            origin_arcsec=None):

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
            ax.grid(color='grey', ls='--')
            ax.set_xlabel('Right Ascension (RA)')
            ax.set_ylabel('Declination (Dec)')

        if origin_arcsec != None:
            ax.set_xlabel('$\Delta$X (arcsec)')
            ax.set_ylabel('$\Delta$Y (arcsec)')

        fig.tight_layout()
        if save:
            fig.savefig('image.png', dpi=300)
        plt.show()
