"""
The Spec object stores information about spectra that can be extracted 
from JWST data cubes. The spectrum is in the form of two 1D lists containing 
wavelength values in µm and surface brightness (flux per unit 'angular area') 
or flux values. By default, JWST data are given in MJy/sr, i.e. surface 
brightness. However, other units can be specified for the y-axis of the spectrum. 

The three preceding elements: wavelength axis, value axis and units are 
attributes of the Spec object. The increment value on the wavelength axis 
is the object's fourth argument. 


Parameters
----------
wvs : array_like
    The x-axis of the spectrum. The axis corresponds to wavelengths in 
    µm. This is the default unit for JWST data.
values : array_like 
    The y axis of the spectrum. By default in JWST data, the values in 
    each spaxel are surface brightness, given in MJy/sr. It is possible 
    to give other units, in flux or surface brightness. The number of 
    points on the y axis must be identical to the number of points on the 
    wavelength axis. 
units : str, optional
    Unit of the y axis of the spectrum. Default values in JWST data are 
    surface brightness in MJy/sr. 

Attributes
----------
wvs : array_like 
    The wavelength values of the spectrum. 
values : array_like 
    The brightness surface or flux values of the spectrum. 
    The unit must match that of the 'units' argument.
units : str
    The unit of the y axis values of the spectrum. Possible units 
    are: 'MJy/sr', 'Jy', 'erg s-1 cm-2 Hz-1', 'erg s-1 cm-2 um-1', 
    'erg s-1 cm-2 um-1 sr-1'.
dwvs : float
    Interval between two points on the wavelength axis, in µm.
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.optimize import curve_fit
from scipy.special import wofz

c_sp = 299792458            # Speed of light (m/s)

class Spec:
    def __init__(self, wvs, values, units='MJy/sr'):

        if wvs.shape[0] != values.shape[0]:
            raise AssertionError('The number of wavelength elements must be identical to the number of spectral values.')
        else:
            self.wvs = wvs                              # Wavelength grid (µm)
            self.values = values                        # Spectrum values
            self.units = units                          # Spectrum values units
            self.dwvs = np.nanmean(np.diff(self.wvs))   # Wavelength step size (µm)

    def convert(self, units, px_area=1):
        """Convert spectrum values into another unit.

        Parameters
        ----------
        units : str
            The unit of the spectrum after conversion. Possible units are : 
            MJy/sr, Jy, erg s-1 cm-2 Hz-1, erg s-1 cm-2 um-1, erg s-1 cm-2 um-1 sr-1.
        px_area : float, optional
            Spatial area over which the spectrum was extracted. Used to 
            convert surface brightness to flux density. The value must 
            be given in steradian.
            
        Returns
        ----------
        Spec 
            The initial spectrum with the y-axis converted to the chosen unit. 
        """

        all_units = ['MJy/sr', 'Jy', 'erg s-1 cm-2 Hz-1', 'erg s-1 cm-2 um-1', 'erg s-1 cm-2 um-1 sr-1']

        if not isinstance(units, str) or units not in all_units:
            raise TypeError('The input units are invalid. They must be one of the following units: MJy/sr, Jy, erg s-1 cm-2 Hz-1, erg s-1 cm-2 um-1, erg s-1 cm-2 um-1 sr-1')
        elif not isinstance(px_area, (int, float)) or px_area < 0:
            raise AssertionError('The input pixel area is invalid. It must be a positive float in steradian.')
        else:

            if units == self.units:                     # Input units as same as self.units
                new_values = np.copy(self.values)

            elif self.units == all_units[0]:            # MJy/sr                                ! Surface Brightness !
                new_values = np.copy(self.values)

                if units == all_units[1]:               # -> Jy                                 ! Flux Density F_nu !
                    new_values *= 1e6 * px_area

                elif units == all_units[2]:             # -> erg s-1 cm-2 Hz-1                  ! Flux Density F_nu !
                    new_values *= 1e6 * px_area * 1e-23

                elif units == all_units[3]:             # -> erg s-1 cm-2 micron-1              ! Flux Density F_lbd !
                    new_values *= 1e6 * px_area
                    new_values *= 1e-23
                    new_values *= (c_sp*1e6 / (self.wvs)**2)

                elif units == all_units[4]:             # -> erg s-1 cm-2 micron-1 sr-1         ! Surface Brightness !
                    new_values *= 1e6
                    new_values *= 1e-23
                    new_values *= (c_sp*1e6 / (self.wvs)**2)

            elif self.units == all_units[1]:        # Jy                                        ! Flux Density F_nu !
                new_values = np.copy(self.values)

                if units == all_units[0]:           # -> MJy/sr                                 ! Surface Brightness !
                    new_values /= (1e6 * px_area)

                elif units == all_units[2]:         # -> erg s-1 cm-2 Hz-1                      ! Flux Density F_nu !
                    new_values *= 1e-23

                elif units == all_units[3]:         # -> erg s-1 cm-2 micron-1                  ! Flux Density F_lbd !
                    new_values *= 1e-23
                    new_values *= (c_sp * 1e6 / (self.wvs) ** 2)

            elif self.units == all_units[2]:        # erg s-1 cm-2 Hz-1                         ! Flux Density F_nu !
                new_values = np.copy(self.values)

                if units == all_units[0]:           # -> MJy/sr                                 ! Surface Brightness !
                    new_values /= (1e-23 * 1e6 * px_area)

                elif units == all_units[1]:         # -> Jy                                     ! Flux Density F_nu !
                    new_values /= 1e-23

                elif units == all_units[3]:         # -> erg s-1 cm-2 micron-1                  ! Flux Density F_lbd !
                    new_values *= (c_sp * 1e6 / (self.wvs) ** 2)

            elif self.units == all_units[3]:        # erg s-1 cm-2 micron-1                     ! Flux Density F_lbd !
                new_values = np.copy(self.values)

                if units == all_units[0]:           # -> MJy/sr                                 ! Surface Brightness !
                    new_values /= (c_sp * 1e6 / (self.wvs) ** 2)
                    new_values /= (1e-23 * 1e6 * px_area)

                elif units == all_units[1]:         # -> Jy                                     ! Flux Density F_nu !
                    new_values /= (c_sp * 1e6 / (self.wvs) ** 2)
                    new_values /= 1e-23

                elif units == all_units[2]:         # -> erg s-1 cm-2 Hz-1                      ! Flux Density F_nu !
                    new_values /= (c_sp * 1e6 / (self.wvs) ** 2)

            self.values = new_values
            self.units  = units

        return

    def cut(self, min, max, units='wav', wv_ref=None):
        """Extract a part of the spectrum using limits of the requiered interval.

        Parameters 
        ----------
        min : float 
            Value of the lower limit of the cut spectrum interval. If no unit is specified or the specified unit is 'wav', then the value is a wavelength in µm.
        max : float
            Value of the upper limit of the cut spectrum interval. If no unit is specified or the specified unit is 'wav', then the value is a wavelength in µm. 
        units : str, optional
            Unit of the x-axis of the spectrum. If the unit is 'wav', then the x-axis values are wavelengths, and the interval limits must be given in wavelengths. If the unit is 'vel', then the values of the x-axis are radial velocities, and the limits of the interval must be given in km/s. The reference wavelength must then be specified to convert wavelengths into radial velocities (via the Doppler shift relation).
        wv_ref : Optional
            In the case of the 'vel' unit, a wavelength reference value must be given for conversion to radial velocity (via the Doppler shift relation), in µm.

        Returns
        ----------
        Spec 
            A Spec object cut by considering the limits of the interval given as input parameters.
        """

        all_units = ['wav', 'vel']

        if not isinstance(units, str) or units not in all_units:
            raise AssertionError("The units input are invalid. They must correspond to wavelengths or velocities. The following character strings must be specified: 'wav' or 'vel'.")

        elif units == all_units[1]:
            if wv_ref is None:
                raise AssertionError("Bounds are given in velocities. No reference wavelength has been specified for the 'wv_ref' parameter.")
            elif not isinstance(wv_ref, (int, float)) or wv_ref <= 0:
                raise AssertionError("The reference wavelength is invalid. It must be an integer float and given in km/s.")
            else:
                rvs = (c_sp * (self.wvs - wv_ref) / wv_ref) / 1000  # km/s
                idxs_cut = np.where(np.logical_and(rvs >= min, rvs <= max))
                wvs_cut = self.wvs[idxs_cut]
                values_cut = self.values[idxs_cut]
                spec_cut = Spec(wvs_cut, values_cut, units=self.units)

        else:
            if units == all_units[0]:
                idxs_cut = np.where(np.logical_and(self.wvs >= min, self.wvs <= max))
                wvs_cut = self.wvs[idxs_cut]
                values_cut = self.values[idxs_cut]
                spec_cut = Spec(wvs_cut, values_cut, units=self.units)

        return spec_cut

    def sub_baseline(self, wv_line, mask_rv=200, deg=1, control_plot=False):
        """Subtracts a baseline by fitting a polynomial around an emission line. 

        Parameters 
        ----------
        wv_line : float
             Wavelength in vacuum and at rest of the emission line, in µm.
        mask_rv : float, optional
            Half-width of the interval used to exclude spectral pixels in baseline fitting. Can be interpreted as line half-width. The value should be given in km/s. 
        deg : int, optional
            Degree of the polynomial used to adjust the baseline around the line. 
        control_plot : bool, optional
            If True, show the different stages of the baseline subtraction.

        Returns
        ----------
        Spec 
           The initial spectrum, subtracted from its baseline.
        """

        dwvs = np.nanmean(np.diff(self.wvs))

        rvs = (c_sp * (self.wvs - wv_line) / wv_line) / 1000  # km/s

        idx_max = np.where(self.values == np.nanmax(self.values))
        v_max = rvs[idx_max]
        v_min_line, v_max_line = v_max - mask_rv, v_max + mask_rv

        idxs_line = np.where(np.logical_and(rvs >= v_min_line, rvs <= v_max_line))
        idxs_cont = np.concatenate([np.arange(0,idxs_line[0][0]+1), np.arange(idxs_line[0][-1], np.shape(self.values)[0])])

        wv_min, wv_max= self.wvs[idxs_line[0][0]], self.wvs[idxs_line[0][-1]]

        wvs_cont = self.wvs[idxs_cont]
        spec_cont = self.values[idxs_cont]

        params = np.polyfit(wvs_cont, spec_cont, deg=deg)
        baseline = np.poly1d(params)(self.wvs)

        spec_cont_sub = self.values - baseline

        if control_plot:
            fig, ax = plt.subplots(figsize=(9,4))

            ax.step(self.wvs+self.dwvs/2, self.values, color='blue', label='Not subtracted')
            ax.step(self.wvs+self.dwvs/2, baseline, color='red', label='baseline')
            ax.step(self.wvs+self.dwvs/2, spec_cont_sub, color='black', label='subtracted')

            ax.axvline(wv_min, color='green', linestyle='--', label='masking limits')
            ax.axvline(wv_max, color='green', linestyle='--')

            ax.set_xlabel('Wavelength (µm)')

            ax.legend()

            fig.tight_layout()
            plt.show()

        return Spec(self.wvs, spec_cont_sub, units=self.units)

    def line_integ(self, wv_line, profile=None, line_width=400, control_plot=False):
        """Computes the integrated flux or surface brightness of a line.  

        Parameters 
        ----------
        wv_line : float
            Wavelength in vacuum and at rest of the emission line, in µm.
        profile : str, optional
            Type of line profile. If no profile is specified, integration is performed using a simpson method. If a profile is specified, it must be one of the following: 'gaus', 'voigt', 'lorentz' and 'moffat'. 
        line_width : float, optional 
            Full-width of the line. The value should be given in km/s.
        control_plot : bool, optional
            If True, show the different steps of the line integration.

        Returns
        ----------
        float, float 
            The first value corresponds to the flux or integrated surface brightness of the line, the second to the error on the first value. 
        """

        if profile == 'gaus':

            def gaussian(x, a, x0, sigma):
                return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

            # Line fit
            p0 = [np.nanmax(self.values), wv_line, 1e-3]
            popt, pcov = curve_fit(gaussian, self.wvs, self.values, p0=p0)

            flux_int = np.sqrt(2*np.pi) * popt[0] * popt[2]
            err_flux_int = flux_int * (np.sqrt(pcov[0][0] / (popt[0]) ** 2 + pcov[2][2] / (popt[2]) ** 2))

            # Plot
            if control_plot:

                x_new = np.linspace(self.wvs[0], self.wvs[-1], np.shape(self.wvs)[0]*4)
                line_fit = gaussian(x_new, *popt)

                fig, ax = plt.subplots(figsize=(8,4))

                ax.step(self.wvs+self.dwvs/2, self.values, c='black', label='Data')
                ax.plot(x_new, line_fit, c='red', label='Best Gaussian fit')
                ax.fill_between(x_new, line_fit, 0, color='red', alpha=0.3)

                ax.legend()

                fig.tight_layout()
                plt.show()

            return flux_int, err_flux_int

        elif profile == 'voigt':

            def voigt(x, A, x0, gamma, sigma):
                return A * np.real(wofz((x - x0 + 1j*gamma) / sigma / np.sqrt(2))) / sigma / np.sqrt(2*np.pi)

            # Error on continuum
            idx_peak_line = np.where(self.values == np.nanmax(self.values))[0]
            rvs = (c_sp * (self.wvs - wv_line) / wv_line) / 1000  # km/s
            rv_peak_line = rvs[idx_peak_line]
            idxs_line = np.where(np.logical_and(rvs >= rv_peak_line - line_width/2, rvs <= rv_peak_line + line_width/2))

            n_bins = np.shape(idxs_line)[1]

            idxs_cont_l = np.where(rvs <= rv_peak_line - line_width/2)
            idxs_cont_r = np.where(rvs >= rv_peak_line + line_width/2)
            idxs_cont = np.concatenate([idxs_cont_l[0], idxs_cont_r[0]])
            std_cont = np.nanstd(self.values[idxs_cont])

            # Line fit
            p0 = [np.nanmax(self.values), wv_line, 1e-3, 1e-3]
            popt, pcov = curve_fit(voigt, self.wvs, self.values, p0=p0, sigma=[std_cont]*np.shape(self.wvs)[0])

            wvs_new = np.linspace(self.wvs[0], self.wvs[-1], np.shape(self.wvs)[0]*4)
            line_fit = voigt(wvs_new, *popt)

            flux_int = integrate.simpson(line_fit, wvs_new)

            # Error on flux value
            err_flux_int = 0

            # Plot
            if control_plot:

                x_new = np.linspace(self.wvs[0], self.wvs[-1], np.shape(self.wvs)[0]*4)
                line_fit = voigt(x_new, *popt)

                fig, ax = plt.subplots(figsize=(8,4))

                ax.step(self.wvs+self.dwvs/2, self.values, c='black', label='Data')
                ax.plot(x_new, line_fit, c='red', label='Best Voigt fit')
                ax.fill_between(x_new, line_fit, 0, color='red', alpha=0.3)

                ax.legend()

                fig.tight_layout()
                plt.show()


            return flux_int, err_flux_int

        elif profile == 'lorentz':

            def lorentz(x, A, x0, gamma):
                return A * gamma / np.pi / ((x - x0)**2 + gamma**2)

            # Error on continuum
            idx_peak_line = np.where(self.values == np.nanmax(self.values))[0]
            rvs = (c_sp * (self.wvs - wv_line) / wv_line) / 1000  # km/s
            rv_peak_line = rvs[idx_peak_line]
            idxs_line = np.where(
                np.logical_and(rvs >= rv_peak_line - line_width / 2, rvs <= rv_peak_line + line_width / 2))

            n_bins = np.shape(idxs_line)[1]

            idxs_cont_l = np.where(rvs <= rv_peak_line - line_width / 2)
            idxs_cont_r = np.where(rvs >= rv_peak_line + line_width / 2)
            idxs_cont = np.concatenate([idxs_cont_l[0], idxs_cont_r[0]])

            std_cont = np.nanstd(self.values[idxs_cont])

            # Line Fit
            p0 = [np.nanmax(self.values), wv_line, 1e-3]
            popt, pcov = curve_fit(lorentz, self.wvs, self.values, p0=p0, sigma=[std_cont]*np.shape(self.wvs)[0])

            wvs_new = np.linspace(self.wvs[0], self.wvs[-1], np.shape(self.wvs)[0] * 4)
            line_fit = lorentz(wvs_new, *popt)

            # Line integration
            flux_int = integrate.simpson(line_fit, wvs_new)

            # Error on flux value
            err_flux_int = 0


            if control_plot:

                fig, ax = plt.subplots(figsize=(8,4))

                ax.step(self.wvs+self.dwvs/2, self.values, c='black', label='Data')
                ax.plot(wvs_new, line_fit, c='red', label='Best Lorentz fit')
                ax.fill_between(wvs_new, line_fit, 0, color='red', alpha=0.3)

                ax.legend()

                fig.tight_layout()
                plt.show()


            return flux_int, err_flux_int

        elif profile == 'moffat':

            def moffat(x, A, x0, alpha, beta):
                return A * (1 + ((x-x0) / alpha)**2 )**(-beta)

            # Error on continuum
            idx_peak_line = np.where(self.values == np.nanmax(self.values))[0]
            rvs = (c_sp * (self.wvs - wv_line) / wv_line) / 1000  # km/s
            rv_peak_line = rvs[idx_peak_line]
            idxs_line = np.where(
                np.logical_and(rvs >= rv_peak_line - line_width / 2, rvs <= rv_peak_line + line_width / 2))

            n_bins = np.shape(idxs_line)[1]

            idxs_cont_l = np.where(rvs <= rv_peak_line - line_width / 2)
            idxs_cont_r = np.where(rvs >= rv_peak_line + line_width / 2)
            idxs_cont = np.concatenate([idxs_cont_l[0], idxs_cont_r[0]])

            std_cont = np.nanstd(self.values[idxs_cont])


            # Line fit
            p0 = [np.nanmax(self.values), wv_line, 1e-3, 1]
            popt, pcov = curve_fit(moffat, self.wvs, self.values, p0=p0, sigma=[std_cont]*np.shape(self.wvs)[0])

            wvs_new = np.linspace(self.wvs[0], self.wvs[-1], np.shape(self.wvs)[0] * 4)
            line_fit = moffat(wvs_new, *popt)

            # Line integration
            flux_int = integrate.simpson(line_fit, wvs_new)

            # Error on flux value
            err_flux_int = 0

            # Plot
            if control_plot:

                fig, ax = plt.subplots(figsize=(8, 4))

                ax.step(self.wvs + self.dwvs / 2, self.values, c='black', label='Data')
                ax.plot(wvs_new, line_fit, c='red', label='Best Moffat fit')
                ax.fill_between(wvs_new, line_fit, 0, color='red', alpha=0.3)

                ax.legend()

                fig.tight_layout()
                plt.show()

            return flux_int, err_flux_int

        elif profile is None:

            # Error on continuum
            idx_peak_line = np.where(self.values == np.nanmax(self.values))[0]
            rvs = (c_sp * (self.wvs - wv_line) / wv_line) / 1000  # km/s
            rv_peak_line = rvs[idx_peak_line]
            idxs_line = np.where(
                np.logical_and(rvs >= rv_peak_line - line_width / 2, rvs <= rv_peak_line + line_width / 2))

            n_bins = np.shape(idxs_line)[1]

            idxs_cont_l = np.where(rvs <= rv_peak_line - line_width / 2)
            idxs_cont_r = np.where(rvs >= rv_peak_line + line_width / 2)
            idxs_cont = np.concatenate([idxs_cont_l[0], idxs_cont_r[0]])

            std_cont = np.nanstd(self.values[idxs_cont])

            wvs_line = self.wvs[idxs_line]
            spec_line = self.values[idxs_line]

            # Line integration
            flux_int = integrate.simpson(spec_line, wvs_line)

            # Error on flux value
            err_flux_int = 0

            # Plot
            if control_plot:
                fig, ax = plt.subplots(figsize=(8, 4))

                ax.step(self.wvs + self.dwvs / 2, self.values, c='black', label='Data')
                ax.axvline(np.nanmin(wvs_line), color='red', linestyle='--', label='Integration interval')
                ax.axvline(np.nanmax(wvs_line), color='red', linestyle='--')
                ax.fill_between(self.wvs, self.values, 0, where = (self.wvs >= np.nanmin(wvs_line)) & (self.wvs <= np.nanmax(wvs_line)), color='red', alpha=0.3)

                ax.legend()

                fig.tight_layout()
                plt.show()

            return flux_int, err_flux_int

        return

    def line_velocity(self, wv_line, line_width=400, control_plot=False):
        """Calculates the Doppler shift of an emission line using Gaussian profile fitting.

        Parameters 
        ----------
        wv_line : float
            Wavelength in vacuum and at rest of the emission line, in µm.
        line_width : float, optional 
            Full-width of the line. The value should be given in km/s.
        control_plot : bool, optional
            If True, show the line fitting.

        Returns
        ----------
        float, float 
            The first value corresponds to the Doppler shift of the line, in terms of radial velocity given in km/s, the second the error associated with the first value. 
        """

        def gaussian(x, a, x0, sigma):
            return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

        rvs = (c_sp * (self.wvs - wv_line) / wv_line) / 1000  # km/s
        drvs = np.nanmean(np.diff(rvs))

        # Error on continuum
        idx_peak_line = np.where(self.values == np.nanmax(self.values))[0]
        rv_peak_line = rvs[idx_peak_line]
        idxs_line = np.where(
            np.logical_and(rvs >= rv_peak_line - line_width / 2, rvs <= rv_peak_line + line_width / 2))

        n_bins = np.shape(idxs_line)[1]

        idxs_cont_l = np.where(rvs <= rv_peak_line - line_width / 2)
        idxs_cont_r = np.where(rvs >= rv_peak_line + line_width / 2)
        idxs_cont = np.concatenate([idxs_cont_l[0], idxs_cont_r[0]])

        std_cont = np.nanstd(self.values[idxs_cont])

            # Line fit
        p0 = [np.nanmax(self.values), 0, line_width/2]
        popt, pcov = curve_fit(gaussian, rvs, self.values, p0=p0, sigma=[std_cont]*np.shape(self.wvs)[0])

        if control_plot:

            x_new = np.linspace(rvs[0], rvs[-1], np.shape(rvs)[0]*4)
            line_fit = gaussian(x_new, *popt)

            fig, ax = plt.subplots(figsize=(8,4))

            ax.step(rvs + drvs/2, self.values, c='black', label='Data')
            ax.plot(x_new, line_fit, c='red', label='Best Gaussian fit')
            ax.fill_between(x_new, line_fit, 0, color='red', alpha=0.3)

            ax.legend()

            fig.tight_layout()
            plt.show()

        return popt[1], np.sqrt(pcov[1][1])

    def copy(self):
        """Copy a spectrum.

        Parameters 
        ----------

        Returns
        ----------
        Spec
            The initial spectrum copied.
        """

        new_wvs = np.copy(self.wvs)
        new_values = np.copy(self.values)
        return Spec(new_wvs, new_values, units=self.units)




"""
- convolve : 
    Convolution du spectre avec un profil Gaussien, changer la résolution spectral
- fwhm : 
    Retourne la largeur d'une raie d'émission
- plot : 
    Affiche le spectre
"""



"""
def sum(spec1, spec2, units=None, px_area=1):

    all_units = ['MJy/sr', 'Jy', 'erg s-1 cm-2 Hz-1', 'erg s-1 cm-2 um-1']

    wvs1, wvs2 = spec1.get_wvs(), spec2.get_wvs()
    units1, units2 = spec1.get_units(), spec2.get_units()
    values1, values2 = spec1.get_values(), spec2.get_values()
    values1 = np.array(values1)
    values2 = np.array(values2)

    if wvs1.all() == wvs2.all():                # Same Wavelength grid

        if units1 == units2:        # Same units
            new_values = values1 + values2

            if units != None:
                if not isinstance(units, str) or units not in all_units:
                    raise TypeError(
                        'The input units are invalid. They must be one of the following units: MJy/sr, Jy, erg s-1 cm-2 Hz-1, erg s-1 cm-2 micron-1.')
                else:
                    new_spec = Spec(wvs1, new_values, units=units1).copy()
                    new_spec.convert(units=units, px_area=px_area)
                    return new_spec
            else:
                new_spec = Spec(wvs1, new_values, units=units1)
                return new_spec

        else:                       # Different units

            if units == None:
                print('The units of the two spectra are different. You must specify units for the output spectrum.')
            else:

                if units1 != units:
                    spec1_convert = spec1.copy()
                    spec1.convert(units=units, px_area=px_area)
                if units2 != units:
                    spec2_convert = spec2.copy()
                    spec2.convert(units=units, px_area=px_area)

                new_values = np.array(spec1_convert.get_values()) + np.array(spec2_convert.get_values())
                new_spec = Spec(wvs1, new_values, units=units)
                return new_spec


    else:                       # Different wavelength grid
        print("The spectra do not have the same wavelength grid.")

"""
