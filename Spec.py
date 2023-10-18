import numpy as np



c_sp = 299792458            # Speed of light (m/s)

class Spec:

    def __init__(self, wvs, values, units='MJy/sr'):

        if wvs.shape[0] != values.shape[0]:
            raise AssertionError('The number of wavelength elements must be identical to the number of spectral values.')
        else:
            self.wvs = wvs                      # Wavelength grid (µm)
            self.values = values                # Spectrum values
            self.units = units                  # Spectrum values units

    def convert(self, units, px_area=1):

        all_units = ['MJy/sr', 'Jy', 'erg s-1 cm-2 Hz-1', 'erg s-1 cm-2 micron-1']

        if not isinstance(units, str) or units not in all_units:
            raise TypeError('The input units are invalid. They must be one of the following units: MJy/sr, Jy, erg s-1 cm-2 Hz-1, erg s-1 cm-2 micron-1.')
        elif not isinstance(px_area, (int, float)) or px_area < 0:
            raise AssertionError('The input pixel area is invalid. It must be a positive float in steradian.')
        else:

            if units == self.units:                     # Input units as same as self.units
                pass

            elif self.units == all_units[0]:            # MJy/sr
                new_values = np.copy(self.values)

                if units == all_units[1]:               # -> Jy
                    new_values *= 1e6 * px_area

                elif units == all_units[2]:             # -> erg s-1 cm-2 Hz-1
                    new_values *= 1e6 * px_area * 1e-23

                elif units == all_units[3]:             # -> erg s-1 cm-2 micron-1
                    new_values *= 1e6 * px_area
                    new_values *= 1e-23
                    new_values *= (c_sp*1e6 / (self.wvs)**2)

                elif self.units == all_units[1]:        # Jy
                    new_values = np.copy(self.values)

                    if units == all_units[0]:           # -> MJy/sr
                        new_values /= (1e6 * px_area)

                    elif units == all_units[2]:         # -> erg s-1 cm-2 Hz-1
                        new_values *= 1e-23

                    elif units == all_units[3]:         # -> erg s-1 cm-2 micron-1
                        new_values *= 1e-23
                        new_values *= (c_sp * 1e6 / (self.wvs) ** 2)

                elif self.units == all_units[2]:        # erg s-1 cm-2 Hz-1
                    new_values = np.copy(self.values)

                    if units == all_units[0]:           # -> MJy/sr
                        new_values /= (1e-23 * 1e6 * px_area)

                    elif units == all_units[1]:         # -> Jy
                        new_values /= 1e-23

                    elif units == all_units[3]:         # -> erg s-1 cm-2 micron-1
                        new_values *= (c_sp * 1e6 / (self.wvs) ** 2)

                elif self.units == all_units[3]:        # erg s-1 cm-2 micron-1
                    new_values = np.copy(self.values)

                    if units == all_units[0]:           # -> MJy/sr
                        new_values /= (c_sp * 1e6 / (self.wvs) ** 2)
                        new_values /= (1e-23 * 1e6 * px_area)

                    elif units == all_units[1]:         # -> Jy
                        new_values /= (c_sp * 1e6 / (self.wvs) ** 2)
                        new_values /= 1e-23

                    elif units == all_units[2]:         # -> erg s-1 cm-2 Hz-1
                        new_values /= (c_sp * 1e6 / (self.wvs) ** 2)

            self.values = new_values





        return
    def cut(self, min, max, units='wav', wv_ref = None):

        all_units = ['wav', 'vel']

        if not isinstance(units, str) or units not in all_units:
            raise AssertionError("The units input are invalid. They must correspond to wavelengths or velocities. The following character strings must be specified: 'wav' or 'vel'.")

        elif units == all_units[1]:
            if wv_ref is None:
                raise AssertionError("Bounds are given in velocities. No reference wavelength has been specified for the 'wv_ref' parameter.")
            elif not isinstance(wv_ref, (int, float)) or wv_ref <= 0:
                raise AssertionError("The reference wavelength is invalid. It must be an integer float and given in km/s.")

        else:
            print('caca')
            if units == all_units[0]:
                idxs_cut = np.where(np.logical_and(self.wvs >= min, self.wvs <= max))
                wvs_cut = self.wvs[idxs_cut]
                values_cut = self.values[idxs_cut]

                return Spec(wvs_cut, values_cut)

            elif units == all_units[1]:

                rvs = (c_sp * (self.wvs - wv_ref) / wv_ref) / 1000          # km/s
                idxs_cut = np.where(np.logical_and(rvs >= min, rvs <= max))
                wvs_cut = self.wvs[idxs_cut]
                values_cut = self.values[idxs_cut]

                return Spec(wvs_cut, values_cut)


