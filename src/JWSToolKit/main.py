from JWSToolKit.Cube import Cube
from JWSToolKit.Spec import Spec
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from Cube import Cube
from Spec import Spec

file = "/Users/delabrov/Desktop/obs_finales/NIRSpec-JWST/DGTauB/jw01644_nirspec_g140h-f100lp_s3d.fits"
cube = Cube(file)

wvs_values = cube.get_wvs()

cube.info()

#int_map = cube.line_emission_map(wv_line = 1.64355271, map_units='erg s-1 cm-2 sr-1', control_plot=False)


WV_LINE = 1.64355271

spectrum_values = cube.extract_spec_circ_aperture(radius=4, position=(25,25), units='Jy')

spectrum = Spec(wvs_values, spectrum_values, units='Jy')

spectrum_cutted = spectrum.cut(-2000, 2000, units='vel', wv_ref=WV_LINE)

spectrum_baseline_sub = spectrum_cutted.sub_baseline(wv_line=WV_LINE)

integrated_intensity = spectrum_baseline_sub.line_integ(wv_line=WV_LINE, profile='gaus', control_plot=True)

doppler_shift_line = spectrum_baseline_sub.line_velocity(wv_line=WV_LINE, control_plot=False)




















"""
spec = cube.extract_spec_circ_aperture(4, [27,27], units='Jy')      # Spectrum extraction

fig, ax = plt.subplots(figsize=(9,5))
ax.step(wvs, spec, color='black')
ax.set_xlabel('Wavelength (µm)')
ax.set_ylabel('Flux density (Jy)')
#plt.show()
plt.close()

spectrum = Spec(wvs, spec, units='Jy')
spectrum.convert(units='erg s-1 cm-2 um-1')                         # Conversion
spectrum_red = spectrum.cut(-2000, 2000, units='vel', wv_ref=2.12)
spectrum_baseline_sub = spectrum_red.sub_baseline(wv_line=2.1218, control_plot=False)
flux_line, err_flux_line = spectrum_baseline_sub.line_integ(wv_line=2.1218, profile='gaus', control_plot=False)
vel_line, err_vel_line = spectrum_baseline_sub.line_velocity(wv_line=2.1218, control_plot=True)
"""
