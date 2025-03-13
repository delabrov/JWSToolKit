from JWSToolKit.Cube import Cube
from JWSToolKit.Spec import Spec
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

file = "/Users/delabrov/Desktop/obs_finales/NIRSpec-JWST/DGTauB/jw01644_nirspec_g140h-f100lp_s3d.fits"
cube = Cube(file)

wvs = cube.get_wvs()

cube.info()

int_map = cube.line_emission_map(wv_line = 1.64355271, map_units='erg s-1 cm-2 sr-1', control_plot=False)




# !!! Mettre toutes les méthodes sous le format EstUneFonction !!!
# !!! Ajouter les annotations de type dans les méthodes !!!



























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
