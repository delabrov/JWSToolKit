import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import sys
import os 

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from JWSToolKit.Cube import Cube
from JWSToolKit.Image import Image
from JWSToolKit.Spec import Spec

DGTAUB_POSITION = [66.76071774, 26.09171944]

MIRI_image_file     = '/Users/delabrov/Documents/Data_Obs/JWST/MIRI/new_reduction/Images/pipeline_outputs/Image3_outputs/jw01644002001_02101_outputs3/jw01644-o002_t001_miri_f1800w_i2d.fits'
NIRCAM_image_file   = "/Users/delabrov/Desktop/obs_finales/NIRCam-JWST/jw01644004001_04101_nrcb_nircam_i2d.fits"
NIRSpec_cube_file   = '/Users/delabrov/Desktop/obs_finales/NIRSpec-JWST/DGTauB/jw01644005001_05101_nirspec_g235h-f170lp_s3d.fits'


cube = Cube(NIRSpec_cube_file)

spec = cube.extract_spec_circ_aperture(radius=4, position=[25,25], units='Jy')








image = Image(NIRCAM_image_file)
#image = Image(MIRI_image_file)
image.info()

x_px, y_px = image.get_px_coords(DGTAUB_POSITION)


# FIGURE
#image.plot(scale='sqrt', use_wcs=False, abs_transform=True, colorbar=True, lims=[0.5, 430], draw_compass=False)#, origin_arcsec=[x_px,y_px])
#plt.scatter(x_px, y_px, marker='*', color='white')
#plt.show()


image_cropped = image.crop(400, 400, center=[x_px, y_px+50])

#x_px, y_px = image_cropped.get_px_coords(DGTAUB_POSITION)
#image_cropped.plot(scale='sqrt', abs_transform='True', colorbar=False, lims=[0.5, 430], draw_compass=True)#, origin_arcsec=[x_px,y_px])
#plt.show()

#image_cropped.save_as_fits()

image_cropped_rotated = image_cropped.rotate(angle=65, control_plot=False)

#image_cropped_rotated.plot(scale='sqrt', abs_transform='True', lims=[0.5, 430], draw_compass=True)
#plt.show()

image_convolved = image_cropped.convolve(fwhm=0.8, psf='gaussian', control_plot=True)

image_convolved.plot(scale='sqrt', abs_transform='True', colorbar=False, draw_compass=True)#, origin_arcsec=[x_px,y_px])
plt.show()






#plt.imshow(abs(image_cropped), cmap='inferno', origin='lower')
#plt.show()





"""
WV_LINE = 2.121833725
ROT_ANGLE = 295

cube = Cube(file)
wvs_values = cube.get_wvs()
cube.info()
"""


"""
spectrum_values = cube.extract_spec_circ_aperture(radius=4, position=(25,25), units='Jy')
spectrum = Spec(wvs_values, spectrum_values, units='Jy')
spectrum_cutted = spectrum.cut(-2000, 2000, units='vel', wv_ref=WV_LINE)
spectrum_baseline_sub = spectrum_cutted.sub_baseline(wv_line=WV_LINE)
integrated_intensity = spectrum_baseline_sub.line_integ(wv_line=WV_LINE, profile='gaus', control_plot=True)
doppler_shift_line = spectrum_baseline_sub.line_velocity(wv_line=WV_LINE, control_plot=False)
"""




















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
