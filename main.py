from Cube import Cube
import numpy as np
import matplotlib.pyplot as plt

file = "/Users/delabrov/Desktop/obs_finales/NIRSpec/jw01644005001_05101_nirspec_g235h-f170lp_s3d.fits"
cube = Cube(file)
wvs = cube.wvs(units='um')
spec = cube.extract_spec(radius = 6, position = [28,28], units='erg s-1 cm-2 micron-1')

plt.plot(wvs, spec, color='black')
plt.show()



