from Cube import Cube
from Spec import Spec
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

file = "/Users/delabrov/Desktop/obs_finales/NIRSpec/jw01644005001_05101_nirspec_g235h-f170lp_s3d.fits"
cube = (Cube(file))

wvs = cube.wvs()
spec = cube.extract_spec(radius = 6, position = (28,28), units = 'Jy')

spec = Spec(wvs, spec, units = 'Jy')
spec_red = spec.cut(-150, 150, units = 'vel', wv_ref = 2.12)

plt.plot(spec_red.wvs, spec_red.values, color='black')
plt.show()




