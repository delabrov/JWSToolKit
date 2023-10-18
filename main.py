from Cube import Cube
from Spec import Spec
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

file = "/Users/delabrov/Desktop/obs_finales/NIRSpec/jw01644005001_05101_nirspec_g235h-f170lp_s3d.fits"
cube = Cube(file)
wvs = cube.wvs()
spec = cube.extract_spec(radius = 6, position = (28,28), units = 'erg s-1 cm-2 micron-1')
spec = Spec(wvs, spec, units = 'erg s-1 cm-2 micron-1')
spec_red = spec.cut(min = -1000, max = 1000, units = 'vel', wv_ref = 2.121833725)
spec_contSub = spec_red.sub_baseline(wv_line = 2.121833725)





plt.plot(spec_red.wvs, spec_red.values, color='black')
plt.plot(spec_contSub.wvs, spec_contSub.values, color='green')
plt.show()




