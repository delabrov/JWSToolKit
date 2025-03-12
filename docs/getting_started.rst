Getting Started
===============

The purpose of this section is to help users get familiar with the JWSToolKit tools. 

Add JWSToolKit to your Python code 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The package is divided into several modules corresponding to the different Classes. The example below shows how to import the 'Cube' and 'Spec' classes into Python code. 

.. code-block:: python

    from Cube import Cube 
    from Spec import Spec

Open a JWST 'data cube' file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

JWST observations with a cube structure (from the NIRSpec-IFU or MIRI-MRS instruments, for example) have a '.s3d' suffix in their filename, 
reflecting the fact that the observations are in the form of a 3D array. By default, pipeline output files are in FITS format, 
which can be opened from the *astropy.io* package. 

.. code-block:: python 

    file = "DataCube_s3d.fits"
    cube = Cube(file)

    cube.info()


By using the .info() method on the *Cube* object, a description of the data information is printed in the console. 
This method reads the two main headers of JWST observations: the 'primary header' and the 'science header'. 

.. code-block:: console

        __________ DATA CUBE INFORMATION __________
    Data file name:DataCube_s3d.fits
    Program PI: Dougados, Catherine, for the project: A cornerstone study of the jet/outflow connexion: the remarkable DG Tau B system
    Program ID: 01644
    Target: DG TAU B
    Telescope: JWST \ Instrument: NIRSPEC
    Configuration: G140H + F100LP
    Number of integrations, groups and frames: 2, 20, 1
    Dither strategy: True
    Dither patern type: 4-POINT-DITHER

    Date and time of observations: 2022-09-05 | 14:04:52.091999
    Target position in the sky: RA(J2000) = 66.76086125 , Dec(J2000) = 26.09164722222222
    Effecive Exposure Time: 9336.896 s
    Total Exposure Time (with overheads): 9803.728 s

    Data type and shape: Data Cube | 77, 67, 3915 (x, y, wvs)
    Spatial pixel sizes in deg (dx, dy): 2.77777781916989e-05, 2.77777781916989e-05
    Spatial pixel sizes in arcsec (dx, dy): 0.1, 0.1
    Spectral pixel size (µm): 0.000235
    Spectral range of observations (µm): 0.9700000000000001 - 1.89
    Units of spectral pixel values: MJy/sr


In this example, we see information such as the name of the project's principal investigator, the name of the source, the instrument and the instrument configuration. 
The second section gives the temporal information of the data: date and time of observations, integration time and target position. 
The last section gives physical information on the data, such as the unit of values stored in the cube, cube size, increment values, etc. 

For example, here the cube has 3915 points on the wavelength axis, for 77 pixels on the horizontal axis and 67 pixels on the vertical axis. 


Example of how to use the methods 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example, we will extract a spectrum from the data cube, inside a circular aperture, and create a *Spec* object. 
The extract_spec_circ_aperture(*params*) method is used to retrieve a 1D vector containing the summed intensity or flux values in a circular aperture: 

.. code-block:: python

    radius      = 5         # unit: px
    position    = [25, 29]  # unit : px

    # This function extracts the integrated spectrum
    spectrum_values = cube.extract_spec_circ_aperture(radius, position, units='Jy')

    # This function retrieves the wavelength axis
    wvs_values = cube.get_wvs()

    # Building the Spec object
    spectrum = Spec(wvs_values, spectrum_values, units='Jy')


