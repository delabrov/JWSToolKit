import pytest
import numpy as np
from astropy.io.fits import Header
from astropy.io import fits
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from JWSToolKit.Cube import Cube



def test_cube_initialization(mocker):
    
    # Création d'un mock pour fits.open pour éviter d'ouvrir un vrai fichier
    mock_hdul = [mocker.Mock(header={"PIXAR_SR": 1.0, "BUNIT": "MJy/sr"}),
                 mocker.Mock(header={"PIXAR_SR": 1.0, "BUNIT": "MJy/sr"}, data=np.zeros((5, 5, 5))),
                 mocker.Mock(data=np.zeros((5, 5, 5)))]

    mocker.patch("astropy.io.fits.open", return_value=mock_hdul)

    cube = Cube("fake_file.fits")

    assert cube.file_name == "fake_file.fits"
    assert cube.data.shape == (5, 5, 5)
    assert cube.units == "MJy/sr"
    assert cube.px_area == 1.0


def test_cube_invalid_filename():
    with pytest.raises(TypeError):
        Cube(123)  


def test_get_wvs(mocker):

    mock_hdul = [
        mocker.Mock(header={}),  # primary_hdu
        mocker.Mock(header={"CRPIX3": 1, "CRVAL3": 2.0, "CDELT3": 0.1, "NAXIS3": 5, 'PIXAR_SR': 1.0, 'BUNIT': 'MJy/sr'}, data=np.zeros((5, 5, 5))),  # sci_hdu
        mocker.Mock(data=np.zeros((5, 5, 5))) 
    ]

    mocker.patch("astropy.io.fits.open", return_value=mock_hdul)

    cube = Cube("fake_file.fits")  

    wvs = cube.get_wvs(units="um")

    assert len(wvs) == 5
    assert np.isclose(wvs[0], 2.0)
    assert np.isclose(wvs[1], 2.1)

    with pytest.raises(Exception):
        cube.get_wvs(units="invalid_unit")


def test_extract_spec_circ_aperture(mocker):

    fake_data = np.random.rand(5, 5, 5)  

    mock_hdul = [
        mocker.Mock(header={}),  # primary_hdu
        mocker.Mock(header={"PIXAR_SR": 1.0, "BUNIT": "MJy/sr"}), 
        mocker.Mock(data=fake_data),
        mocker.Mock(errs=fake_data)
    ]

    mocker.patch("astropy.io.fits.open", return_value=mock_hdul)

    cube = Cube("fake_file.fits")
    cube.data = fake_data
    cube.errs = fake_data
    cube.size = np.shape(fake_data)


    spectrum = cube.extract_spec_circ_aperture(radius=1, position=[2, 2])
    assert len(spectrum) == 5  

    with pytest.raises(TypeError):
        cube.extract_spec_circ_aperture(radius="invalid", position=[2, 2])  

    with pytest.raises(Exception):
        cube.extract_spec_circ_aperture(radius=1, position=[2])  

"""
def test_cube_rotation(mocker):

    fake_header = Header({
    "PIXAR_SR": 1.0,
    "BUNIT": "MJy/sr",
    "CDELT1": 0.1,
    "CDELT2": 0.1,
    "CUNIT1": "deg",
    "CUNIT2": "deg",
    "CRPIX3": 1, 
    "CRVAL3": 2.0, 
    "CDELT3": 0.1, 
    "NAXIS3": 5,
    "PC1_1": -1,
    "PC1_2": 0,
    "PC1_3": 0,
    "PC2_1": 0,
    "PC2_2": 1,
    "PC2_3": 0,
    "PC3_1": 0,
    "PC3_2": 0,
    "PC3_3": 1
    })

    fake_data = np.random.rand(5, 5, 5)
    
    mock_hdul = [mocker.Mock(header=fake_header),
                 mocker.Mock(data=fake_data),
                 mocker.Mock(errs=fake_data)
                 ]

    mocker.patch("astropy.io.fits.open", return_value=mock_hdul)

    cube = Cube("fake_file.fits")
    cube.data_header = fake_header.copy()
    cube.data = fake_data
    cube.errs = fake_data
    cube.size = np.shape(fake_data)
    cube.px_area = 1.0
    rotated_cube = cube.rotate(angle=45)

    assert rotated_cube.size == cube.size  
    assert isinstance(rotated_cube, Cube)  
"""






def test_cube_info(mocker, capsys):
    
    fake_header = {'PI_NAME': 'test', 
                   'PIXAR_SR': 1.0, 
                   'BUNIT': 'MJy/sr',
                   'TITLE': 'test',
                   'PROGRAM': 'test',
                   'TARGNAME': 'test',
                   'TELESCOP': 'test',
                   'INSTRUME': 'test',
                   'GRATING': 'test',
                   'FILTER': 'test',
                   'NINTS': 'test',
                   'NGROUPS': 'test',
                   'NFRAMES': 'test',
                   'NUMDTHPT': 2,
                   'PATTTYPE': 'test',
                   'DATE-OBS': 'test',
                   'TIME-OBS': 'test',
                   'TARG_RA': 4,
                   'TARG_DEC': 21,
                   'EFFEXPTM': 31,
                   'DURATION': 3000,
                   }
    fake_data = np.random.rand(5, 5, 5)

    mock_hdul = [
        mocker.Mock(header={}),  # primary_hdu
        mocker.Mock(header=fake_header, data=fake_data),  # sci_hdu
        mocker.Mock(data=fake_data) 
    ]

    mocker.patch("astropy.io.fits.open", return_value=mock_hdul)

    cube = Cube("fake_file.fits")  
    print(type(cube.data_header))
    cube.data_header = fake_header
    cube.info()

    captured = capsys.readouterr()
    assert "DATA CUBE INFORMATION" in captured.out  # Vérifie que du texte est bien affiché


"""
def test_cube_full_pipeline(mocker):
    
    mock_hdul = [mocker.Mock(header={"PIXAR_SR": 1.0, "BUNIT": "MJy/sr", "CRPIX3": 1, "CRVAL3": 2.0, "CDELT3": 0.1, "NAXIS3": 5}),
                 mocker.Mock(data=np.random.rand(5, 5, 5))]

    mocker.patch("astropy.io.fits.open", return_value=mock_hdul)

    cube = Cube("fake_file.fits")

    # Extraire un spectre
    spectrum = cube.extract_spec_circ_aperture(radius=1, position=[2, 2])
    assert len(spectrum) == 5  # Vérifie la taille du spectre

    # Rotation
    rotated_cube = cube.rotate(angle=45)
    assert isinstance(rotated_cube, Cube)  # Vérifie que la rotation fonctionne

    # Vérifier les longueurs d'onde
    wvs = rotated_cube.get_wvs()
    assert len(wvs) == 5  # Vérifie la taille

"""
