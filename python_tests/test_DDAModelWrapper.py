import dda_model
import numpy as np
from scipy.special import sici
import pytest


def _construct_test_model():
    sym_axis = [10.5, 10.5]
    geometry_shape = [22, 22, 10]
    pixel_size = 15.0
    light_direction = [0, 0, 1]
    light_polarization = [1, 0, 0]
    wavelength = 542
    initialization = np.loadtxt("data/dielectrics.txt")
    dielectric_constants = [1.01 + 0j, 5.96282 + 3.80423e-7j]

    model = dda_model.DDAModelWrapper(
        geometry_shape,
        pixel_size,
        initialization,
        dielectric_constants=dielectric_constants,
        light_direction=light_direction,
        light_polarization=light_polarization,
        light_wavelength_nm=wavelength,
        symmetry_axes=sym_axis,
    )
    return model


def test_wrapper_objective():
    # This tests at the lower-level C++ API (not the wrapped API)
    model = _construct_test_model()
    objective_value = model.objective()
