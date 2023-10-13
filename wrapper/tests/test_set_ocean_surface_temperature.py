import unittest
import os
import sys
import numpy as np
import shield.wrapper
from copy import deepcopy
from mpi4py import MPI
from util import (
    get_default_config,
    get_state_single_variable,
    main,
)


test_dir = os.path.dirname(os.path.abspath(__file__))


def select_ocean_values(*fields):
    is_ocean = np.isclose(get_state_single_variable("land_sea_mask"), 0.0)
    return (field[is_ocean] for field in fields)


def replace_ocean_state_with_random_values(names):
    is_ocean = np.isclose(get_state_single_variable("land_sea_mask"), 0.0)
    old_state = shield.wrapper.get_state(names=names)
    replace_state = deepcopy(old_state)
    for name, quantity in replace_state.items():
        values = np.random.uniform(size=quantity.extent)
        values = np.where(is_ocean, values, quantity.view[:])
        quantity.view[:] = values
    shield.wrapper.set_state(replace_state)
    return replace_state


class PrescribeSSTTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(PrescribeSSTTests, self).__init__(*args, **kwargs)

    def setUp(self):
        pass

    def tearDown(self):
        MPI.COMM_WORLD.barrier()

    def test_prescribing_sst_changes_model_state(self):
        checkpoint_state = shield.wrapper.get_state(shield.wrapper.get_restart_names())

        # If we do not set the sea surface temperature and
        # use_climatological_sst is set to .false., the sea surface temperature
        # will remain at what it was set to in the initial conditions for the
        # duration of the run.
        shield.wrapper.step()
        air_temperature_from_default_ocean_temperature = get_state_single_variable(
            "air_temperature"
        )

        shield.wrapper.set_state(checkpoint_state)
        replace_ocean_state_with_random_values(["surface_temperature"])
        shield.wrapper.step()
        air_temperature_from_prescribed_ocean_temperature = get_state_single_variable(
            "air_temperature"
        )

        assert not np.allclose(
            air_temperature_from_default_ocean_temperature,
            air_temperature_from_prescribed_ocean_temperature,
        )

    def test_prescribing_sst_changes_surface_temperature_diagnostic(self):
        replaced_state = replace_ocean_state_with_random_values(["surface_temperature"])
        prescribed_sst = replaced_state["surface_temperature"].view[:]
        shield.wrapper.step()
        surface_temperature_diagnostic = shield.wrapper.get_diagnostic_by_name(
            "tsfc", module_name="gfs_sfc"
        ).view[:]

        result, expected = select_ocean_values(
            surface_temperature_diagnostic, prescribed_sst
        )
        np.testing.assert_allclose(result, expected)


if __name__ == "__main__":
    config = get_default_config()

    # Clear diag_table for these tests, since outputting interval-averaged
    # physics diagnostics leads to unrelated errors; see
    # ai2cm/fv3gfs-fortran#384 for more context.
    config["diag_table"] = "no_output"

    # Turn off the mixed layer ocean model for now.  We can think about how to
    # prescribe the reference SST for the mixed layer ocean model later.
    config["namelist"]["gfs_physics_nml"]["do_ocean"] = False
    main(test_dir, config)
