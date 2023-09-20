import unittest
import os
import numpy as np
import shield.wrapper
from shield.wrapper._properties import FLAGSTRUCT_PROPERTIES
from mpi4py import MPI

from util import get_default_config, generate_data_dict, main

test_dir = os.path.dirname(os.path.abspath(__file__))
FORTRAN_TO_PYTHON_TYPE = {"integer": int, "real": float, "logical": bool}


class FlagsTest(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(FlagsTest, self).__init__(*args, **kwargs)
        self.flagstruct_data = generate_data_dict(FLAGSTRUCT_PROPERTIES)
        self.mpi_comm = MPI.COMM_WORLD

    def setUp(self):
        pass

    def tearDown(self):
        self.mpi_comm.barrier()

    def test_flagstruct_properties_present_in_metadata(self):
        """Test that some small subset of flagstruct names are in the data dictionary"""
        for name in ["do_adiabatic_init", "ptop"]:
            self.assertIn(name, self.flagstruct_data.keys())

    def test_get_all_flagstruct_properties(self):
        self._get_all_properties_helper(self.flagstruct_data)

    def _get_all_properties_helper(self, properties):
        for name, data in properties.items():
            with self.subTest(name):
                result = getattr(shield.wrapper.flags, name)
                expected_type = FORTRAN_TO_PYTHON_TYPE[data["type_fortran"]]
                self.assertIsInstance(result, expected_type)

    def test_ptop(self):
        """Test that getting a real flag produces its expected result."""
        result = shield.wrapper.flags.ptop
        expected = 64.247
        np.testing.assert_allclose(result, expected)

    def test_n_split(self):
        """Test that getting an integer flag produces its expected result."""
        result = shield.wrapper.flags.n_split
        expected = 6
        self.assertEqual(result, expected)

    def test_dt_atmos(self):
        result = shield.wrapper.flags.dt_atmos
        expected = 900
        self.assertEqual(result, expected)


if __name__ == "__main__":
    config = get_default_config()
    main(test_dir, config)
