import unittest
import os
import shield.wrapper
import pace.util
from mpi4py import MPI

from util import get_default_config, main

test_dir = os.path.dirname(os.path.abspath(__file__))


class DiagnosticTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.mpi_comm = MPI.COMM_WORLD

    def setUp(self):
        pass

    def tearDown(self):
        self.mpi_comm.barrier()

    def test_get_diag_info(self):
        output = shield.wrapper._get_diagnostic_info()
        assert len(output) > 0
        for index, item in output.items():
            self.assertIsInstance(item.axes, int)
            self.assertIsInstance(item.module_name, str)

            self.assertIsInstance(item.name, str)
            self.assertIsNot(item.name, "")

            self.assertIsInstance(item.description, str)
            self.assertIsInstance(item.unit, str)

    def test_get_diagnostic_data(self):
        names_to_get = ["eta_shal", "u10m"]
        for name in names_to_get:
            quantity = shield.wrapper.get_diagnostic_by_name(
                name, module_name="gfs_phys"
            )
            info = shield.wrapper.get_diagnostic_metadata_by_name(
                name, module_name="gfs_phys"
            )
            self.assertIsInstance(quantity, pace.util.Quantity)
            assert quantity.view[:].ndim == info.axes
            assert quantity.units == info.unit
            if name == "eta_shal":
                assert quantity.view[:].ndim == 3
            elif name == "u10m":
                assert quantity.view[:].ndim == 2
            else:
                raise ValueError(f"Testing only implemented for eta_shal and u10m")


if __name__ == "__main__":
    config = get_default_config()
    main(test_dir, config)
