import hashlib
import os
import subprocess
import typing

import fv3config
import pytest

from pathlib import Path


TEST_DIR = Path(__file__).parent
CONFIG_DIR = TEST_DIR / "config"
EXECUTABLE = Path("/SHiELD/SHiELD_build/Build/bin/SHiELD_nh.prod.64bit.gnu.x")


def get_config(filename):
    config_filename = CONFIG_DIR / filename
    with open(config_filename, "r") as f:
        return fv3config.load(f)


def checksum_file(path: Path) -> str:
    sum = hashlib.md5()
    BUFFER_SIZE = 1024 * 1024
    with open(path, "rb") as f:
        while True:
            buf = f.read(BUFFER_SIZE)
            if not buf:
                break
            sum.update(buf)
    return sum.hexdigest()


def get_rundir_netcdfs(rundir: Path) -> typing.List[Path]:
    restart_directory = rundir / "RESTART"
    diagnostics_files = sorted(rundir.glob("*.nc"))
    restart_files = sorted(restart_directory.glob("*.nc"))
    return diagnostics_files + restart_files    


def checksum_rundir_to_dict(rundir: Path) -> typing.Dict[str, str]:
    rundir_netcdfs = get_rundir_netcdfs(rundir)
    return {file.name: checksum_file(file) for file in rundir_netcdfs}
        

def checksum_rundir_to_file(rundir: Path, file: Path):
    """checksum rundir storing output in file"""
    rundir_netcdfs = get_rundir_netcdfs(rundir)
    for path in rundir_netcdfs:
        print(path, checksum_file(path), file=file)


def _run_simulation(config: dict, rundir: Path, command: typing.List[str]):
    fv3config.write_run_directory(config, rundir)
    n_processes = fv3config.config.get_n_processes(config)
    command = ["mpirun", "-n", str(n_processes)] + command
    completed_process = subprocess.run(command, cwd=rundir, capture_output=True)
    if completed_process.returncode != 0:
        print("Tail of Stderr:")
        print(completed_process.stderr[-2000:].decode())
        print("Tail of Stdout:")
        print(completed_process.stdout[-2000:].decode())
        pytest.fail()
    return completed_process
    
        
def run_fortran_executable(config: dict, rundir: Path):
    command = [EXECUTABLE.absolute().as_posix()]
    return _run_simulation(config, rundir, command)


def run_python_wrapper(config: dict, rundir: Path):
    command = ["python3", "-m", "mpi4py", "-m", "shield.wrapper.run"]
    return _run_simulation(config, rundir, command)
    

@pytest.mark.parametrize("config_filename", ["default.yml"])
def test_regression(config_filename: str, tmp_path: Path, regtest):
    config = get_config(config_filename)
    rundir_fortran = tmp_path / "rundir-fortran"
    rundir_wrapper = tmp_path / "rundir-wrapper"

    run_fortran_executable(config, rundir_fortran)
    run_python_wrapper(config, rundir_wrapper)

    checksums_fortran = checksum_rundir_to_dict(rundir_fortran)
    checksums_wrapper = checksum_rundir_to_dict(rundir_wrapper) 
    
    assert checksums_wrapper == checksums_fortran
    checksum_rundir_to_file(rundir_wrapper, file=regtest)
