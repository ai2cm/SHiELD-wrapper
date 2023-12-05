import pace.util
from ._wrapper import (
    _get_diagnostic_data,
    _get_diagnostic_info,
    initialize,
    step,
    step_dynamics,
    step_radiation,
    step_pre_radiation,
    step_post_radiation_physics,
    step_physics,
    save_intermediate_restart_if_enabled,
    save_fortran_restart,
    cleanup,
    get_state,
    set_state,
    get_n_ghost_cells,
    get_step_count,
    get_tracer_metadata,
    compute_physics,
    apply_physics,
    flags,
    DiagnosticInfo,
    transform_agrid_winds_to_dgrid_winds,
    transform_dgrid_winds_to_agrid_winds
)
from ._restart import get_restart_names, open_restart
from . import examples

from .thermodynamics import set_state_mass_conserving


def get_diagnostic_by_name(
    name: str, module_name: str = "gfs_phys"
) -> pace.util.Quantity:
    """Get a diagnostic field as a Quantity

    Currently, only supports diagnostics defined in the FV3GFS_io.F90
    """
    metadata = get_diagnostic_metadata_by_name(name, module_name)
    return _get_diagnostic_data(metadata.diag_manager_controlled, metadata.index)


def get_diagnostic_metadata_by_name(
    name: str, module_name: str = "gfs_phys"
) -> DiagnosticInfo:
    """Get diagnostic metadata by name

    Currently, only supports diagnostics defined in the FV3GFS_io.F90
    """
    info = _get_diagnostic_info()

    key = module_name, name
    if key in info:
        return info[key]
    else:
        raise ValueError(f"There is no diagnostic {name} in module {module_name}.")


__version__ = "0.1.0"

__all__ = list(key for key in locals().keys() if not key.startswith("_"))
