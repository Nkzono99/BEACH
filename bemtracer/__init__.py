"""bem_charge_tracer package."""

from .sim import (
    BEMElement,
    BEMMesh,
    Particle,
    SimConfig,
    TestParticleSimulator,
    calc_electric_field,
    BEMField,
    ZeroB,
    UniformB,
    AbsorptionInteraction,
    InsulatorSurfaceCharge,
)
from .geometry import (
    make_plane,
    make_box,
    make_cylinder,
    make_sphere,
    mesh_from_vertices_faces,
    merge_meshes,
)
from .importers import MeshImporter, OBJImporter, FreeCADImporter, load_mesh
from .fortran_results import FortranRunResult, load_fortran_result, list_fortran_runs, plot_charge_mesh, plot_charges

from .injection import (
    PositionSampler,
    VelocitySampler,
    UniformPositionSampler,
    ShiftedMaxwellVelocitySampler,
    RandomBeamInjector,
    FixedBeamInjector,
)

__all__ = [
    "BEMElement",
    "BEMMesh",
    "Particle",
    "SimConfig",
    "TestParticleSimulator",
    "calc_electric_field",
    "BEMField",
    "ZeroB",
    "UniformB",
    "AbsorptionInteraction",
    "InsulatorSurfaceCharge",
    "make_plane",
    "make_box",
    "make_cylinder",
    "make_sphere",
    "mesh_from_vertices_faces",
    "merge_meshes",
    "MeshImporter",
    "OBJImporter",
    "FreeCADImporter",
    "load_mesh",
    "PositionSampler",
    "VelocitySampler",
    "UniformPositionSampler",
    "ShiftedMaxwellVelocitySampler",
    "RandomBeamInjector",
    "FixedBeamInjector",
    "FortranRunResult",
    "load_fortran_result",
    "list_fortran_runs",
    "plot_charges",
    "plot_charge_mesh",
]
