"""bem_charge_tracer package."""

from .core import (
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
]
