import numpy as np

from beach import (
    BEMElement,
    BEMMesh,
    OBJImporter,
    FreeCADImporter,
    make_box,
    make_cylinder,
    make_plane,
    make_sphere,
    merge_meshes,
)


def test_plane_triangle_count_and_area() -> None:
    mesh = make_plane(size_x=2.0, size_y=3.0, nx=5, ny=4)
    assert mesh.nelem == 2 * 5 * 4
    assert np.isclose(mesh.areas.sum(), 6.0)


def test_box_triangle_count() -> None:
    nx, ny, nz = 2, 3, 4
    mesh = make_box(size=(1.0, 1.0, 1.0), nx=nx, ny=ny, nz=nz)
    expected = 4 * (nx * ny + ny * nz + nx * nz)
    assert mesh.nelem == expected


def test_cylinder_triangle_count() -> None:
    n_theta, n_z = 16, 5
    mesh = make_cylinder(radius=1.0, height=2.0, n_theta=n_theta, n_z=n_z, cap=True)
    expected_side = 2 * n_theta * n_z
    expected_caps = 2 * n_theta
    assert mesh.nelem == expected_side + expected_caps


def test_sphere_triangle_count() -> None:
    n_lon, n_lat = 20, 10
    mesh = make_sphere(radius=1.0, n_lon=n_lon, n_lat=n_lat)
    # top + bottom: 2*n_lon, middle strips: 2*n_lon*(n_lat-2)
    expected = 2 * n_lon * (n_lat - 1)
    assert mesh.nelem == expected


def test_sphere_normals_point_outward() -> None:
    mesh = make_sphere(radius=1.0, n_lon=20, n_lat=10, center=(0.0, 0.0, 0.0))
    assert np.all(np.einsum("ij,ij->i", mesh.normals, mesh.centers) > 0.0)


def test_merge_meshes_copies_elements() -> None:
    base = BEMElement(
        v0=np.array([0.0, 0.0, 0.0]),
        v1=np.array([1.0, 0.0, 0.0]),
        v2=np.array([0.0, 1.0, 0.0]),
    )
    mesh_a = BEMMesh([base])
    mesh_b = merge_meshes([mesh_a])
    mesh_b.elements[0].q = 5.0
    assert mesh_a.elements[0].q == 0.0


def test_obj_importer(tmp_path) -> None:
    p = tmp_path / "shape.obj"
    p.write_text(
        "\n".join(
            [
                "v 0 0 0",
                "v 1 0 0",
                "v 1 1 0",
                "v 0 1 0",
                "f 1 2 3 4",
            ]
        ),
        encoding="utf-8",
    )

    mesh = OBJImporter().load(p)
    assert mesh.nelem == 2
    assert np.isclose(mesh.areas.sum(), 1.0)


def test_freecad_importer_hook() -> None:
    with np.testing.assert_raises(NotImplementedError):
        FreeCADImporter().load("dummy.fcstd")
