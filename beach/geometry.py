from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np

from .sim import BEMElement, BEMMesh


ArrayLike3 = Sequence[float] | np.ndarray


def mesh_from_vertices_faces(vertices: np.ndarray, faces: np.ndarray) -> BEMMesh:
    """Build a ``BEMMesh`` from vertex/triangle-index arrays."""
    vv = np.asarray(vertices, dtype=float)
    ff = np.asarray(faces, dtype=int)

    if vv.ndim != 2 or vv.shape[1] != 3:
        raise ValueError(f"vertices must have shape (N, 3), got {vv.shape}")
    if ff.ndim != 2 or ff.shape[1] != 3:
        raise ValueError(f"faces must have shape (M, 3), got {ff.shape}")
    if ff.size and (ff.min() < 0 or ff.max() >= vv.shape[0]):
        raise ValueError("faces contain out-of-range vertex index")

    elements = [
        BEMElement(v0=vv[i0].copy(), v1=vv[i1].copy(), v2=vv[i2].copy(), q=0.0)
        for i0, i1, i2 in ff
    ]
    return BEMMesh(elements)


def make_plane(
    *,
    size_x: float = 1.0,
    size_y: float = 1.0,
    nx: int = 1,
    ny: int = 1,
    center: ArrayLike3 = (0.0, 0.0, 0.0),
) -> BEMMesh:
    """Create a rectangular plane in the XY plane, split into triangles."""
    if nx <= 0 or ny <= 0:
        raise ValueError("nx and ny must be positive")

    cx, cy, cz = np.asarray(center, dtype=float)
    x0, x1 = cx - 0.5 * size_x, cx + 0.5 * size_x
    y0, y1 = cy - 0.5 * size_y, cy + 0.5 * size_y

    xs = np.linspace(x0, x1, nx + 1)
    ys = np.linspace(y0, y1, ny + 1)
    elements: list[BEMElement] = []

    for ix in range(nx):
        for iy in range(ny):
            p00 = np.array([xs[ix], ys[iy], cz], dtype=float)
            p10 = np.array([xs[ix + 1], ys[iy], cz], dtype=float)
            p01 = np.array([xs[ix], ys[iy + 1], cz], dtype=float)
            p11 = np.array([xs[ix + 1], ys[iy + 1], cz], dtype=float)
            elements.append(BEMElement(v0=p00, v1=p10, v2=p11))
            elements.append(BEMElement(v0=p00, v1=p11, v2=p01))

    return BEMMesh(elements)


def make_box(
    *,
    size: ArrayLike3 = (1.0, 1.0, 1.0),
    center: ArrayLike3 = (0.0, 0.0, 0.0),
    nx: int = 1,
    ny: int = 1,
    nz: int = 1,
) -> BEMMesh:
    """Create a closed box surface mesh as triangles."""
    sx, sy, sz = np.asarray(size, dtype=float)
    cx, cy, cz = np.asarray(center, dtype=float)
    hx, hy, hz = 0.5 * sx, 0.5 * sy, 0.5 * sz

    faces: list[BEMElement] = []

    # z = +/- hz planes
    for sign in (-1.0, 1.0):
        z = cz + sign * hz
        plane = make_plane(size_x=sx, size_y=sy, nx=nx, ny=ny, center=(cx, cy, z))
        for e in plane.elements:
            if sign < 0:
                faces.append(BEMElement(v0=e.v0, v1=e.v2, v2=e.v1))
            else:
                faces.append(e)

    # x = +/- hx planes
    ys = np.linspace(cy - hy, cy + hy, ny + 1)
    zs = np.linspace(cz - hz, cz + hz, nz + 1)
    for sign in (-1.0, 1.0):
        x = cx + sign * hx
        for iy in range(ny):
            for iz in range(nz):
                p00 = np.array([x, ys[iy], zs[iz]], dtype=float)
                p10 = np.array([x, ys[iy + 1], zs[iz]], dtype=float)
                p01 = np.array([x, ys[iy], zs[iz + 1]], dtype=float)
                p11 = np.array([x, ys[iy + 1], zs[iz + 1]], dtype=float)
                if sign < 0:
                    faces.append(BEMElement(v0=p00, v1=p11, v2=p10))
                    faces.append(BEMElement(v0=p00, v1=p01, v2=p11))
                else:
                    faces.append(BEMElement(v0=p00, v1=p10, v2=p11))
                    faces.append(BEMElement(v0=p00, v1=p11, v2=p01))

    # y = +/- hy planes
    xs = np.linspace(cx - hx, cx + hx, nx + 1)
    zs = np.linspace(cz - hz, cz + hz, nz + 1)
    for sign in (-1.0, 1.0):
        y = cy + sign * hy
        for ix in range(nx):
            for iz in range(nz):
                p00 = np.array([xs[ix], y, zs[iz]], dtype=float)
                p10 = np.array([xs[ix + 1], y, zs[iz]], dtype=float)
                p01 = np.array([xs[ix], y, zs[iz + 1]], dtype=float)
                p11 = np.array([xs[ix + 1], y, zs[iz + 1]], dtype=float)
                if sign < 0:
                    faces.append(BEMElement(v0=p00, v1=p10, v2=p11))
                    faces.append(BEMElement(v0=p00, v1=p11, v2=p01))
                else:
                    faces.append(BEMElement(v0=p00, v1=p11, v2=p10))
                    faces.append(BEMElement(v0=p00, v1=p01, v2=p11))

    return BEMMesh(faces)


def make_cylinder(
    *,
    radius: float = 0.5,
    height: float = 1.0,
    n_theta: int = 24,
    n_z: int = 1,
    cap: bool = True,
    center: ArrayLike3 = (0.0, 0.0, 0.0),
) -> BEMMesh:
    """Create a cylinder aligned with Z axis."""
    if n_theta < 3 or n_z <= 0:
        raise ValueError("n_theta >= 3 and n_z > 0 are required")

    cx, cy, cz = np.asarray(center, dtype=float)
    zvals = np.linspace(cz - 0.5 * height, cz + 0.5 * height, n_z + 1)
    thetas = np.linspace(0.0, 2 * np.pi, n_theta, endpoint=False)

    elems: list[BEMElement] = []
    for iz in range(n_z):
        z0, z1 = zvals[iz], zvals[iz + 1]
        for it in range(n_theta):
            t0 = thetas[it]
            t1 = thetas[(it + 1) % n_theta]
            p00 = np.array([cx + radius * np.cos(t0), cy + radius * np.sin(t0), z0])
            p10 = np.array([cx + radius * np.cos(t1), cy + radius * np.sin(t1), z0])
            p01 = np.array([cx + radius * np.cos(t0), cy + radius * np.sin(t0), z1])
            p11 = np.array([cx + radius * np.cos(t1), cy + radius * np.sin(t1), z1])
            elems.append(BEMElement(v0=p00, v1=p10, v2=p11))
            elems.append(BEMElement(v0=p00, v1=p11, v2=p01))

    if cap:
        z_bottom = zvals[0]
        z_top = zvals[-1]
        c_bottom = np.array([cx, cy, z_bottom])
        c_top = np.array([cx, cy, z_top])
        for it in range(n_theta):
            t0 = thetas[it]
            t1 = thetas[(it + 1) % n_theta]
            b0 = np.array(
                [cx + radius * np.cos(t0), cy + radius * np.sin(t0), z_bottom]
            )
            b1 = np.array(
                [cx + radius * np.cos(t1), cy + radius * np.sin(t1), z_bottom]
            )
            t0p = np.array([cx + radius * np.cos(t0), cy + radius * np.sin(t0), z_top])
            t1p = np.array([cx + radius * np.cos(t1), cy + radius * np.sin(t1), z_top])
            elems.append(BEMElement(v0=c_bottom, v1=b1, v2=b0))
            elems.append(BEMElement(v0=c_top, v1=t0p, v2=t1p))

    return BEMMesh(elems)


def make_sphere(
    *,
    radius: float = 0.5,
    n_lon: int = 24,
    n_lat: int = 12,
    center: ArrayLike3 = (0.0, 0.0, 0.0),
) -> BEMMesh:
    """Create a UV sphere triangulation."""
    if n_lon < 3 or n_lat < 2:
        raise ValueError("n_lon >= 3 and n_lat >= 2 are required")

    c = np.asarray(center, dtype=float)
    lon = np.linspace(0.0, 2 * np.pi, n_lon, endpoint=False)
    lat = np.linspace(0.0, np.pi, n_lat + 1)

    def sph(theta: float, phi: float) -> np.ndarray:
        x = radius * np.sin(phi) * np.cos(theta)
        y = radius * np.sin(phi) * np.sin(theta)
        z = radius * np.cos(phi)
        return c + np.array([x, y, z], dtype=float)

    elems: list[BEMElement] = []
    for ilat in range(n_lat):
        p0 = lat[ilat]
        p1 = lat[ilat + 1]
        for ilon in range(n_lon):
            t0 = lon[ilon]
            t1 = lon[(ilon + 1) % n_lon]
            a = sph(t0, p0)
            b = sph(t1, p0)
            c0 = sph(t0, p1)
            d = sph(t1, p1)

            if ilat == 0:
                elems.append(BEMElement(v0=a, v1=c0, v2=d))
            elif ilat == n_lat - 1:
                elems.append(BEMElement(v0=a, v1=d, v2=b))
            else:
                elems.append(BEMElement(v0=a, v1=d, v2=b))
                elems.append(BEMElement(v0=a, v1=c0, v2=d))

    return BEMMesh(elems)


def merge_meshes(meshes: Iterable[BEMMesh]) -> BEMMesh:
    """Merge multiple ``BEMMesh`` instances into one."""
    elems: list[BEMElement] = []
    for m in meshes:
        elems.extend(
            BEMElement(v0=e.v0.copy(), v1=e.v1.copy(), v2=e.v2.copy(), q=e.q)
            for e in m.elements
        )
    return BEMMesh(elems)
