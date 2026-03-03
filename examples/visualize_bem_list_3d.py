import numpy as np

from bemtracer import (
    BEMElement,
    BEMMesh,
    Particle,
    SimConfig,
    TestParticleSimulator,
    ZeroB,
)


def build_plate_mesh(nx: int = 8, ny: int = 8, size: float = 1.0) -> BEMMesh:
    """Create a square plate at z=0, triangulated into 2*nx*ny elements."""
    elements: list[BEMElement] = []
    xs = np.linspace(-0.5 * size, 0.5 * size, nx + 1)
    ys = np.linspace(-0.5 * size, 0.5 * size, ny + 1)

    for ix in range(nx):
        for iy in range(ny):
            p00 = np.array([xs[ix], ys[iy], 0.0])
            p10 = np.array([xs[ix + 1], ys[iy], 0.0])
            p01 = np.array([xs[ix], ys[iy + 1], 0.0])
            p11 = np.array([xs[ix + 1], ys[iy + 1], 0.0])

            elements.append(BEMElement(v0=p00, v1=p10, v2=p11, q=0.0))
            elements.append(BEMElement(v0=p00, v1=p11, v2=p01, q=0.0))

    return BEMMesh(elements)


def build_beam(n: int = 2000) -> list[Particle]:
    """Generate a focused electron beam toward the plate."""
    qe = -1.602176634e-19
    me = 9.10938356e-31

    rng = np.random.default_rng(42)
    xy = rng.normal(loc=0.0, scale=0.08, size=(n, 2))

    particles: list[Particle] = []
    for i in range(n):
        particles.append(
            Particle(
                x=np.array([xy[i, 0], xy[i, 1], 0.5]),
                v=np.array([0.0, 0.0, -2.5e5]),
                q=qe,
                m=me,
                w=1.0,
            )
        )
    return particles


def plot_bem_charges_3d(bem_list: list[BEMMesh], *, title: str = "BEM element charge distribution") -> None:
    """3D plot of charge on each triangle in bem_list."""
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    faces: list[np.ndarray] = []
    charges: list[float] = []
    for mesh in bem_list:
        for elem in mesh.elements:
            faces.append(np.stack([elem.v0, elem.v1, elem.v2], axis=0))
            charges.append(elem.q)

    q = np.asarray(charges, dtype=float)
    q_abs_max = np.max(np.abs(q)) if q.size else 1.0
    if q_abs_max == 0.0:
        q_abs_max = 1.0

    norm = colors.TwoSlopeNorm(vmin=-q_abs_max, vcenter=0.0, vmax=q_abs_max)
    cmap = cm.get_cmap("coolwarm")

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")

    poly = Poly3DCollection(faces, linewidths=0.2, edgecolors="k")
    poly.set_facecolor(cmap(norm(q)))
    ax.add_collection3d(poly)

    vertices = np.concatenate(faces, axis=0)
    mins = vertices.min(axis=0)
    maxs = vertices.max(axis=0)
    center = 0.5 * (mins + maxs)
    radius = 0.55 * np.max(maxs - mins + 1e-12)
    ax.set_xlim(center[0] - radius, center[0] + radius)
    ax.set_ylim(center[1] - radius, center[1] + radius)
    ax.set_zlim(center[2] - radius, center[2] + radius)

    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")
    ax.set_title(title)

    mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array(q)
    cbar = fig.colorbar(mappable, ax=ax, pad=0.1)
    cbar.set_label("element charge q [C]")

    plt.tight_layout()
    plt.show()


def main() -> None:
    mesh = build_plate_mesh(nx=10, ny=10, size=1.0)
    bem_list = [mesh]

    particles = build_beam(n=2500)

    cfg = SimConfig(
        dt=1e-9,
        npcls_per_step=100,
        max_step=5000,
        tol_rel=1e-6,
        use_hybrid=True,
        r_switch_factor=3.0,
        n_sub=2,
        softening_factor=0.1,
    )

    sim = TestParticleSimulator(bem_list, cfg, B_model=ZeroB())
    stats = sim.run(particles)

    print("stats:", stats)
    print(
        "charge min/max [C]:",
        float(mesh.charges().min()),
        float(mesh.charges().max()),
    )
    print("total charge [C]:", float(mesh.charges().sum()))

    plot_bem_charges_3d(bem_list, title="Plate charge after absorption")


if __name__ == "__main__":
    main()
