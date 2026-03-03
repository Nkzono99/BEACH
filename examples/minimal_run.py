import numpy as np
from beach import (
    BEMElement,
    BEMMesh,
    SimConfig,
    TestParticleSimulator,
    ZeroB,
    FixedBeamInjector,
)


def main() -> None:
    # Simple single triangle (z=0 plane)
    tri = BEMElement(
        v0=np.array([0.0, 0.0, 0.0]),
        v1=np.array([1.0, 0.0, 0.0]),
        v2=np.array([0.0, 1.0, 0.0]),
        q=0.0,
    )
    mesh = BEMMesh([tri])
    bem_list = [mesh]

    # Electron beam toward the triangle
    qe = -1.602176634e-19
    me = 9.10938356e-31

    injector = FixedBeamInjector(
        x0=np.array([0.25, 0.25, 0.5]),
        v0=np.array([0.0, 0.0, -2e5]),
        q=qe,
        m=me,
        w=1.0,
    )
    particles = injector.sample(500)

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
    print("element charges [C]:", mesh.charges())


if __name__ == "__main__":
    main()
