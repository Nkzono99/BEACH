import numpy as np

from beach import (
    BEMElement,
    BEMMesh,
    RandomBeamInjector,
    ShiftedMaxwellVelocitySampler,
    SimConfig,
    TestParticleSimulator,
    UniformPositionSampler,
    ZeroB,
)


def main() -> None:
    tri = BEMElement(
        v0=np.array([0.0, 0.0, 0.0]),
        v1=np.array([1.0, 0.0, 0.0]),
        v2=np.array([0.0, 1.0, 0.0]),
        q=0.0,
    )
    mesh = BEMMesh([tri])

    qe = -1.602176634e-19
    me = 9.10938356e-31

    rng = np.random.default_rng(123)
    pos_sampler = UniformPositionSampler(
        low=np.array([0.1, 0.1, 0.45]),
        high=np.array([0.4, 0.4, 0.55]),
        rng=rng,
    )
    vel_sampler = ShiftedMaxwellVelocitySampler(
        drift_velocity=np.array([0.0, 0.0, -2.0e5]),
        temperature_K=2.0e4,
        rng=rng,
    )
    injector = RandomBeamInjector(
        q=qe,
        m=me,
        w=1.0,
        position_sampler=pos_sampler,
        velocity_sampler=vel_sampler,
    )
    particles = injector.sample(500)

    cfg = SimConfig(dt=1e-9, npcls_per_step=100, max_step=5000, tol_rel=1e-6)
    sim = TestParticleSimulator([mesh], cfg, B_model=ZeroB())
    stats = sim.run(particles)

    print("stats:", stats)
    print("element charges [C]:", mesh.charges())


if __name__ == "__main__":
    main()
