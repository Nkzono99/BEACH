import numpy as np
import pytest

from bemtracer import (
    FixedBeamInjector,
    RandomBeamInjector,
    ShiftedMaxwellVelocitySampler,
    UniformPositionSampler,
)


def test_uniform_position_sampler_bounds() -> None:
    rng = np.random.default_rng(0)
    low = np.array([-1.0, 2.0, -3.0])
    high = np.array([1.0, 4.0, -1.0])
    sampler = UniformPositionSampler(low=low, high=high, rng=rng)

    x = sampler.sample(1000)
    assert x.shape == (1000, 3)
    assert np.all(x >= low[None, :])
    assert np.all(x <= high[None, :])


def test_shifted_maxwell_mean_close_to_drift() -> None:
    rng = np.random.default_rng(1)
    drift = np.array([1.0e4, -2.0e4, 3.0e4])
    sampler = ShiftedMaxwellVelocitySampler(
        drift_velocity=drift,
        thermal_speed=5.0e3,
        rng=rng,
    )

    v = sampler.sample(100000, m=9.10938356e-31)
    mean_v = v.mean(axis=0)
    assert np.all(np.abs(mean_v - drift) < 120.0)


def test_random_beam_injector_seed_reproducible() -> None:
    qe = -1.602176634e-19
    me = 9.10938356e-31
    drift = np.array([0.0, 0.0, -2.0e5])
    low = np.array([0.0, 0.0, 0.4])
    high = np.array([0.2, 0.2, 0.6])

    rng1 = np.random.default_rng(42)
    rng2 = np.random.default_rng(42)

    injector1 = RandomBeamInjector(
        q=qe,
        m=me,
        w=1.0,
        position_sampler=UniformPositionSampler(low=low, high=high, rng=rng1),
        velocity_sampler=ShiftedMaxwellVelocitySampler(
            drift_velocity=drift,
            thermal_speed=1.0e3,
            rng=rng1,
        ),
    )
    injector2 = RandomBeamInjector(
        q=qe,
        m=me,
        w=1.0,
        position_sampler=UniformPositionSampler(low=low, high=high, rng=rng2),
        velocity_sampler=ShiftedMaxwellVelocitySampler(
            drift_velocity=drift,
            thermal_speed=1.0e3,
            rng=rng2,
        ),
    )

    p1 = injector1.sample(32)
    p2 = injector2.sample(32)

    x1 = np.stack([p.x for p in p1], axis=0)
    x2 = np.stack([p.x for p in p2], axis=0)
    v1 = np.stack([p.v for p in p1], axis=0)
    v2 = np.stack([p.v for p in p2], axis=0)

    assert np.allclose(x1, x2)
    assert np.allclose(v1, v2)


def test_fixed_beam_injector_replicates_values() -> None:
    injector = FixedBeamInjector(
        x0=np.array([0.25, 0.25, 0.5]),
        v0=np.array([0.0, 0.0, -2e5]),
        q=-1.0,
        m=2.0,
        w=3.0,
    )
    particles = injector.sample(4)

    assert len(particles) == 4
    assert all(np.allclose(p.x, [0.25, 0.25, 0.5]) for p in particles)
    assert all(np.allclose(p.v, [0.0, 0.0, -2e5]) for p in particles)
    assert all(p.q == -1.0 and p.m == 2.0 and p.w == 3.0 for p in particles)


def test_fixed_beam_injector_rejects_negative_n() -> None:
    injector = FixedBeamInjector(
        x0=np.array([0.25, 0.25, 0.5]),
        v0=np.array([0.0, 0.0, -2e5]),
        q=-1.0,
        m=2.0,
    )

    with pytest.raises(ValueError, match="n must be non-negative"):
        injector.sample(-1)
