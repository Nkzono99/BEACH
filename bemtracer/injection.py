from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Protocol

import numpy as np

from .core import Particle

KB = 1.380649e-23  # [J/K]


class PositionSampler(Protocol):
    """Protocol for sampling particle positions [m]."""

    def sample(self, n: int) -> np.ndarray:
        """Return positions with shape (n, 3)."""
        ...


class VelocitySampler(Protocol):
    """Protocol for sampling particle velocities [m/s]."""

    def sample(self, n: int, m: float) -> np.ndarray:
        """Return velocities with shape (n, 3)."""
        ...


@dataclass
class UniformPositionSampler:
    """Uniform position sampler in an axis-aligned box [m]."""

    low: np.ndarray
    high: np.ndarray
    rng: Optional[np.random.Generator] = None

    def __post_init__(self) -> None:
        self.low = np.asarray(self.low, dtype=float)
        self.high = np.asarray(self.high, dtype=float)
        if self.low.shape != (3,) or self.high.shape != (3,):
            raise ValueError("low/high must have shape (3,)")
        if np.any(self.high < self.low):
            raise ValueError("high must be >= low for all axes")
        if self.rng is None:
            self.rng = np.random.default_rng()

    def sample(self, n: int) -> np.ndarray:
        if n < 0:
            raise ValueError("n must be non-negative")
        return self.rng.uniform(self.low, self.high, size=(n, 3))


@dataclass
class ShiftedMaxwellVelocitySampler:
    """Shifted Maxwellian velocity sampler [m/s].

    Parameters
    ----------
    drift_velocity:
        Mean drift velocity vector [m/s].
    temperature_K:
        Isotropic temperature [K]. If set, thermal speed is computed as
        sigma = sqrt(k_B T / m) for each component.
    thermal_speed:
        Optional direct standard deviation per component [m/s]. If provided,
        this is used instead of temperature_K.
    rng:
        Optional random number generator.
    """

    drift_velocity: np.ndarray
    temperature_K: Optional[float] = None
    thermal_speed: Optional[float] = None
    rng: Optional[np.random.Generator] = None

    def __post_init__(self) -> None:
        self.drift_velocity = np.asarray(self.drift_velocity, dtype=float)
        if self.drift_velocity.shape != (3,):
            raise ValueError("drift_velocity must have shape (3,)")
        if self.temperature_K is None and self.thermal_speed is None:
            raise ValueError("either temperature_K or thermal_speed must be provided")
        if self.temperature_K is not None and self.temperature_K < 0:
            raise ValueError("temperature_K must be >= 0")
        if self.thermal_speed is not None and self.thermal_speed < 0:
            raise ValueError("thermal_speed must be >= 0")
        if self.rng is None:
            self.rng = np.random.default_rng()

    def _sigma(self, m: float) -> float:
        if m <= 0:
            raise ValueError("particle mass m must be > 0")
        if self.thermal_speed is not None:
            return float(self.thermal_speed)
        assert self.temperature_K is not None
        return float(np.sqrt(KB * self.temperature_K / m))

    def sample(self, n: int, m: float) -> np.ndarray:
        if n < 0:
            raise ValueError("n must be non-negative")
        sigma = self._sigma(m)
        vth = self.rng.normal(loc=0.0, scale=sigma, size=(n, 3))
        return vth + self.drift_velocity[None, :]


@dataclass
class RandomBeamInjector:
    """Composable injector that creates Particle objects from samplers.

    Units: x[m], v[m/s], q[C], m[kg], w[-].
    """

    q: float
    m: float
    w: float
    position_sampler: PositionSampler
    velocity_sampler: VelocitySampler

    def sample(self, n: int) -> list[Particle]:
        x = self.position_sampler.sample(n)
        v = self.velocity_sampler.sample(n, self.m)
        return [
            Particle(x=x[i].copy(), v=v[i].copy(), q=self.q, m=self.m, w=self.w)
            for i in range(n)
        ]


@dataclass
class FixedBeamInjector:
    """Deterministic injector with fixed x/v for all particles."""

    x0: np.ndarray
    v0: np.ndarray
    q: float
    m: float
    w: float = 1.0

    def __post_init__(self) -> None:
        self.x0 = np.asarray(self.x0, dtype=float)
        self.v0 = np.asarray(self.v0, dtype=float)
        if self.x0.shape != (3,) or self.v0.shape != (3,):
            raise ValueError("x0 and v0 must have shape (3,)")

    def sample(self, n: int) -> list[Particle]:
        return [
            Particle(x=self.x0.copy(), v=self.v0.copy(), q=self.q, m=self.m, w=self.w)
            for _ in range(n)
        ]
