from __future__ import annotations

import numpy as np


def boris_push(x, v, q, m, dt, E, B):
    """Standard Boris pusher for E+B (E-only works as well)."""
    qm = q / m
    v_minus = v + qm * E * (0.5 * dt)
    t = qm * B * (0.5 * dt)
    t2 = np.dot(t, t)
    s = 2.0 * t / (1.0 + t2)
    v_prime = v_minus + np.cross(v_minus, t)
    v_plus = v_minus + np.cross(v_prime, s)
    v_new = v_plus + qm * E * (0.5 * dt)
    x_new = x + v_new * dt
    return x_new, v_new
