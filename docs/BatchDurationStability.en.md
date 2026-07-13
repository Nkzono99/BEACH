title: `batch_duration` Stability and Steady Value

Lang: [English](BatchDurationStability.en.md) | [日本語](BatchDurationStability.md)

# `batch_duration` Stability and Steady Value

This section organizes the theoretical relationship between `sim.batch_duration` in the BEACH batch loop, or the one-batch physical time determined from `sim.batch_duration_step`, and the validity and stability of the converged wall-charge distribution.
In the current implementation, when `batch_duration_step` is used, `sim.batch_duration = sim.dt * sim.batch_duration_step`; in `reservoir_face` injection, the physical inflow per batch determines the macro-particle count or weight.

Implementation entry points:

- Batch procedure: [Computational model overview](Algorithms.en.html)
- Parameter definitions: [Parameters](Parameters.en.html) for `sim.batch_duration` / `sim.batch_duration_step`
- Injection usage: `src/particles/bem_injection.f90` (`reservoir_face` / `photo_raycast`)
- Batch generation and weight resolution: `src/config/bem_app_config_runtime.f90`

## Practical Guide

For first runs, choose `batch_duration` empirically before relying on the theory below.

1. Start with a conservative, small `sim.batch_duration_step`.
2. Inspect `charge_history.csv`, `last_rel_change` in `summary.txt`, and absorbed / escaped counts.
3. Run 2x and 1/2x `batch_duration_step` comparisons and compare final charge distributions and history shapes.
4. If the final charge and history are nearly unchanged and no oscillation or divergence appears, keep that value.
5. If fluctuations dominate, tune `target_macro_particles_per_batch`, `w_particle`, `batch_count`, and `history_stride` in addition to `batch_duration`.

| Symptom | Check | Typical response |
| --- | --- | --- |
| Charge history oscillates strongly by batch | `charge_history.csv` | Lower `batch_duration_step` |
| Final charge changes strongly with step size | 1/2x and 2x comparison | Recompute with smaller `batch_duration_step` |
| History is too noisy to read | `target_macro_particles_per_batch`, `w_particle` | Adjust macro-particle count or weight |
| Run stops before settling | `batches` in `summary.txt` | Increase `batch_count` |

Treat `batch_duration` as the deterministic explicit-update time step, and particle count / weight as the Monte Carlo noise controls.

### 1. Reduction to a continuous-time model

Let $q_j(t)$ be the accumulated charge of insulator wall element `j`, and let $J_j(\mathbf q)$ be the incident charge flux per unit wall area at that charge state.
The absorption-only model becomes:

$$
\frac{dq_j}{dt} \;=\; J_j(\mathbf q)\, A_j
$$

where $A_j$ is the element area. Since $J$ depends on the field created by wall charge, it is generally **nonlinear**.

One BEACH batch can be viewed as an **explicit update that freezes the field at the start of the batch**. In expectation:

$$
\mathbf q^{n+1} \;=\; \mathbf q^n \;+\; \Delta t_b \cdot \mathbf J(\mathbf q^n)\,\mathbf A \;+\; \boldsymbol\eta^n
$$

where:

- $\Delta t_b = $ `sim.batch_duration`
- $\mathbf A$ is the element-area vector
- $\boldsymbol\eta^n$ is Monte Carlo sampling error within the batch

The implementation follows this picture: particles in a batch see the same field $E(\mathbf q^n)$, and the charge delta is applied to the wall at the end of the batch.
Thus `batch_duration` is the time step of this explicit update.

### 2. Validity of the steady value

Write the mean update map as:

$$
F_{\Delta t_b}(\mathbf q) \;=\; \mathbf q \;+\; \Delta t_b\, \mathbf J(\mathbf q)\,\mathbf A
$$

Its fixed point $\mathbf q^{\ast}$ satisfies:

$$
F_{\Delta t_b}(\mathbf q^{\ast}) = \mathbf q^{\ast}
\quad\Longleftrightarrow\quad
\mathbf J(\mathbf q^{\ast}) = 0
$$

Therefore, **the fixed point of the mean model itself does not depend on $\Delta t_b$**.

In this sense, if the iteration converges stably and Monte Carlo error is sufficiently averaged, changing `batch_duration` does not change the targeted continuous-time steady solution.

However, this statement applies only to the **fixed point of the mean model**. Actual runs include:

- finite-sample error per batch
- fluctuation in monitoring quantities used to judge convergence
- residual error from stopping at a finite batch count

So the observed converged value can retain weak `batch_duration` dependence. The safe statement is:

> The mean fixed point of the iteration is independent of `batch_duration`, but finite-sample and finite-time calculations can show small step-size dependence.

### 3. Linear stability

Near a fixed point $\mathbf q^{\ast}$, define perturbations $\delta\mathbf q^n = \mathbf q^n - \mathbf q^{\ast}$.
The linearized mean update is:

$$
\delta \mathbf q^{n+1} \;=\; \bigl(I + \Delta t_b\, M\bigr)\,\delta \mathbf q^n,
\qquad
M_{ij} \;\equiv\; \frac{\partial (J_i A_i)}{\partial q_j}\bigg|_{\mathbf q^{\ast}}
$$

The stability condition for a general multi-degree-of-freedom system is the spectral-radius condition:

$$
\rho\!\left(I + \Delta t_b\, M\right) < 1
$$

For each eigenvalue $\lambda_k$:

$$
|1 + \Delta t_b\, \lambda_k| < 1
$$

This is the essential BEACH stability condition.

As an insulator wall accumulates charge, it tends to attract fewer particles of the same sign, so the dominant eigenvalues of $M$ are expected to be real negative, $\mathrm{Re}(\lambda_k) < 0$.
Only under this **real-negative dominant mode assumption**, using response time scale $\tau_k \equiv 1/|\lambda_k|$, the fastest mode avoids divergence when:

$$
0 \;<\; \Delta t_b \;<\; \frac{2}{|\lambda_{\max}|} \;=\; 2\,\tau_{\min}
$$

and converges monotonically, or overdamped, when:

$$
0 \;<\; \Delta t_b \;<\; \frac{1}{|\lambda_{\max}|} \;=\; \tau_{\min}
$$

Practical interpretation:

- $\Delta t_b < 2\,\tau_{\min}$: non-divergence under the real-negative dominant mode assumption
- $\Delta t_b < \tau_{\min}$: monotone convergence under the same assumption
- in a general coupled system, the precise condition is $\rho(I + \Delta t_b\, M) < 1$

Thus the $2\tau$ / $\tau$ rule is better described as an **explicit-Euler stability guide under a real-negative dominant mode assumption**, not as a strict BEACH CFL condition.

### 4. Relation to Monte Carlo noise

For a one-mode approximation with noise:

$$
\delta q^{n+1} \;=\; \left(1 - \frac{\Delta t_b}{\tau}\right)\,\delta q^n \;+\; \xi^n
$$

the steady variance depends on the variance of $\xi^n$.
The key point is that **the $\Delta t_b$ dependence of $\mathrm{Var}(\xi^n)$ depends on injection normalization**.
BEACH `reservoir_face` has two modes.

#### 4.1 Fixed `w_particle`

When `w_particle` is specified directly, the physical inflow count changes in proportion to $\Delta t_b$, so the expected macro-particle count per batch follows:

$$
N_\text{macro} \;\propto\; \Delta t_b
$$

The shot-noise variance of the batch charge increment can be regarded roughly as proportional to $\Delta t_b$:

$$
\mathrm{Var}(\xi^n) \;\approx\; \alpha\, \Delta t_b
$$

In the limit $\Delta t_b \ll \tau$, the steady variance does not depend strongly on `batch_duration`.

#### 4.2 Fixed `target_macro_particles_per_batch`

When `w_particle` is solved from `target_macro_particles_per_batch`, the weight is determined as in `src/config/bem_app_config_runtime.f90:644`:

$$
w_\text{particle} \;\propto\; \frac{\Gamma\, A\, \Delta t_b}{N_\text{target}}
$$

so the noise dependence on $\Delta t_b$ differs from the simple $\mathrm{Var}(\xi^n) \propto \Delta t_b$ of §4.1.
The macro-particle count is fixed, while the contribution per particle is proportional to $\Delta t_b$.

#### 4.3 Practical interpretation

The useful separation is:

- `batch_duration` mainly controls **deterministic stability**
- the main knobs for statistical noise are **`w_particle` or `target_macro_particles_per_batch`**

In particular, neither of these is generally true:

> Making `batch_duration` smaller always lowers noise.
> Making `batch_duration` larger leaves noise almost unchanged.

The answer depends on injection normalization.

### 5. Physical estimate of $\tau_{\min}$

$\tau_{\min}$ is the fastest effective response time that controls numerical stability.
It is hard to express with one general physical formula because it depends on geometry, potential distribution, upstream distribution function, and injection model.
In practice, estimate two different quantities.

#### 5.1 Charging / sheath relaxation time

A natural estimate uses an effective capacitance $C_\text{eff}$ and effective conductance $G_\text{eff}$:

$$
\tau_\text{charge} \;\sim\; \frac{C_\text{eff}}{G_\text{eff}}
$$

or a typical potential change $\Delta\phi$ and effective current $I_\text{eff}$:

$$
\tau_\text{charge} \;\sim\; \frac{C_\text{eff}\,\Delta\phi}{I_\text{eff}}
$$

This is a relatively slow charging timescale affected by geometry and shielding.

#### 5.2 Inverse plasma frequency

Another fast reference is:

$$
\tau_{pe} \;=\; \omega_{pe}^{-1} \;=\; \sqrt{\frac{\varepsilon_0\, m_e}{n_e\, e^2}}
$$

This is the microscopic fast timescale of an electron plasma and is useful as a reference for how sharply the system can respond.

However, treating $\omega_{pe}^{-1}$ directly as an upper bound on $\tau_{\min}$ is too strong.
It is better viewed as a **fast-side physical reference**. The effective time constant that limits `batch_duration` often comes from $\tau_\text{charge}`, including geometry and incoming-flux limits.

#### 5.3 Practical choice

For $\tau_{\min}$, estimate:

- $\omega_{pe}^{-1}$: microscopic fast reference
- $\tau_\text{charge}$: system-specific charging / sheath relaxation timescale

Then refine with numerical experiments.

> $\omega_{pe}^{-1}$ is only a fast reference; the actual stability limit is set by an effective response time that often includes $\tau_\text{charge}$.

### 6. Practical usage

1. Estimate both $\omega_{pe}^{-1}$ and $\tau_\text{charge}$ as physical scales.
2. Start with a conservative, smaller `batch_duration` if oscillation should be avoided.
3. Compare charge history and monitoring quantities with `batch_duration` multiplied by 1/2 and 2, as a **step-size sensitivity check**.
4. If the converged values nearly agree and no oscillation or divergence is visible, the `batch_duration` is practically adequate.
5. If noise is large, first adjust `w_particle` or `target_macro_particles_per_batch`. Do not try to solve noise only by changing `batch_duration`.
6. Oscillation in `charge_history.csv` `last_rel_change`, or jitter in element charge time series, is a useful diagnostic. This is better called a **step-size sensitivity check** than strict Richardson extrapolation, because it does not assume a power law of the error.

### 7. Summary

| Item | Conclusion |
|---|---|
| Validity of steady value | The fixed point of the mean update does not depend on `batch_duration` |
| Exact stability condition | $\rho(I + \Delta t_b\, M) < 1$ |
| $2\tau$, $\tau$ rules | Explicit-Euler approximate guide under real-negative dominant modes |
| Role of $\omega_{pe}^{-1}$ | Microscopic fast reference, not generally a direct stability upper bound |
| Noise and `batch_duration` | Dependence is set by injection normalization |
| Main noise-reduction knobs | Adjust `w_particle` or `target_macro_particles_per_batch` |
| Practical check | Step-size sensitivity check by varying `batch_duration` |

It is theoretically clean to say that the steady value of the mean model does not depend on how `batch_duration` is chosen.
The general stability condition $\rho(I + \Delta t_b M) < 1$ follows directly from classical stability analysis.
The remaining uncertainty is the value of $\tau_{\min}$ itself, which must be narrowed down with both case-specific physical estimates and numerical experiments.

### Related documents

- [Fortran parameter file specification](Parameters.en.html) — how to set `sim.batch_duration` / `sim.batch_duration_step`
- [Fortran-centered workflow](Workflow.en.html) — batch-loop execution control
- [Computational model overview](Algorithms.en.html) — relation between physical models and numerical methods
