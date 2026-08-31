title: Theory of `batch_duration`

Lang: [English](BatchDurationTheory.en.md) | [日本語](BatchDurationTheory.md)

# Theory of `batch_duration`

This page explains how `batch_duration` affects charging-update stability and the observed steady value.
For the procedure used to choose a value, see [How to choose `batch_duration`](BatchDurationStability.en.html).

In one sentence, `batch_duration` is the explicit update width of the mean charging model. The fixed point of a mean model
that converges stably does not depend on the width, but discrete-update stability, finite-time error, and Monte Carlo error do.

## Mean charging model

Let $q_j$ be the accumulated charge of insulator element $j$, $A_j$ its area, and $J_j(\mathbf q)$ the net surface-charge
flux for the current charge of all elements, $\mathbf q$. $J_j$ includes incident and emitted contributions with their
signs. Define the mean charging rate as

$$
F_j(\mathbf q)=J_j(\mathbf q)A_j.
$$

The continuous-time mean model is then

$$
\frac{dq_j}{dt}=F_j(\mathbf q).
$$

This equation represents the response obtained by sufficiently averaging particle statistics in the field produced by the
same surface charge.

BEACH freezes the field at its batch-start state throughout one batch. One batch can therefore be viewed conceptually as

$$
\mathbf q^{n+1}
=\mathbf q^n+\Delta t_b\,\mathbf F(\mathbf q^n)+\boldsymbol\eta^n,
$$

where $\Delta t_b$ is `batch_duration` and $\boldsymbol\eta^n$ is Monte Carlo error from a finite number of
macro-particles.

## When the fixed point is independent of width

A fixed point $\mathbf q^\ast$ of the mean update satisfies

$$
\mathbf F(\mathbf q^\ast)=\mathbf 0.
$$

When every element has positive area, this is the same condition as $\mathbf J(\mathbf q^\ast)=\mathbf 0$.
Because $\Delta t_b>0$ does not enter the fixed-point equation, the fixed point itself is independent of `batch_duration`
as long as the mean model converges stably to the same fixed point.

An actual run still includes finite-sample $\boldsymbol\eta^n$, a finite physical end time, and nonlinear response.
Consequently, an identical theoretical fixed point does not imply that observed final charge is independent of width.

## Linear stability near a fixed point

Assume that $\mathbf F$ is differentiable near the fixed point and define its Jacobian as

$$
M=\left.\frac{\partial\mathbf F}{\partial\mathbf q}\right|_{\mathbf q^\ast}.
$$

Temporarily neglect Monte Carlo error and set
$\delta\mathbf q^n=\mathbf q^n-\mathbf q^\ast$. The linearized update is

$$
\delta\mathbf q^{n+1}
\simeq (I+\Delta t_b M)\delta\mathbf q^n.
$$

The general linear stability condition is therefore

$$
\rho(I+\Delta t_b M)<1,
$$

where $\rho$ is the spectral radius.

Only when all dominant eigenvalues are real and negative and the fastest response can be represented by $\tau_{\min}$,

$$
\Delta t_b<2\tau_{\min}
$$

is a non-divergence guide, while

$$
\Delta t_b<\tau_{\min}
$$

is a monotone-convergence guide. A single-time-scale rule cannot determine stability with complex eigenvalues, non-normal
coupling, or strong nonlinearity. These guides are not general BEACH CFL conditions.

## Physical time scales

When selecting an initial candidate, the inverse electron plasma frequency $\omega_{pe}^{-1}$ and the charging time from
effective capacitance and conductance,

$$
\tau_\text{charge}\sim\frac{C_\text{eff}}{G_\text{eff}},
$$

are distinct scales. The first characterizes a fast plasma response, while the second estimates the slower change of
surface charge. The actual upper limit also depends on geometry, potential, inflow distribution, and surface response, so
neither scale alone is a general upper bound for `batch_duration`.

## Coupling to Monte Carlo error

The distribution of $\boldsymbol\eta^n$ depends on macro-particle count and weight. For a source whose particle weight
changes with time width, a width comparison includes both time-discretization differences and sampling-variance differences.

Keep the RNG seed and parallel layout fixed in a fixed-width comparison, and record absorbed/escaped counts,
macro-particle count, and weight in addition to surface-charge norms. If noise dominates, improve particle statistics before
drawing a conclusion about linear stability.

## What the adaptive $k\ne0$ condition controls

Adaptive progression accepts a width when the $k\ne0$ potential produced by the difference $\Delta\mathbf q$ between
candidate and batch-start charge satisfies, at every panel centroid $\mathbf x_i$,

$$
\max_i\left|\Delta\phi_{k\ne0}(\mathbf x_i;\Delta\mathbf q)\right|
\le \Phi_{\max}.
$$

$\Phi_{\max}$ is `max_nonzero_mode_potential_step`.

This condition is a trust bound that prevents the frozen $k\ne0$ field from changing too much in one update.
It is not a local-truncation-error estimator and does not establish an order of global accuracy. It also does not control
the $k=0$ update or the applicability range of a response table or particle sampling.

## Assumptions and applicability limits

- The mean model assumes that particle response at the same $\mathbf q$ can be statistically averaged.
- The linear condition applies only where $\mathbf F$ is differentiable near the fixed point of interest.
- Fixed-point analysis alone cannot predict the final state with multiple fixed points, hysteresis, or a time-dependent attractor.
- A stability condition that omits $\boldsymbol\eta^n$ does not guarantee the fluctuation amplitude of a finite-macro-particle run.
- The adaptive $k\ne0$ condition is not a stability condition for the complete system, including matching-plane inner iteration or implicit zero-mode updates.
- This explanation covers v1.0 insulator accumulation and does not extend to unimplemented resistive or dielectric response.

Select the final value with a step-size sensitivity check at the same physical time, not from theoretical scales alone.
Treat execution completion, numerical convergence, and physical validity as separate decisions.

## Related documents

- [How to choose `batch_duration`](BatchDurationStability.en.html) — fixed-width comparison and adaptive configuration
- [BEACH computational cycle](Algorithms.en.html) — update order with a frozen field within each batch
- [Input parameter reference](Parameters.en.html) — configuration contract
- [Validate Results](ValidationGuide.en.html) — numerical convergence and physical validity
