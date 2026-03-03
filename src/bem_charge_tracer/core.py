from __future__ import annotations
from dataclasses import dataclass
from typing import Protocol, Optional, List, Tuple, Iterable
import numpy as np

# -----------------------------
# Constants
# -----------------------------
EPS0 = 8.8541878128e-12
K_COULOMB = 1.0 / (4.0 * np.pi * EPS0)


# -----------------------------
# Geometry / Mesh
# -----------------------------
@dataclass
class BEMElement:
    """Triangular boundary element with total charge q [C]."""
    v0: np.ndarray  # (3,)
    v1: np.ndarray  # (3,)
    v2: np.ndarray  # (3,)
    q: float = 0.0  # [C] element total charge

    def centroid(self) -> np.ndarray:
        return (self.v0 + self.v1 + self.v2) / 3.0

    def area(self) -> float:
        return 0.5 * np.linalg.norm(np.cross(self.v1 - self.v0, self.v2 - self.v0))

    def normal(self) -> np.ndarray:
        n = np.cross(self.v1 - self.v0, self.v2 - self.v0)
        nn = np.linalg.norm(n)
        return n / nn if nn > 0 else n


class BEMMesh:
    """
    Triangular surface mesh with per-element charge.
    Provides collision (segment-triangle) and neighborhood queries.
    """

    def __init__(self, elements: List[BEMElement]):
        self.elements = elements
        self.nelem = len(elements)

        v0 = np.stack([e.v0 for e in elements], axis=0) if self.nelem else np.zeros((0, 3), dtype=float)
        v1 = np.stack([e.v1 for e in elements], axis=0) if self.nelem else np.zeros((0, 3), dtype=float)
        v2 = np.stack([e.v2 for e in elements], axis=0) if self.nelem else np.zeros((0, 3), dtype=float)
        self.v0 = v0
        self.v1 = v1
        self.v2 = v2

        self.centers = np.stack([e.centroid() for e in elements], axis=0) if self.nelem else np.zeros((0, 3), dtype=float)
        self.areas = np.array([e.area() for e in elements], dtype=float) if self.nelem else np.zeros((0,), dtype=float)
        self.normals = np.stack([e.normal() for e in elements], axis=0) if self.nelem else np.zeros((0, 3), dtype=float)

        # Axis-aligned bounding boxes (broad-phase collision)
        self.bb_min = np.minimum(np.minimum(v0, v1), v2) if self.nelem else np.zeros((0, 3), dtype=float)
        self.bb_max = np.maximum(np.maximum(v0, v1), v2) if self.nelem else np.zeros((0, 3), dtype=float)

        # Representative length per element (for heuristics)
        self.h_elem = np.sqrt(np.maximum(self.areas, 0.0)) if self.nelem else np.zeros((0,), dtype=float)

    def charges(self) -> np.ndarray:
        return np.array([e.q for e in self.elements], dtype=float)

    def add_charges(self, dq: np.ndarray) -> None:
        """Add dq [C] to each element. dq shape must be (nelem,)."""
        if dq.shape != (self.nelem,):
            raise ValueError(f"dq must have shape ({self.nelem},), got {dq.shape}")
        for i, d in enumerate(dq):
            if d != 0.0:
                self.elements[i].q += float(d)

    def near_elements_by_center(self, x: np.ndarray, r_switch: float) -> np.ndarray:
        """
        Cheap near query: centroid distance < r_switch (+ element size margin).
        """
        if self.nelem == 0:
            return np.zeros((0,), dtype=int)
        dr = self.centers - x[None, :]
        d2 = np.einsum("ij,ij->i", dr, dr)
        # margin using element scale to avoid missing near triangles
        margin = (0.5 * self.h_elem + 1e-15)
        return np.where(d2 <= (r_switch + margin) ** 2)[0]


# -----------------------------
# Collision
# -----------------------------
@dataclass
class HitInfo:
    mesh_id: int
    elem_idx: int
    t: float
    pos: np.ndarray     # (3,)
    normal: np.ndarray  # (3,)


def _segment_triangle_intersect(p0, p1, v0, v1, v2, eps=1e-12):
    """
    Möller–Trumbore with segment constraint t in [0,1].
    Returns (t, hitpos) or None.
    """
    d = p1 - p0
    e1 = v1 - v0
    e2 = v2 - v0
    h = np.cross(d, e2)
    a = np.dot(e1, h)
    if abs(a) < eps:
        return None
    f = 1.0 / a
    s = p0 - v0
    u = f * np.dot(s, h)
    if u < 0.0 or u > 1.0:
        return None
    q = np.cross(s, e1)
    v = f * np.dot(d, q)
    if v < 0.0 or (u + v) > 1.0:
        return None
    t = f * np.dot(e2, q)
    if t < 0.0 or t > 1.0:
        return None
    hit = p0 + t * d
    return float(t), hit


def find_first_hit_in_mesh(mesh: BEMMesh, p0: np.ndarray, p1: np.ndarray) -> Tuple[Optional[int], Optional[float], Optional[np.ndarray]]:
    """
    Returns (elem_idx, t, hitpos) for earliest intersection, or (None, None, None).
    """
    if mesh.nelem == 0:
        return None, None, None

    seg_min = np.minimum(p0, p1)
    seg_max = np.maximum(p0, p1)
    cand = np.where(
        np.all(mesh.bb_max >= seg_min, axis=1) &
        np.all(mesh.bb_min <= seg_max, axis=1)
    )[0]

    best_t = None
    best_hit = None
    best_idx = None
    for i in cand:
        out = _segment_triangle_intersect(p0, p1, mesh.v0[i], mesh.v1[i], mesh.v2[i])
        if out is None:
            continue
        t, hit = out
        if (best_t is None) or (t < best_t):
            best_t, best_hit, best_idx = t, hit, int(i)

    return best_idx, best_t, best_hit


def find_first_hit(bem_list: List[BEMMesh], p0: np.ndarray, p1: np.ndarray) -> Optional[HitInfo]:
    """
    Checks all meshes and returns earliest hit among them.
    """
    best: Optional[HitInfo] = None
    for mid, mesh in enumerate(bem_list):
        idx, t, hit = find_first_hit_in_mesh(mesh, p0, p1)
        if idx is None:
            continue
        if (best is None) or (t < best.t):
            best = HitInfo(
                mesh_id=mid,
                elem_idx=idx,
                t=t,
                pos=hit,
                normal=mesh.normals[idx].copy() if mesh.nelem else np.zeros(3),
            )
    return best


# -----------------------------
# Field evaluation (Electric / Magnetic)
# -----------------------------
def _triangle_subcentroids(v0, v1, v2, n_sub: int) -> np.ndarray:
    """
    Uniform barycentric subdivision into ~n_sub^2 sub-triangles.
    Returns array of sub-triangle centroids (M, 3).
    Charge per sub-triangle will be q / M.
    """
    if n_sub <= 1:
        return ((v0 + v1 + v2) / 3.0)[None, :]

    e1 = v1 - v0
    e2 = v2 - v0
    pts = []
    n = n_sub
    for i in range(n):
        for j in range(n - i):
            p00 = v0 + (i / n) * e1 + (j / n) * e2
            p10 = v0 + ((i + 1) / n) * e1 + (j / n) * e2
            p01 = v0 + (i / n) * e1 + ((j + 1) / n) * e2
            pts.append((p00 + p10 + p01) / 3.0)

            if j < (n - 1 - i):
                p11 = v0 + ((i + 1) / n) * e1 + ((j + 1) / n) * e2
                pts.append((p10 + p11 + p01) / 3.0)

    return np.asarray(pts, dtype=float)


def _E_point_charges(r: np.ndarray, centers: np.ndarray, charges: np.ndarray, softening: float) -> np.ndarray:
    """
    Electric field at r from point charges located at centers.
    r: (3,) or (...,3)
    centers: (N,3), charges: (N,)
    returns E same shape as r
    """
    r = np.asarray(r, dtype=float)
    if centers.shape[0] == 0:
        return np.zeros_like(r, dtype=float)

    if r.ndim == 1:
        dr = r[None, :] - centers
        r2 = np.einsum("ij,ij->i", dr, dr) + softening * softening
        inv_r3 = 1.0 / (np.sqrt(r2) * r2)
        E = K_COULOMB * np.einsum("i,ij->j", charges * inv_r3, dr)
        return E
    else:
        dr = r[..., None, :] - centers[None, ...]        # (...,N,3)
        r2 = np.einsum("...ij,...ij->...i", dr, dr) + softening * softening
        inv_r3 = 1.0 / (np.sqrt(r2) * r2)
        E = K_COULOMB * np.einsum("...i,...ij->...j", charges[None, ...] * inv_r3, dr)
        return E


class MagneticFieldModel(Protocol):
    def B(self, r: np.ndarray, t: float) -> np.ndarray:
        ...


class ZeroB:
    def B(self, r: np.ndarray, t: float) -> np.ndarray:
        return np.zeros(3, dtype=float)


class UniformB:
    def __init__(self, B0: np.ndarray):
        self.B0 = np.asarray(B0, dtype=float).reshape(3)

    def B(self, r: np.ndarray, t: float) -> np.ndarray:
        return self.B0


class BEMField:
    """
    Electric field from bem_list charges.
    - far: point-charge at element centroid (fast, crude near boundary)
    - hybrid near correction: replace near elements with subdivided surface approximation

    NOTE: far part is still O(Nelem) per query in this prototype.
    """

    def __init__(
        self,
        bem_list: List[BEMMesh],
        use_hybrid: bool = True,
        r_switch_factor: float = 3.0,   # r_switch = factor * h_ref
        n_sub: int = 2,                 # near element subdivision level
        softening_factor: float = 0.1,  # eps = factor * h_ref
    ):
        self.bem_list = bem_list
        self.use_hybrid = use_hybrid
        self.r_switch_factor = float(r_switch_factor)
        self.n_sub = int(n_sub)
        self.softening_factor = float(softening_factor)

        # simple caches for speed (recomputed on demand)
        self._centers_all: Optional[np.ndarray] = None
        self._charges_all: Optional[np.ndarray] = None
        self._softening: Optional[float] = None

    def rebuild_cache(self) -> None:
        centers_all = []
        charges_all = []
        h_refs = []

        for mesh in self.bem_list:
            if mesh.nelem == 0:
                continue
            centers_all.append(mesh.centers)
            charges_all.append(mesh.charges())
            h_refs.append(np.median(mesh.h_elem) if mesh.nelem > 0 else 1.0)

        self._centers_all = np.vstack(centers_all) if centers_all else np.zeros((0, 3), dtype=float)
        self._charges_all = np.concatenate(charges_all) if charges_all else np.zeros((0,), dtype=float)

        h_ref = float(np.median(np.array(h_refs))) if len(h_refs) > 0 else 1.0
        self._softening = self.softening_factor * h_ref

    def E(self, r: np.ndarray) -> np.ndarray:
        """
        Electric field at position r (3,).
        """
        if self._centers_all is None or self._charges_all is None or self._softening is None:
            self.rebuild_cache()

        centers_all = self._centers_all
        charges_all = self._charges_all
        soft = float(self._softening)

        # Far (approx, all elements)
        E_far = _E_point_charges(r, centers_all, charges_all, softening=soft)

        if not self.use_hybrid:
            return E_far

        # Hybrid near correction: for each mesh, near indices -> add (exact_subdiv - approx_centroid)
        E_corr = np.zeros(3, dtype=float)
        for mesh in self.bem_list:
            if mesh.nelem == 0:
                continue

            # pick reference size from this mesh
            h_ref = float(np.median(mesh.h_elem))
            r_switch = self.r_switch_factor * h_ref
            idxs = mesh.near_elements_by_center(r, r_switch=r_switch)
            if idxs.size == 0:
                continue

            # approximate term for those elements (centroid point charges)
            centers = mesh.centers[idxs]
            q = mesh.charges()[idxs]
            E_approx_near = _E_point_charges(r, centers, q, softening=soft)

            # "exact" via subdivided surface (still approximate, but better near boundary)
            E_exact_near = np.zeros(3, dtype=float)
            n_sub = max(1, self.n_sub)
            for ei in idxs:
                qj = float(mesh.elements[int(ei)].q)
                if qj == 0.0:
                    continue
                pts = _triangle_subcentroids(mesh.v0[ei], mesh.v1[ei], mesh.v2[ei], n_sub=n_sub)
                # distribute charge equally to sub-triangles
                qk = (qj / pts.shape[0]) * np.ones((pts.shape[0],), dtype=float)
                E_exact_near += _E_point_charges(r, pts, qk, softening=soft)

            E_corr += (E_exact_near - E_approx_near)

        return E_far + E_corr


def calc_electric_field(x, y, z, bem_list: List[BEMMesh], **kwargs):
    """
    Public interface:
      calc_electric_field(x, y, z, bem_list) -> (Ex, Ey, Ez)

    kwargs forwarded to BEMField:
      use_hybrid, r_switch_factor, n_sub, softening_factor
    """
    field = BEMField(bem_list, **kwargs)
    r = np.stack([np.asarray(x), np.asarray(y), np.asarray(z)], axis=-1)

    # vectorize via loop (prototype). For bulk points, optimize with batch evaluation later.
    if r.ndim == 1:
        E = field.E(r)
        return E[0], E[1], E[2]

    Eout = np.zeros_like(r, dtype=float)
    it = np.nditer(r[..., 0], flags=["multi_index"])
    while not it.finished:
        idx = it.multi_index
        Eout[idx] = field.E(r[idx])
        it.iternext()
    return Eout[..., 0], Eout[..., 1], Eout[..., 2]


# -----------------------------
# Particle / Interaction / Charge model
# -----------------------------
@dataclass
class Particle:
    x: np.ndarray   # (3,)
    v: np.ndarray   # (3,)
    q: float        # [C] particle charge
    m: float        # [kg]
    w: float = 1.0  # superparticle weight
    alive: bool = True


class InteractionModel(Protocol):
    def on_collision(self, p: Particle, hit: HitInfo) -> Tuple[bool, float]:
        """
        Returns:
          alive: continue particle or not
          dq_to_surface: charge [C] to deposit to hit element (already includes weight if desired)
        """
        ...


class AbsorptionInteraction:
    """Absorb particle at collision. Deposit particle charge to surface."""
    def on_collision(self, p: Particle, hit: HitInfo) -> Tuple[bool, float]:
        dq = p.q * p.w
        return False, dq


class SurfaceChargeModel(Protocol):
    def begin_batch(self, nelem_total: int) -> None: ...
    def deposit(self, hit: HitInfo, dq: float) -> None: ...
    def commit_batch(self, bem_list: List[BEMMesh]) -> Tuple[float, float]:
        """
        Returns (norm_dq, norm_q_after) for convergence check.
        """
        ...


class InsulatorSurfaceCharge:
    """Insulator: simply accumulate dq on hit elements; no redistribution."""
    def __init__(self):
        self._dq_per_mesh: List[np.ndarray] = []

    def begin_batch(self, nelem_total: int) -> None:
        # keep interface; buffers are allocated lazily once meshes are known
        self._dq_per_mesh = []

    def _ensure_buffers(self, bem_list: List[BEMMesh]) -> None:
        if self._dq_per_mesh:
            return
        self._dq_per_mesh = [np.zeros(mesh.nelem, dtype=float) for mesh in bem_list]

    def deposit(self, hit: HitInfo, dq: float) -> None:
        self._dq_per_mesh[hit.mesh_id][hit.elem_idx] += dq

    def commit_batch(self, bem_list: List[BEMMesh]) -> Tuple[float, float]:
        self._ensure_buffers(bem_list)

        dq_all = np.concatenate(self._dq_per_mesh) if self._dq_per_mesh else np.zeros((0,), dtype=float)
        norm_dq = float(np.linalg.norm(dq_all))

        for mesh, dq in zip(bem_list, self._dq_per_mesh):
            mesh.add_charges(dq)

        q_all = np.concatenate([mesh.charges() for mesh in bem_list]) if bem_list else np.zeros((0,), dtype=float)
        norm_q = float(np.linalg.norm(q_all))

        # reset buffers for next batch
        self._dq_per_mesh = [np.zeros(mesh.nelem, dtype=float) for mesh in bem_list]
        return norm_dq, norm_q


# -----------------------------
# Boris pusher
# -----------------------------
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


# -----------------------------
# Simulator
# -----------------------------
@dataclass
class SimConfig:
    dt: float
    npcls_per_step: int
    max_step: int
    tol_rel: float = 1e-4
    q_floor: float = 1e-30

    # Field options
    use_hybrid: bool = True
    r_switch_factor: float = 3.0
    n_sub: int = 2
    softening_factor: float = 0.1


class TestParticleSimulator:
    def __init__(
        self,
        bem_list: List[BEMMesh],
        config: SimConfig,
        interaction: Optional[InteractionModel] = None,
        surface_charge: Optional[SurfaceChargeModel] = None,
        B_model: Optional[MagneticFieldModel] = None,
    ):
        self.bem_list = bem_list
        self.cfg = config
        self.interaction = interaction if interaction is not None else AbsorptionInteraction()
        self.surface = surface_charge if surface_charge is not None else InsulatorSurfaceCharge()
        self.B_model = B_model if B_model is not None else ZeroB()

        self.field = BEMField(
            bem_list,
            use_hybrid=config.use_hybrid,
            r_switch_factor=config.r_switch_factor,
            n_sub=config.n_sub,
            softening_factor=config.softening_factor,
        )

        # charge model buffer init
        self.surface.begin_batch(sum(m.nelem for m in bem_list))
        # for v0.1 prototype, allocate buffers once at init
        if isinstance(self.surface, InsulatorSurfaceCharge):
            self.surface._ensure_buffers(bem_list)

    def _rebuild_field_cache(self):
        # called after charge commit
        self.field.rebuild_cache()

    def run(self, particles: Iterable[Particle]) -> dict:
        """
        particles: iterable of Particle (treated as a stream).
        Each particle is advanced for at most max_step timesteps (or until absorbed).

        Batch update interpretation in v0.1:
          npcls_per_step = number of processed particles.
        """
        stats = {
            "processed_particles": 0,
            "absorbed": 0,
            "escaped": 0,
            "batches": 0,
            "last_rel_change": None,
        }

        batch_count = 0

        for p in particles:
            if not p.alive:
                continue

            # advance particle up to max_step timesteps (or until collision absorbs it)
            for _ in range(self.cfg.max_step):
                x0 = p.x.copy()
                v0 = p.v.copy()

                E = self.field.E(x0)
                B = self.B_model.B(x0, t=0.0)  # time input placeholder
                x1, v1 = boris_push(x0, v0, p.q, p.m, self.cfg.dt, E, B)

                # collision check along segment x0->x1
                hit = find_first_hit(self.bem_list, x0, x1)
                if hit is not None:
                    alive, dq = self.interaction.on_collision(p, hit)
                    self.surface.deposit(hit, dq)
                    p.alive = alive
                    if not alive:
                        stats["absorbed"] += 1
                        break
                else:
                    p.x = x1
                    p.v = v1

            stats["processed_particles"] += 1
            batch_count += 1

            # batch commit
            if batch_count >= self.cfg.npcls_per_step:
                norm_dq, norm_q = self.surface.commit_batch(self.bem_list)
                self._rebuild_field_cache()
                rel = norm_dq / max(norm_q, self.cfg.q_floor)
                stats["batches"] += 1
                stats["last_rel_change"] = rel
                batch_count = 0

                # stop condition (change rate)
                if rel < self.cfg.tol_rel:
                    break

        return stats
