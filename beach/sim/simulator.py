from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from .collision import find_first_hit
from .field import BEMField, MagneticFieldModel, ZeroB
from .mesh import BEMMesh
from .particles import (
    AbsorptionInteraction,
    InsulatorSurfaceCharge,
    InteractionModel,
    SurfaceChargeModel,
)
from .pusher import boris_push


@dataclass
class SimConfig:
    dt: float
    npcls_per_step: int
    max_step: int
    tol_rel: float = 1e-4
    q_floor: float = 1e-30
    use_hybrid: bool = True
    r_switch_factor: float = 3.0
    n_sub: int = 2
    softening_factor: float = 0.1


class TestParticleSimulator:
    def __init__(
        self,
        bem_list: list[BEMMesh],
        config: SimConfig,
        interaction: Optional[InteractionModel] = None,
        surface_charge: Optional[SurfaceChargeModel] = None,
        B_model: Optional[MagneticFieldModel] = None,
    ):
        self.bem_list = bem_list
        self.cfg = config
        self.interaction = (
            interaction if interaction is not None else AbsorptionInteraction()
        )
        self.surface = (
            surface_charge if surface_charge is not None else InsulatorSurfaceCharge()
        )
        self.B_model = B_model if B_model is not None else ZeroB()

        self.field = BEMField(
            bem_list,
            use_hybrid=config.use_hybrid,
            r_switch_factor=config.r_switch_factor,
            n_sub=config.n_sub,
            softening_factor=config.softening_factor,
        )

        self.surface.begin_batch(sum(m.nelem for m in bem_list))
        if isinstance(self.surface, InsulatorSurfaceCharge):
            self.surface._ensure_buffers(bem_list)

    def _rebuild_field_cache(self):
        self.field.rebuild_cache()

    def run(self, particles) -> dict:
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

            for _ in range(self.cfg.max_step):
                x0 = p.x.copy()
                v0 = p.v.copy()

                E = self.field.E(x0)
                B = self.B_model.B(x0, t=0.0)
                x1, v1 = boris_push(x0, v0, p.q, p.m, self.cfg.dt, E, B)

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

            if batch_count >= self.cfg.npcls_per_step:
                norm_dq, norm_q = self.surface.commit_batch(self.bem_list)
                self._rebuild_field_cache()
                rel = norm_dq / max(norm_q, self.cfg.q_floor)
                stats["batches"] += 1
                stats["last_rel_change"] = rel
                batch_count = 0

                if rel < self.cfg.tol_rel:
                    break

        return stats
