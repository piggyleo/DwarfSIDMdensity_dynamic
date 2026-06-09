"""Parameter containers for the Hayashi-style model."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class HayashiParameters:
    """Model parameters matching the Hayashi et al. notation.

    Lengths are in pc, density is in Msun / pc^3, velocities are km/s.
    """

    q_halo: float
    b_halo_pc: float
    rho0_msun_pc3: float
    beta_z: float
    alpha: float
    beta: float
    gamma: float
    inclination_rad: float
    systemic_velocity_kms: float

    def halo_kwargs(self) -> dict[str, float]:
        return {
            "q": self.q_halo,
            "b_pc": self.b_halo_pc,
            "rho0_msun_pc3": self.rho0_msun_pc3,
            "alpha": self.alpha,
            "beta": self.beta,
            "gamma": self.gamma,
        }

