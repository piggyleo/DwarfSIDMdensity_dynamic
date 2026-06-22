"""Data loading and Hayashi-equivalent sample selection."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .constants import DEG_TO_RAD, KPC_TO_PC

LOW_RISK_GALAXIES: tuple[str, ...] = (
    "Canes Venatici I",
    "Canes Venatici II",
    "Coma Berenices",
    "Eridanus II",
    "Horologium I",
    "Hydra II",
    "Leo T",
    "Pisces II",
    "Ursa Major I",
    "Ursa Major II",
    "Willman 1",
)


@dataclass(frozen=True)
class GalaxyObservables:
    galaxy: str
    n_sample_table1: int
    distance_kpc: float
    b_star_pc: float
    qprime: float
    center_ra_deg: float
    center_dec_deg: float


@dataclass(frozen=True)
class GalaxyData:
    observables: GalaxyObservables
    stars: pd.DataFrame

    @property
    def x_pc(self) -> np.ndarray:
        return self.stars["x_pc"].to_numpy(float)

    @property
    def y_pc(self) -> np.ndarray:
        return self.stars["y_pc"].to_numpy(float)

    @property
    def velocity_kms(self) -> np.ndarray:
        return self.stars["v_los_kms"].to_numpy(float)

    @property
    def velocity_err_kms(self) -> np.ndarray:
        return self.stars["v_los_err_kms"].to_numpy(float)


def load_galaxy_data(
    galaxy_csv: str | Path,
    *,
    model_member_flags: tuple[int, ...] = (1, 2),
    require_velocity: bool = True,
    structural_centers_csv: str | Path | None = None,
) -> GalaxyData:
    """Load one per-galaxy CSV and select the default modeling sample.

    The project CSV files contain one ``global`` row and many ``star`` rows. The
    preferred selector uses ``model_member_flag`` when present:

    - ``1``: baseline clean modeling member.
    - ``2``: included by default but tagged for sensitivity tests.
    - ``3``: excluded from the simple Jeans baseline.

    Files that have not yet been migrated to ``model_member_flag`` should fail
    fast rather than silently falling back to a historical selector.
    """

    path = Path(galaxy_csv)
    table = pd.read_csv(path, comment="#")
    if "row_kind" not in table:
        raise ValueError(f"{path} is missing row_kind")

    global_rows = table.loc[table["row_kind"] == "global"]
    if global_rows.empty:
        raise ValueError(f"{path} has no global row")
    global_row = global_rows.iloc[0]

    stars = table.loc[table["row_kind"] == "star"].copy()
    if "model_member_flag" not in stars.columns:
        raise ValueError(
            f"{path} is missing model_member_flag; migrate the galaxy table before modeling"
        )
    stars["model_member_flag"] = pd.to_numeric(stars["model_member_flag"], errors="coerce")
    stars = stars[stars["model_member_flag"].isin(model_member_flags)].copy()

    numeric_cols = [
        "ra_deg",
        "dec_deg",
        "v_los_kms",
        "v_los_err_down_kms",
        "v_los_err_up_kms",
    ]
    for col in numeric_cols:
        stars[col] = pd.to_numeric(stars[col], errors="coerce")

    if require_velocity:
        stars = stars.dropna(subset=["ra_deg", "dec_deg", "v_los_kms", "v_los_err_down_kms", "v_los_err_up_kms"])

    stars = stars.copy()
    stars["v_los_err_kms"] = 0.5 * (stars["v_los_err_down_kms"] + stars["v_los_err_up_kms"])

    distance_kpc = float(global_row["De_kpc"])
    center_ra, center_dec = _load_structural_center(
        str(global_row["galaxy"]),
        path,
        structural_centers_csv=structural_centers_csv,
    )
    x_pc, y_pc = _project_tangent_plane_pc(
        stars["ra_deg"].to_numpy(float),
        stars["dec_deg"].to_numpy(float),
        center_ra,
        center_dec,
        distance_kpc,
    )
    stars["x_pc"] = x_pc
    stars["y_pc"] = y_pc

    obs = GalaxyObservables(
        galaxy=str(global_row["galaxy"]),
        n_sample_table1=int(float(global_row["n_sample_table1"])),
        distance_kpc=distance_kpc,
        b_star_pc=float(global_row["b_star_pc"]),
        qprime=float(global_row["qprime"]),
        center_ra_deg=center_ra,
        center_dec_deg=center_dec,
    )
    return GalaxyData(observables=obs, stars=stars.reset_index(drop=True))


def load_low_risk_galaxies(
    data_dir: str | Path,
    *,
    structural_centers_csv: str | Path | None = None,
) -> dict[str, GalaxyData]:
    """Load all low-risk modeling galaxies listed in PROJECT_MEMORY.md."""

    data_path = Path(data_dir)
    loaded: dict[str, GalaxyData] = {}
    for csv_path in sorted(data_path.glob("*.csv")):
        if csv_path.name.startswith("_"):
            continue
        frame = pd.read_csv(csv_path, comment="#", nrows=1)
        if frame.empty or frame.loc[0, "galaxy"] not in LOW_RISK_GALAXIES:
            continue
        galaxy = load_galaxy_data(csv_path, structural_centers_csv=structural_centers_csv)
        loaded[galaxy.observables.galaxy] = galaxy
    return loaded


def _load_structural_center(
    galaxy_name: str,
    galaxy_csv: Path,
    *,
    structural_centers_csv: str | Path | None,
) -> tuple[float, float]:
    centers_path = Path(structural_centers_csv) if structural_centers_csv is not None else _default_centers_path(galaxy_csv)
    centers = pd.read_csv(centers_path, comment="#")
    row = centers.loc[centers["galaxy"] == galaxy_name]
    if row.empty:
        raise ValueError(f"{centers_path} has no structural center for {galaxy_name!r}")
    center = row.iloc[0]
    ra = pd.to_numeric(center["ra_center_deg"], errors="coerce")
    dec = pd.to_numeric(center["dec_center_deg"], errors="coerce")
    if not np.isfinite(ra) or not np.isfinite(dec):
        raise ValueError(f"{centers_path} has an invalid structural center for {galaxy_name!r}")
    return float(ra), float(dec)


def _default_centers_path(galaxy_csv: Path) -> Path:
    data_root = galaxy_csv.resolve().parent.parent
    centers_path = data_root / "processed" / "galaxy_structural_centers.csv"
    if centers_path.exists():
        return centers_path
    return Path("data/processed/galaxy_structural_centers.csv")


def _project_tangent_plane_pc(
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
    center_ra_deg: float,
    center_dec_deg: float,
    distance_kpc: float,
) -> tuple[np.ndarray, np.ndarray]:
    ra0 = center_ra_deg * DEG_TO_RAD
    dec0 = center_dec_deg * DEG_TO_RAD
    ra = ra_deg * DEG_TO_RAD
    dec = dec_deg * DEG_TO_RAD
    distance_pc = distance_kpc * KPC_TO_PC
    x = (ra - ra0) * np.cos(dec0) * distance_pc
    y = (dec - dec0) * distance_pc
    return x, y
