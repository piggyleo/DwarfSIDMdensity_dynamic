#!/usr/bin/env python
"""Compare the local PSIDM-25 density and mass profile to the source model."""

from __future__ import annotations

import argparse
import ast
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for path in (PROJECT_ROOT, PROJECT_ROOT / "src"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from hayashi_jeans.halos import SIDMPSIDM25Halo


def load_reference_model(source_root: Path):
    """Load only the source density class, avoiding unrelated JAX lensing imports."""

    source_path = source_root / "lib" / "SIDM_Parametric_Model_jax.py"
    tree = ast.parse(source_path.read_text(encoding="utf-8"), filename=str(source_path))
    class_node = next(
        node
        for node in tree.body
        if isinstance(node, ast.ClassDef) and node.name == "SIDM_parametric_simple"
    )
    module = ast.Module(body=[class_node], type_ignores=[])
    namespace = {"np": np}
    exec(compile(ast.fix_missing_locations(module), str(source_path), "exec"), namespace)
    return namespace["SIDM_parametric_simple"]()


def enclosed_mass(density, radius_pc: float) -> float:
    value = quad(lambda r: 4.0 * np.pi * r * r * float(density(r)), 0.0, radius_pc, epsrel=2e-9, limit=300)[0]
    return float(value)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source-root",
        default="/Users/wangkaihao/sidm/SIDM_Lensing_Model-main",
    )
    parser.add_argument(
        "--output",
        default=str(PROJECT_ROOT / "outputs/diagnostics/sidm_reference_profile_validation.csv"),
    )
    args = parser.parse_args()

    reference = load_reference_model(Path(args.source_root))
    taus = [0.0, 0.01, 0.1, 0.35, 0.5, 0.9, 1.05, 1.08]
    x_values = np.geomspace(1.0e-5, 1.0e3, 300)
    mass_x_values = np.geomspace(1.0e-4, 1.0e2, 30)
    rows = []

    for tau in taus:
        halo = SIDMPSIDM25Halo(q=1.0, rs0_pc=1.0, rho_s0_msun_pc3=1.0, tau=tau)
        reference_density = np.asarray(reference.get_density(x_values, current_tr=tau, rhoss=1.0, rss=1.0))
        local_density = np.asarray(halo.density_at_ellipsoidal_radius(x_values))
        density_rel = np.abs(local_density - reference_density) / np.maximum(np.abs(reference_density), 1.0e-300)

        mass_rel = []
        for radius in mass_x_values:
            local_mass = enclosed_mass(halo.density_at_ellipsoidal_radius, radius)
            reference_mass = enclosed_mass(
                lambda r: reference.get_density(r, current_tr=tau, rhoss=1.0, rss=1.0),
                radius,
            )
            mass_rel.append(abs(local_mass - reference_mass) / max(abs(reference_mass), 1.0e-300))

        rows.append(
            {
                "tau": tau,
                "density_rel_error_median": float(np.median(density_rel)),
                "density_rel_error_p95": float(np.percentile(density_rel, 95.0)),
                "density_rel_error_max": float(np.max(density_rel)),
                "mass_rel_error_median": float(np.median(mass_rel)),
                "mass_rel_error_p95": float(np.percentile(mass_rel, 95.0)),
                "mass_rel_error_max": float(np.max(mass_rel)),
            }
        )

    result = pd.DataFrame(rows)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)
    print(result.to_string(index=False))
    print(f"wrote {output}")
    if result["density_rel_error_max"].max() >= 1.0e-10:
        raise SystemExit("SIDM density reference validation failed")
    if result["mass_rel_error_max"].max() >= 1.0e-6:
        raise SystemExit("SIDM enclosed-mass reference validation failed")


if __name__ == "__main__":
    main()
