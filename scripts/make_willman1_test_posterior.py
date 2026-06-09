#!/usr/bin/env python
"""Create a small Willman 1 halo-sample file for Figure 1 plotting tests.

This is only a plotting fixture. It is not a Hayashi reproduction chain.
Replace it with samples from `run_inference` for scientific comparisons.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="outputs/willman1_test_posterior_samples.csv")
    parser.add_argument("--n-samples", type=int, default=512)
    parser.add_argument("--seed", type=int, default=20260515)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    samples = pd.DataFrame(
        {
            "q_halo": np.clip(rng.normal(1.0, 0.08, args.n_samples), 0.4, 2.5),
            "log10_b_halo_pc": rng.normal(3.2, 0.18, args.n_samples),
            "log10_rho0_msun_pc3": rng.normal(-1.5, 0.16, args.n_samples),
            "alpha": np.clip(rng.normal(1.8, 0.25, args.n_samples), 0.2, 4.0),
            "beta": np.clip(rng.normal(6.4, 0.5, args.n_samples), 2.1, 10.0),
            "gamma": np.clip(rng.normal(1.2, 0.12, args.n_samples), 0.0, 2.0),
        }
    )

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    samples.to_csv(output, index=False)
    print(f"wrote {output}")
    print(f"n_samples={len(samples)}")


if __name__ == "__main__":
    main()
