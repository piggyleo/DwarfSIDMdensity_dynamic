from argparse import Namespace
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.run_galaxy_nautilus import default_output_name, numeric_galaxy_slug, validate_output_paths


def output_args(galaxy_slug: str, filename_stem: str, *, galaxy: str = "unused") -> Namespace:
    return Namespace(
        allow_output_galaxy_mismatch=False,
        galaxy_slug=galaxy_slug,
        galaxy=galaxy,
        chain_output=f"outputs/{filename_stem}_chain.csv",
        checkpoint_output=f"outputs/diagnostics/{filename_stem}_sampler.h5",
        figure_output=f"outputs/figures/{filename_stem}_density_profile.png",
        corner_output=f"outputs/figures/{filename_stem}_corner.png",
    )


@pytest.mark.parametrize(
    "filename_stem",
    [
        "willman_1_weighted_mcrscatter",
        "willman_i_weighted_mcrscatter",
        "willman1_weighted_mcrscatter",
    ],
)
def test_validate_output_paths_accepts_willman_aliases(filename_stem: str) -> None:
    validate_output_paths(output_args("willman_1", filename_stem))


@pytest.mark.parametrize(
    "filename_stem",
    [
        "eridanus_ii_weighted",
        "eridanus_2_weighted",
        "eridanus2_weighted",
    ],
)
def test_validate_output_paths_accepts_eridanus_aliases(filename_stem: str) -> None:
    validate_output_paths(output_args("eridanus_2", filename_stem))


def test_validate_output_paths_rejects_another_galaxy_alias() -> None:
    with pytest.raises(ValueError, match="recognized filename alias"):
        validate_output_paths(output_args("eridanus_2", "willman_i_weighted", galaxy="willman1"))


@pytest.mark.parametrize(
    ("galaxy_name", "expected_slug"),
    [
        ("Bootes I", "bootes_1"),
        ("Eridanus II", "eridanus_2"),
        ("Leo IV", "leo_4"),
        ("Willman 1", "willman_1"),
        ("Antlia 2", "antlia_2"),
        ("Leo T", "leo_t"),
    ],
)
def test_numeric_galaxy_slug_uses_arabic_numerals(galaxy_name: str, expected_slug: str) -> None:
    assert numeric_galaxy_slug(galaxy_name) == expected_slug


@pytest.mark.parametrize(
    ("galaxy_slug", "halo_model", "expected_name"),
    [
        ("eridanus_2", "sidm", "eridanus_2_sidm"),
        ("eridanus_ii", "sidm", "eridanus_2_sidm"),
        ("willman_1", "generalized-hernquist", "willman_1_generalized_hernquist"),
    ],
)
def test_default_output_name_uses_halo_model(
    galaxy_slug: str,
    halo_model: str,
    expected_name: str,
) -> None:
    assert default_output_name(galaxy_slug, halo_model) == expected_name
