"""A collection of heavier tests to run locally.

Useful for debugging or testing new features.

There are three degrees of testing:
- lightweight: BALK
- medium: BENLDE
- heavy: euro34

!!!!!!!!!!!!!!!!!!!!!!!!IMPORTANT!!!!!!!!!!!!!!!!!!!!!!!!
- You normally want to run these tests one case at a time. E.g.:

    `pytest tests/local_test.py::test_full_run[BALK]`

- Do not run this on Github's CI!
"""

import subprocess
from pathlib import Path

import pytest

SCENARIOS = ["low", "medium", "high"]


def build_request_all(shape: str):
    """Construct a full request for the given shapes."""
    return " ".join(
        [
            f"results/{shape}/{scenario}/totals.png"
            for scenario in SCENARIOS
        ]
    )


@pytest.mark.parametrize("shape", ["BALK", "BENLDE", "euro34"])
def test_full_run(user_path: Path, shape: str):
    """Test a full request of module outputs (using images as proxy)."""
    request = build_request_all(shape)

    assert subprocess.run(
        f"snakemake --use-conda --cores 4 --forceall {request}",
        shell=True,
        check=True,
        cwd=user_path.parent.parent,
    )
    assert subprocess.run(
        f"snakemake --use-conda --cores 4 {request} --report results/{shape}/report.html",
        shell=True,
        check=True,
        cwd=user_path.parent.parent,
    )
    assert subprocess.run(
        f"snakemake --use-conda --cores 4 {request} --rulegraph | dot -Tpng > results/{shape}/rulegraph.png",
        shell=True,
        check=True,
        cwd=user_path.parent.parent,
    )
