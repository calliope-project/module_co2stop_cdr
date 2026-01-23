"""CO2Stop file unzipper."""

import sys
import zipfile
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    snakemake: Any


def unzip_to_path(input_path: str, output_path: str, internal_path: str) -> None:
    """Unzip files from a zip archive."""
    with zipfile.ZipFile(input_path, "r") as zfile:
        try:
            data = zfile.read(f"CO2JRC_OpenFormats/{internal_path}")
        except KeyError as e:
            raise FileNotFoundError(f"File {internal_path!r} not found in zip archive") from e
        with open(output_path, "wb") as i:
            i.write(data)


def main() -> None:
    """Main snakemake process."""
    zip_path = snakemake.input.zipfile
    for name, output_path in snakemake.output.items():
        internal_path = snakemake.params.get(name)
        unzip_to_path(zip_path, output_path, internal_path)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    main()
