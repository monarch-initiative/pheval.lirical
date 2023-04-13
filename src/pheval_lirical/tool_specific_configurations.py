from dataclasses import dataclass
from pathlib import Path


@dataclass
class ExomiserConfigurations:
    exomiser_data_path: Path
    exomiser_hg19_data: Path
    exomiser_hg38_data: Path


@dataclass
class LiricalConfigurations:
    mode: str
    lirical_software_directory_name: Path
    exomiser_configurations: ExomiserConfigurations
