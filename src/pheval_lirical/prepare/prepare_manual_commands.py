from dataclasses import dataclass
from pathlib import Path

from phenopackets import PhenotypicFeature


@dataclass
class LiricalManualCommandLineArguments:
    """Minimal arguments required to run LIRICAL manually on the command line."""

    lirical_jar_file: Path
    observed_phenotypes: [PhenotypicFeature]
    negated_phenotypes: [PhenotypicFeature] or None
    assembly: str
    vcf_file_path: Path
    sample_id: str
    lirical_data: Path
    exomiser_data: Path
    output_dir: Path
    output_prefix: str
    exomiser_hg19_data_path: Path = None
    exomiser_hg38_data_path: Path = None
