from dataclasses import dataclass
from pathlib import Path


@dataclass
class LiricalPhenopacketCommandLineArguments:
    """Minimal arguments required to run LIRICAL in phenopacket mode on the command line."""

    lirical_jar_file: Path
    phenopacket_path: Path
    vcf_file_path: Path
    assembly: str
    lirical_data: Path
    exomiser_data: Path
    output_dir: Path
    output_prefix: str
    exomiser_hg19_data_path: Path = None
    exomiser_hg38_data_path: Path = None
