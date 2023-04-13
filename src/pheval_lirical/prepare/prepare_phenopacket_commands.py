from dataclasses import dataclass
from pathlib import Path

from pheval.utils.phenopacket_utils import PhenopacketUtil
from phenopackets import Phenopacket


@dataclass
class LiricalPhenopacketCommandLineArguments:
    lirical_jar_path: Path
    phenopacket_path: Path
    vcf_file_path: Path
    assembly: str
    lirical_data: Path
    exomiser_data_path: Path
    output_dir: Path
    output_prefix: str
