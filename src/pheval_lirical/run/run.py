import os
import subprocess
from pathlib import Path

from pheval.utils.file_utils import all_files

from pheval_lirical.config_parser import LiricalConfig
from pheval_lirical.prepare.prepare_commands import prepare_commands


def prepare_lirical_commands(
    config: LiricalConfig,
    input_dir: Path,
    tool_input_commands_dir: Path,
    raw_results_dir: Path,
    testdata_dir: Path,
):
    """Write commands to run LIRICAL."""
    phenopacket_dir = Path(testdata_dir).joinpath(
        [
            directory
            for directory in os.listdir(str(testdata_dir))
            if "phenopacket" in str(directory)
        ][0]
    )
    vcf_dir = Path(testdata_dir).joinpath(
        [directory for directory in os.listdir(str(testdata_dir)) if "vcf" in str(directory)][0]
    )
    prepare_commands(
        lirical_jar=config.run.path_to_lirical_software_directory,
        input_dir=input_dir.joinpath("data"),
        exomiser_data_dir=config.run.exomiser_configurations.path_to_exomiser_data_directory,
        phenopacket_dir=phenopacket_dir,
        vcf_dir=vcf_dir,
        file_prefix=Path(testdata_dir).name,
        tool_input_commands_dir=tool_input_commands_dir,
        raw_results_dir=raw_results_dir,
    ),


def run_lirical_local(tool_input_commands_dir: Path, testdata_dir: Path):
    """Run Phen2Gene locally."""
    batch_file = [
        file
        for file in all_files(Path(tool_input_commands_dir))
        if file.name.startswith(Path(testdata_dir).name)
    ][0]
    print("running LIRICAL")
    subprocess.run(
        ["bash", str(batch_file)],
        shell=False,
    )
