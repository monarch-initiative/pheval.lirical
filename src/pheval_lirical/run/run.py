import os
import subprocess
from pathlib import Path

from pheval.utils.file_utils import all_files

from pheval_lirical.prepare.prepare_commands import prepare_commands
from pheval_lirical.tool_specific_configurations import (
    ExomiserConfigurations,
    LiricalConfigurations,
)


def prepare_lirical_commands(
    input_dir: Path,
    tool_input_commands_dir: Path,
    raw_results_dir: Path,
    testdata_dir: Path,
    lirical_version: str,
    tool_specific_configurations: LiricalConfigurations,
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
    exomiser_configurations = ExomiserConfigurations(
        **tool_specific_configurations.exomiser_configurations
    )
    prepare_commands(
        lirical_jar=[
            filename
            for filename in all_files(
                input_dir.joinpath(tool_specific_configurations.lirical_software_directory_name)
            )
            if filename.name.endswith(".jar")
        ][0],
        input_dir=input_dir.joinpath("data"),
        exomiser_data_dir=input_dir.joinpath(exomiser_configurations.exomiser_data_path)
        if exomiser_configurations.exomiser_data_path is not None
        else None,
        phenopacket_dir=phenopacket_dir,
        vcf_dir=vcf_dir,
        file_prefix=Path(testdata_dir).name,
        tool_input_commands_dir=tool_input_commands_dir,
        raw_results_dir=raw_results_dir,
        mode=tool_specific_configurations.mode,
        lirical_version=lirical_version,
        exomiser_hg19_data=input_dir.joinpath(exomiser_configurations.exomiser_hg19_data)
        if exomiser_configurations.exomiser_hg19_data is not None
        else None,
        exomiser_hg38_data=input_dir.joinpath(exomiser_configurations.exomiser_hg38_data)
        if exomiser_configurations.exomiser_hg38_data is not None
        else None,
    ),


def run_lirical_local(tool_input_commands_dir: Path, testdata_dir: Path):
    """Run LIRICAL locally."""
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
