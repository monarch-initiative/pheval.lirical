import os
import subprocess
from pathlib import Path

from pheval.utils.file_utils import all_files

from pheval_lirical.config_parser import LiricalConfig
from pheval_lirical.prepare.prepare_commands import prepare_commands


def prepare_lirical_commands(
    config: LiricalConfig, input_dir: Path, output_dir: Path, testdata_dir: Path
):
    Path(output_dir).joinpath(f"lirical_{config.run.version}_{Path(input_dir).name}").mkdir(
        parents=True, exist_ok=True
    )
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
        input_dir=input_dir,
        exomiser_data_dir=config.run.exomiser_configurations.path_to_exomiser_data_directory,
        phenopacket_dir=phenopacket_dir,
        vcf_dir=vcf_dir,
        file_prefix=Path(testdata_dir).name,
        output_dir=Path(output_dir).joinpath(
            f"lirical_{config.run.version}_{Path(input_dir).name}"
        ),
        results_dir=Path(output_dir).joinpath(
            f"lirical_{config.run.version}_{Path(input_dir).name}/{Path(testdata_dir).name}_results/lirical_results"
        ),
    )


def run_lirical_local(config: LiricalConfig, input_dir: Path, testdata_dir: Path, output_dir: Path):
    """Run Phen2Gene locally."""
    Path(output_dir).joinpath(
        f"lirical_{config.run.version}_{Path(input_dir).name}/{Path(testdata_dir).name}_results/lirical_results"
    ).mkdir(parents=True, exist_ok=True)
    batch_file = [
        file
        for file in all_files(
            Path(output_dir).joinpath(
                f"lirical_{config.run.version}_{Path(input_dir).name}/lirical_batch_files"
            )
        )
        if file.name.startswith(os.path.basename(testdata_dir))
    ][0]
    print("running LIRICAL")
    subprocess.run(
        ["bash", str(batch_file)],
        shell=False,
    )
