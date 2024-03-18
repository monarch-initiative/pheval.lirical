import subprocess
from pathlib import Path

from pheval.utils.file_utils import all_files

from pheval_lirical.prepare.prepare_commands import prepare_commands
from pheval_lirical.tool_specific_configuration_parser import LIRICALToolSpecificConfigurations


def prepare_lirical_commands(
    input_dir: Path,
    tool_input_commands_dir: Path,
    raw_results_dir: Path,
    testdata_dir: Path,
    lirical_version: str,
    tool_specific_configurations: LIRICALToolSpecificConfigurations,
    gene_analysis: bool,
    variant_analysis: bool,
):
    """Write commands to run LIRICAL."""
    phenopacket_dir = Path(testdata_dir).joinpath("phenopackets")
    vcf_dir = Path(testdata_dir).joinpath("vcf") if gene_analysis or variant_analysis else None
    prepare_commands(
        lirical_jar=input_dir.joinpath(tool_specific_configurations.lirical_jar_executable),
        input_dir=input_dir.joinpath("data"),
        exomiser_data_dir=(
            input_dir.joinpath(
                tool_specific_configurations.exomiser_db_configurations.exomiser_database
            )
            if tool_specific_configurations.exomiser_db_configurations.exomiser_database is not None
            else None
        ),
        phenopacket_dir=phenopacket_dir,
        vcf_dir=vcf_dir,
        file_prefix=Path(testdata_dir).name,
        tool_input_commands_dir=tool_input_commands_dir,
        raw_results_dir=raw_results_dir,
        mode=tool_specific_configurations.mode,
        lirical_version=lirical_version,
        exomiser_hg19_data=(
            input_dir.joinpath(
                tool_specific_configurations.exomiser_db_configurations.exomiser_hg19_database
            )
            if tool_specific_configurations.exomiser_db_configurations.exomiser_hg19_database
            else None
        ),
        exomiser_hg38_data=(
            input_dir.joinpath(
                tool_specific_configurations.exomiser_db_configurations.exomiser_hg38_database
            )
            if tool_specific_configurations.exomiser_db_configurations.exomiser_hg38_database
            is not None
            else None
        ),
        gene_analysis=gene_analysis,
        variant_analysis=variant_analysis,
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
