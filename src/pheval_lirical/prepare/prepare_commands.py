from dataclasses import dataclass
from pathlib import Path

import click
from phenopackets import Phenopacket, PhenotypicFeature
from pheval.utils.file_utils import all_files
from pheval.utils.phenopacket_utils import PhenopacketUtil, phenopacket_reader


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
    exomiser_data: Path or None
    output_dir: Path
    output_prefix: str


def obtain_negated_phenotypes(phenopacket: Phenopacket) -> PhenotypicFeature:
    # TODO move to PhEval
    """Obtain negated phenotypic features from a Phenopacket."""
    negated_phenotypic_features = []
    phenotypes = PhenopacketUtil(phenopacket).phenotypic_features()
    for phenotype in phenotypes:
        if phenotype.excluded:
            negated_phenotypic_features.append(phenotype)
    return negated_phenotypic_features if negated_phenotypic_features != [] else None


def obtain_sample_id(phenopacket: Phenopacket):
    # TODO move to PhEval
    """Obtain sample ID from a Phenopacket."""
    return phenopacket.subject.id


def create_command_line_arguments(
    phenopacket: Phenopacket,
    lirical_jar: Path,
    input_dir: Path,
    exomiser_data_dir: Path,
    phenopacket_path: Path,
    vcf_dir: Path,
    results_dir: Path,
) -> LiricalManualCommandLineArguments:
    """Create manual command line arguments to run LIRICAL in manual mode."""
    phenopacket_util = PhenopacketUtil(phenopacket)
    vcf_file_data = phenopacket_util.vcf_file_data(
        phenopacket_path=phenopacket_path, vcf_dir=vcf_dir
    )
    return LiricalManualCommandLineArguments(
        lirical_jar_file=lirical_jar,
        observed_phenotypes=[
            hpo.type.id for hpo in phenopacket_util.observed_phenotypic_features()
        ],
        negated_phenotypes=None
        if obtain_negated_phenotypes(phenopacket) is None
        else [hpo.type.id for hpo in obtain_negated_phenotypes(phenopacket)],
        assembly=vcf_file_data.file_attributes["genomeAssembly"],
        vcf_file_path=vcf_file_data.uri,
        lirical_data=input_dir,
        exomiser_data=exomiser_data_dir,
        sample_id=obtain_sample_id(phenopacket),
        output_dir=results_dir,
        output_prefix=phenopacket_path.stem,
    )


class CommandWriter:
    def __init__(self, output_file: Path):
        self.file = open(output_file, "w")

    def write_local_command(
        self, command_line_arguments: LiricalManualCommandLineArguments
    ) -> None:
        """Write command to run LIRICAL in manual mode locally."""
        try:
            self.file.write(
                "java"
                + " -jar "
                + str(command_line_arguments.lirical_jar_file)
                + " R "
                + "--observed-phenotypes "
                + ",".join(command_line_arguments.observed_phenotypes)
            )
            if command_line_arguments.negated_phenotypes is not None:
                self.file.write(
                    " --negated-phenotypes " + command_line_arguments.negated_phenotypes
                )
            self.file.write(
                " --vcf "
                + str(command_line_arguments.vcf_file_path)
                + " --assembly "
                + command_line_arguments.assembly
                + " --data "
                + str(command_line_arguments.lirical_data)
            )
            if command_line_arguments.exomiser_data is not None:
                self.file.write(" --exomiser " + str(command_line_arguments.exomiser_data))
            self.file.write(
                " --sample-id "
                + '"'
                + command_line_arguments.sample_id
                + '"'
                + " --prefix "
                + command_line_arguments.output_prefix
                + " --output-directory "
                + str(command_line_arguments.output_dir)
                + " --output-format "
                + "tsv"
                "\n"
            )
        except IOError:
            print("Error writing ", self.file)

    def close(self) -> None:
        """Close file."""
        try:
            self.file.close()
        except IOError:
            print("Error closing ", self.file)


def write_single_command(
    lirical_jar,
    input_dir: Path,
    exomiser_data_dir: Path,
    phenopacket_path: Path,
    vcf_dir: Path,
    command_writer: CommandWriter,
    results_dir,
) -> None:
    """Write a single command for LIRICAL to run in manual mode."""
    phenopacket = phenopacket_reader(phenopacket_path)
    arguments = create_command_line_arguments(
        lirical_jar=lirical_jar,
        input_dir=input_dir,
        exomiser_data_dir=exomiser_data_dir,
        phenopacket_path=phenopacket_path,
        phenopacket=phenopacket,
        vcf_dir=vcf_dir,
        results_dir=results_dir,
    )
    command_writer.write_local_command(arguments)


def write_local_commands(
    lirical_jar,
    input_dir,
    exomiser_data_dir,
    phenopacket_dir,
    vcf_dir,
    command_file_path,
    results_dir,
) -> None:
    """Write commands to run LIRICAL in manual mode."""
    command_writer = CommandWriter(command_file_path)
    for phenopacket_path in all_files(phenopacket_dir):
        write_single_command(
            lirical_jar,
            input_dir,
            exomiser_data_dir,
            phenopacket_path,
            vcf_dir,
            command_writer,
            results_dir,
        )
    command_writer.close()


def prepare_commands(
    lirical_jar: Path,
    input_dir: Path,
    exomiser_data_dir: Path,
    phenopacket_dir: Path,
    vcf_dir: Path,
    file_prefix: str,
    output_dir: Path,
    results_dir: Path,
):
    """Prepare command batch files to run LIRICAL."""
    output_dir.joinpath("lirical_batch_files").mkdir(parents=True, exist_ok=True)
    command_file_path = output_dir.joinpath(f"lirical_batch_files/{file_prefix}-lirical-batch.txt")
    write_local_commands(
        lirical_jar,
        input_dir,
        exomiser_data_dir,
        phenopacket_dir,
        vcf_dir,
        command_file_path,
        results_dir,
    )


@click.command("prepare-commands")
@click.option("--lirical-jar", "-l", required=True, help="Path to Lirical jar file.", type=Path)
@click.option("--input-dir", "-i", required=True, help="Path to Lirical data directory.", type=Path)
@click.option(
    "--exomiser-data-dir",
    "-exomiser",
    required=False,
    help="Path to exomiser data directory.",
    type=Path,
)
@click.option("--phenopacket-dir", "-p", required=True, help="Path to phenopacket.", type=Path)
@click.option("--vcf-dir", "-v", required=True, help="Path to vcf directory.", type=Path)
@click.option("--file-prefix", "-f", required=True, help="File prefix", type=Path)
@click.option("-output-dir", "-o", required=True, help="Path to output of batch files.", type=Path)
@click.option(
    "-results-dir", "-r", required=True, help="Path to output LIRICAL results.", type=Path
)
def prepare_commands_command(
    lirical_jar: Path,
    input_dir: Path,
    exomiser_data_dir: Path,
    phenopacket_dir: Path,
    vcf_dir: Path,
    file_prefix: str,
    output_dir: Path,
    results_dir: Path,
):
    prepare_commands(
        lirical_jar,
        input_dir,
        exomiser_data_dir,
        phenopacket_dir,
        vcf_dir,
        file_prefix,
        output_dir,
        results_dir,
    )
