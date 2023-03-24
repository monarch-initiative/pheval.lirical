from dataclasses import dataclass
from pathlib import Path

import click
from phenopackets import Phenopacket, PhenotypicFeature
from pheval.utils.file_utils import files_with_suffix
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
    exomiser_data: Path
    output_dir: Path
    output_prefix: str


class CommandCreator:
    def __init__(
        self,
        phenopacket_path: Path,
        phenopacket: Phenopacket,
        lirical_jar: Path,
        input_dir: Path,
        exomiser_data_dir: Path,
        vcf_dir: Path,
        results_dir: Path,
    ):
        self.phenopacket_path = phenopacket_path
        self.lirical_jar = lirical_jar
        self.input_dir = input_dir
        self.exomiser_data_dir = exomiser_data_dir
        self.vcf_dir = vcf_dir
        self.results_dir = results_dir
        self.phenopacket_util = PhenopacketUtil(phenopacket)

    def get_list_negated_phenotypic_features(self):
        """Return list of negated HPO ids if there are any present, otherwise return None."""
        return (
            None
            if self.phenopacket_util.negated_phenotypic_features() == []
            else [hpo.type.id for hpo in self.phenopacket_util.negated_phenotypic_features()]
        )

    def get_list_observed_phenotypic_features(self):
        """Return list of observed HPO ids."""
        return [hpo.type.id for hpo in self.phenopacket_util.observed_phenotypic_features()]

    def get_vcf_path(self):
        """Return the vcf file path."""
        return self.phenopacket_util.vcf_file_data(
            phenopacket_path=self.phenopacket_path, vcf_dir=self.vcf_dir
        ).uri

    def get_vcf_assembly(self):
        """Return the vcf assembly."""
        return self.phenopacket_util.vcf_file_data(
            phenopacket_path=self.phenopacket_path, vcf_dir=self.vcf_dir
        ).file_attributes["genomeAssembly"]

    def add_manual_cli_arguments(self):
        """Return all CLI arguments to run LIRICAL in manual mode."""
        return LiricalManualCommandLineArguments(
            lirical_jar_file=self.lirical_jar,
            observed_phenotypes=self.get_list_observed_phenotypic_features(),
            negated_phenotypes=self.get_list_negated_phenotypic_features(),
            assembly=self.get_vcf_assembly(),
            vcf_file_path=self.get_vcf_path(),
            lirical_data=self.input_dir,
            exomiser_data=self.exomiser_data_dir,
            sample_id=self.phenopacket_util.sample_id(),
            output_dir=self.results_dir,
            output_prefix=self.phenopacket_path.stem,
        )


def create_command_arguments(
    phenopacket_dir: Path,
    lirical_jar: Path,
    input_dir: Path,
    exomiser_data_dir: Path,
    vcf_dir: Path,
    output_dir: Path,
) -> list[LiricalManualCommandLineArguments]:
    """Return a list of LIRICAL command line arguments for a directory of phenopackets."""
    phenopacket_paths = files_with_suffix(phenopacket_dir, ".json")
    commands = []
    for phenopacket_path in phenopacket_paths:
        phenopacket = phenopacket_reader(phenopacket_path)
        commands.append(
            CommandCreator(
                phenopacket_path=phenopacket_path,
                phenopacket=phenopacket,
                lirical_jar=lirical_jar,
                input_dir=input_dir,
                exomiser_data_dir=exomiser_data_dir,
                vcf_dir=vcf_dir,
                results_dir=output_dir,
            ).add_manual_cli_arguments()
        )
    return commands


class CommandWriter:
    def __init__(self, output_file: Path):
        self.file = open(output_file, "w")

    def write_basic_manual_command(self, command_arguments: LiricalManualCommandLineArguments):
        """Write the basic command do run LIRICAL in manual mode."""
        try:
            self.file.write("java" + " -jar " + str(command_arguments.lirical_jar_file) + " R ")
        except IOError:
            print("Error writing ", self.file)

    def write_observed_phenotypic_features(
        self, command_arguments: LiricalManualCommandLineArguments
    ):
        """Write observed HPO ids to command."""
        try:
            self.file.write(
                "--observed-phenotypes " + ",".join(command_arguments.observed_phenotypes)
            )
        except IOError:
            print("Error writing ", self.file)

    def write_negated_phenotypic_features(
        self, command_arguments: LiricalManualCommandLineArguments
    ):
        """Write negated HPO ids to command."""
        try:
            if command_arguments.negated_phenotypes is not None:
                self.file.write(
                    " --negated-phenotypes " + ",".join(command_arguments.negated_phenotypes)
                )
        except IOError:
            print("Error writing ", self.file)

    def write_vcf(self, command_arguments: LiricalManualCommandLineArguments):
        """Write related VCF arguments to command."""
        try:
            self.file.write(
                " --vcf "
                + str(command_arguments.vcf_file_path)
                + " --assembly "
                + command_arguments.assembly
                + " --sample-id "
                + '"'
                + command_arguments.sample_id
                + '"'
            )
        except IOError:
            print("Error writing ", self.file)

    def write_data_dirs(self, command_arguments: LiricalManualCommandLineArguments):
        """Write data directory locations to command."""
        try:
            self.file.write(
                " --data "
                + str(command_arguments.lirical_data)
                + " --exomiser "
                + str(command_arguments.exomiser_data)
            )
        except IOError:
            print("Error writing ", self.file)

    def write_output_parameters(self, command_arguments: LiricalManualCommandLineArguments):
        """Write related output parameter arguments to command."""
        try:
            self.file.write(
                " --prefix "
                + command_arguments.output_prefix
                + " --output-directory "
                + str(command_arguments.output_dir)
                + " --output-format "
                + "tsv"
            )
        except IOError:
            print("Error writing ", self.file)

    def write_manual_command(self, command_arguments: LiricalManualCommandLineArguments):
        """Write LIRICAL command to file to run in manual mode."""
        self.write_basic_manual_command(command_arguments)
        self.write_observed_phenotypic_features(command_arguments)
        self.write_negated_phenotypic_features(command_arguments)
        self.write_vcf(command_arguments)
        self.write_data_dirs(command_arguments)
        self.write_output_parameters(command_arguments)
        self.file.write("\n")

    def close(self) -> None:
        """Close file."""
        try:
            self.file.close()
        except IOError:
            print("Error closing ", self.file)


def write_all_manual_commands(
    command_arguments: [LiricalManualCommandLineArguments],
    tool_input_commands_dir: Path,
    file_prefix: Path,
):
    """Write all commands to file for running LIRICAL in manual mode."""
    command_writer = CommandWriter(
        output_file=tool_input_commands_dir.joinpath(f"{file_prefix}-lirical-commands.txt")
    )
    for command_argument in command_arguments:
        command_writer.write_manual_command(command_argument)
    command_writer.close()


def prepare_commands(
    lirical_jar: Path,
    input_dir: Path,
    exomiser_data_dir: Path,
    phenopacket_dir: Path,
    vcf_dir: Path,
    file_prefix: str,
    tool_input_commands_dir: Path,
    raw_results_dir: Path,
):
    """Prepare command batch files to run LIRICAL."""
    command_arguments = create_command_arguments(
        phenopacket_dir, lirical_jar, input_dir, exomiser_data_dir, vcf_dir, raw_results_dir
    )
    write_all_manual_commands(command_arguments, tool_input_commands_dir, file_prefix)


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
    """Prepare command batch files to run LIRICAL."""
    output_dir.joinpath("tool_input_commands").mkdir(parents=True, exist_ok=True)
    prepare_commands(
        lirical_jar,
        input_dir,
        exomiser_data_dir,
        phenopacket_dir,
        vcf_dir,
        file_prefix,
        output_dir.joinpath("tool_input_commands"),
        results_dir,
    )
