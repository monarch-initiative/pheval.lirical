from pathlib import Path

import click
from packaging import version
from phenopackets import Phenopacket, PhenotypicFeature
from pheval.utils.file_utils import files_with_suffix
from pheval.utils.phenopacket_utils import PhenopacketUtil, phenopacket_reader

from pheval_lirical.prepare.prepare_manual_commands import LiricalManualCommandLineArguments
from pheval_lirical.prepare.prepare_phenopacket_commands import (
    LiricalPhenopacketCommandLineArguments,
)


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
        mode: str,
        exomiser_hg19_data_path: Path,
        exomiser_hg38_data_path: Path,
    ):
        self.phenopacket_path = phenopacket_path
        self.lirical_jar = lirical_jar
        self.input_dir = input_dir
        self.exomiser_data_dir = exomiser_data_dir
        self.vcf_dir = vcf_dir
        self.results_dir = results_dir
        self.mode = mode
        self.exomiser_hg19_data_path = exomiser_hg19_data_path
        self.exomiser_hg38_data_path = exomiser_hg38_data_path
        self.phenopacket_util = PhenopacketUtil(phenopacket)

    def get_list_negated_phenotypic_features(self):
        """Return list of negated HPO ids if there are any present, otherwise return None."""
        return (
            None
            if self.phenopacket_util.negated_phenotypic_features() == []
            else [hpo.type.id for hpo in self.phenopacket_util.negated_phenotypic_features()]
        )

    def get_list_observed_phenotypic_features(self) -> [PhenotypicFeature]:
        """Return list of observed HPO ids."""
        return [hpo.type.id for hpo in self.phenopacket_util.observed_phenotypic_features()]

    def get_vcf_path(self) -> Path:
        """Return the vcf file path."""
        return self.phenopacket_util.vcf_file_data(
            phenopacket_path=self.phenopacket_path, vcf_dir=self.vcf_dir
        ).uri

    def get_vcf_assembly(self) -> str:
        """Return the vcf assembly."""
        return self.phenopacket_util.vcf_file_data(
            phenopacket_path=self.phenopacket_path, vcf_dir=self.vcf_dir
        ).file_attributes["genomeAssembly"]

    def add_manual_cli_arguments(
        self, gene_analysis: bool, variant_analysis: bool
    ) -> LiricalManualCommandLineArguments:
        """Return all CLI arguments to run LIRICAL in manual mode."""
        return LiricalManualCommandLineArguments(
            lirical_jar_file=self.lirical_jar,
            observed_phenotypes=self.get_list_observed_phenotypic_features(),
            negated_phenotypes=self.get_list_negated_phenotypic_features(),
            assembly=self.get_vcf_assembly() if gene_analysis or variant_analysis else None,
            vcf_file_path=self.get_vcf_path() if gene_analysis or variant_analysis else None,
            lirical_data=self.input_dir,
            exomiser_data=self.exomiser_data_dir,
            sample_id=self.phenopacket_util.sample_id(),
            output_dir=self.results_dir,
            output_prefix=self.phenopacket_path.stem,
            exomiser_hg19_data_path=(
                self.exomiser_hg19_data_path if gene_analysis or variant_analysis else None
            ),
            exomiser_hg38_data_path=(
                self.exomiser_hg38_data_path if gene_analysis or variant_analysis else None
            ),
        )

    def add_phenopacket_cli_arguments(
        self, gene_analysis: bool, variant_analysis: bool
    ) -> LiricalPhenopacketCommandLineArguments:
        """Return all CLI arguments to run LIRICAL in phenopacket mode."""
        return LiricalPhenopacketCommandLineArguments(
            lirical_jar_file=self.lirical_jar,
            phenopacket_path=self.phenopacket_path,
            vcf_file_path=self.get_vcf_path() if gene_analysis or variant_analysis else None,
            assembly=self.get_vcf_assembly() if gene_analysis or variant_analysis else None,
            lirical_data=self.input_dir,
            exomiser_data=self.exomiser_data_dir,
            output_dir=self.results_dir,
            output_prefix=self.phenopacket_path.stem,
            exomiser_hg19_data_path=(
                self.exomiser_hg19_data_path if gene_analysis or variant_analysis else None
            ),
            exomiser_hg38_data_path=(
                self.exomiser_hg38_data_path if gene_analysis or variant_analysis else None
            ),
        )

    def add_cli_arguments(
        self, gene_analysis: bool, variant_analysis: bool
    ) -> LiricalManualCommandLineArguments or LiricalManualCommandLineArguments:
        """Return all CLI arguments."""
        if self.mode.lower() == "phenopacket":
            return self.add_phenopacket_cli_arguments(gene_analysis, variant_analysis)
        elif self.mode.lower() == "manual":
            return self.add_manual_cli_arguments(gene_analysis, variant_analysis)


def create_command_arguments(
    phenopacket_dir: Path,
    lirical_jar: Path,
    input_dir: Path,
    exomiser_data_dir: Path,
    vcf_dir: Path,
    output_dir: Path,
    mode: str,
    exomiser_hg19_data: Path,
    exomiser_hg38_data: Path,
    gene_analysis: bool,
    variant_analysis: bool,
) -> list[LiricalManualCommandLineArguments] or list[LiricalPhenopacketCommandLineArguments]:
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
                mode=mode,
                exomiser_hg19_data_path=exomiser_hg19_data,
                exomiser_hg38_data_path=exomiser_hg38_data,
            ).add_cli_arguments(gene_analysis, variant_analysis)
        )
    return commands


class CommandWriter:
    def __init__(self, mode: str, lirical_version: str, output_file: Path):
        self.mode = mode
        self.version = lirical_version
        self.file = open(output_file, "w")

    def write_java_command(
        self,
        command_arguments: LiricalManualCommandLineArguments
        or LiricalPhenopacketCommandLineArguments,
    ) -> None:
        """Write the basic command do run LIRICAL jar file."""
        self.file.write("java" + " -jar " + str(command_arguments.lirical_jar_file))

    def write_mode(self) -> None:
        """Write mode to run LIRICAL"""
        if self.mode.lower() == "phenopacket":
            self.file.write(" P")
        elif self.mode.lower() == "manual":
            self.file.write(" R")

    def write_phenopacket_path(
        self, command_arguments: LiricalPhenopacketCommandLineArguments
    ) -> None:
        """Write the phenopacket path."""
        self.file.write(" --phenopacket " + str(command_arguments.phenopacket_path))

    def write_observed_phenotypic_features(
        self, command_arguments: LiricalManualCommandLineArguments
    ) -> None:
        """Write observed HPO ids to command."""
        self.file.write(" --observed-phenotypes " + ",".join(command_arguments.observed_phenotypes))

    def write_negated_phenotypic_features(
        self, command_arguments: LiricalManualCommandLineArguments
    ) -> None:
        """Write negated HPO ids to command."""
        if command_arguments.negated_phenotypes is not None:
            self.file.write(
                " --negated-phenotypes " + ",".join(command_arguments.negated_phenotypes)
            )

    def write_vcf_file_properties(
        self,
        command_arguments: LiricalManualCommandLineArguments
        or LiricalPhenopacketCommandLineArguments,
    ) -> None:
        """Write related VCF arguments to command."""
        if command_arguments.vcf_file_path is not None:
            self.file.write(
                " --vcf "
                + str(command_arguments.vcf_file_path)
                + " --assembly "
                + command_arguments.assembly
            )

    def write_sample_id(self, command_arguments: LiricalManualCommandLineArguments) -> None:
        """Write the sample id."""
        self.file.write(" --sample-id " + '"' + command_arguments.sample_id + '"')

    def write_lirical_data_dir(
        self,
        command_arguments: LiricalManualCommandLineArguments
        or LiricalPhenopacketCommandLineArguments,
    ) -> None:
        """Write LIRICAL data directory location."""
        self.file.write(" --data " + str(command_arguments.lirical_data))

    def write_exomiser_data_dir(
        self,
        command_arguments: LiricalManualCommandLineArguments
        or LiricalPhenopacketCommandLineArguments,
    ) -> None:
        """Write Exomiser data location, dealing with deprecated parameters."""
        if version.parse(self.version) > version.parse("2.0.0-RC1"):
            if command_arguments.exomiser_hg19_data_path is not None:
                self.file.write(" -e19 " + str(command_arguments.exomiser_hg19_data_path))
            if command_arguments.exomiser_hg38_data_path is not None:
                self.file.write(" -e38 " + str(command_arguments.exomiser_hg38_data_path))
        if version.parse(self.version) < version.parse("2.0.0-RC2"):
            self.file.write(" --exomiser " + str(command_arguments.exomiser_data))

    def write_output_parameters(self, command_arguments: LiricalManualCommandLineArguments) -> None:
        """Write related output parameter arguments to command."""
        self.file.write(
            " --prefix "
            + command_arguments.output_prefix
            + " --output-directory "
            + str(command_arguments.output_dir)
            + " --output-format "
            + "tsv"
        )

    def write_common_arguments(
        self,
        command_arguments: LiricalManualCommandLineArguments
        or LiricalPhenopacketCommandLineArguments,
    ) -> None:
        """Write common CLI parameters."""
        self.write_java_command(command_arguments)
        self.write_mode()
        self.write_vcf_file_properties(command_arguments)
        self.write_lirical_data_dir(command_arguments)
        self.write_exomiser_data_dir(command_arguments)
        self.write_output_parameters(command_arguments)

    def write_manual_command(self, command_arguments: LiricalManualCommandLineArguments) -> None:
        """Write LIRICAL command to file to run in manual mode."""
        self.write_common_arguments(command_arguments)
        self.write_observed_phenotypic_features(command_arguments)
        self.write_negated_phenotypic_features(command_arguments)
        self.write_sample_id(command_arguments)
        self.file.write("\n")

    def write_phenopacket_command(
        self, command_arguments: LiricalPhenopacketCommandLineArguments
    ) -> None:
        """Write LIRICAL command to file to run in phenopacket mode."""
        self.write_common_arguments(command_arguments)
        self.write_phenopacket_path(command_arguments)
        self.file.write("\n")

    def write_command(
        self,
        command_arguments: LiricalManualCommandLineArguments
        or LiricalPhenopacketCommandLineArguments,
    ) -> None:
        """Write LIRICAL command."""
        try:
            (
                self.write_phenopacket_command(command_arguments)
                if self.mode.lower() == "phenopacket"
                else self.write_manual_command(command_arguments)
            )
        except IOError:
            print("Error writing ", self.file)

    def close(self) -> None:
        """Close file."""
        try:
            self.file.close()
        except IOError:
            print("Error closing ", self.file)


def write_all_commands(
    command_arguments: [LiricalManualCommandLineArguments]
    or [LiricalPhenopacketCommandLineArguments],
    tool_input_commands_dir: Path,
    file_prefix: Path,
    mode: str,
    lirical_version: str,
) -> None:
    """Write all commands to file for running LIRICAL."""
    command_writer = CommandWriter(
        mode=mode,
        lirical_version=lirical_version,
        output_file=tool_input_commands_dir.joinpath(f"{file_prefix}-lirical-commands.txt"),
    )
    for command_argument in command_arguments:
        command_writer.write_command(command_argument)
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
    mode: str,
    lirical_version: str,
    exomiser_hg19_data: Path,
    exomiser_hg38_data: Path,
    gene_analysis: bool,
    variant_analysis: bool,
) -> None:
    """Prepare command batch files to run LIRICAL."""
    command_arguments = create_command_arguments(
        phenopacket_dir,
        lirical_jar,
        input_dir,
        exomiser_data_dir,
        vcf_dir,
        raw_results_dir,
        mode,
        exomiser_hg19_data,
        exomiser_hg38_data,
        gene_analysis,
        variant_analysis,
    )
    write_all_commands(
        command_arguments, tool_input_commands_dir, file_prefix, mode, lirical_version
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
@click.option("--output-dir", "-o", required=True, help="Path to output of batch files.", type=Path)
@click.option(
    "--results-dir", "-r", required=True, help="Path to output LIRICAL results.", type=Path
)
@click.option(
    "--mode",
    "-m",
    required=False,
    default="manual",
    help="Mode to run LIRICAL.",
    type=click.Choice(["phenopacket", "manual"]),
)
@click.option("--lirical_version", "-v", required=True, help="Lirical version.", type=str)
@click.option(
    "--exomiser-hg19", "-e19", required=False, help="Exomiser hg19 variant database.", type=Path
)
@click.option(
    "--exomiser-hg38", "-e38", required=False, help="Exomiser hg38 variant database.", type=Path
)
@click.option(
    "--gene-analysis/--no-gene-analysis",
    default=False,
    required=False,
    type=bool,
    show_default=True,
    help="Specify analysis for gene prioritisation",
)
@click.option(
    "--variant-analysis/--no-variant-analysis",
    default=False,
    required=False,
    type=bool,
    show_default=True,
    help="Specify analysis for variant prioritisation",
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
    mode: str,
    lirical_version: str,
    exomiser_hg19: Path,
    exomiser_hg38: Path,
    gene_analysis: bool,
    variant_analysis: bool,
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
        mode,
        lirical_version,
        exomiser_hg19,
        exomiser_hg38,
        gene_analysis,
        variant_analysis,
    )
