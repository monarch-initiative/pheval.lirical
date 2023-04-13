from dataclasses import dataclass
from pathlib import Path

import click
from phenopackets import PhenotypicFeature

from pheval_lirical.prepare.prepare_commands import create_command_arguments, CommandWriter


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
