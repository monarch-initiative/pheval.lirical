from pathlib import Path

from pheval.utils.file_utils import files_with_suffix
from pheval.utils.phenopacket_utils import PhenopacketUtil, phenopacket_reader

from pheval_lirical.prepare.prepare_manual_commands import LiricalManualCommandLineArguments
from phenopackets import Phenopacket

from pheval_lirical.prepare.prepare_phenopacket_commands import LiricalPhenopacketCommandLineArguments
# from packaging import version


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
            mode: str
    ):
        self.phenopacket_path = phenopacket_path
        self.lirical_jar = lirical_jar
        self.input_dir = input_dir
        self.exomiser_data_dir = exomiser_data_dir
        self.vcf_dir = vcf_dir
        self.results_dir = results_dir
        self.mode = mode
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

    def add_phenopacket_cli_arguments(self):
        return LiricalPhenopacketCommandLineArguments(lirical_jar_path=self.lirical_jar,
                                                      phenopacket_path=self.phenopacket_path,
                                                      vcf_file_path=self.get_vcf_path(),
                                                      assembly=self.get_vcf_assembly(),
                                                      lirical_data=self.input_dir,
                                                      exomiser_data_path=self.exomiser_data_dir,
                                                      output_dir=self.results_dir,
                                                      output_prefix=self.phenopacket_path.stem)

    def add_cli_arguments(self):
        if self.mode.lower() == "phenopacket":
            return self.add_phenopacket_cli_arguments()
        elif self.mode.lower() == "manual":
            return self.add_manual_cli_arguments()


def create_command_arguments(
        phenopacket_dir: Path,
        lirical_jar: Path,
        input_dir: Path,
        exomiser_data_dir: Path,
        vcf_dir: Path,
        output_dir: Path,
        mode: str
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
                mode=mode
            ).add_cli_arguments()
        )
    return commands


class CommandWriter:
    def __init__(self, mode: str, output_file: Path):
        self.mode = mode
        self.file = open(output_file, "w")

    def write_java_command(self,
                           command_arguments: LiricalManualCommandLineArguments or LiricalPhenopacketCommandLineArguments):
        """Write the basic command do run LIRICAL in manual mode."""
        self.file.write("java" + " -jar " + str(command_arguments.lirical_jar_file))

    def write_mode(self):
        if self.mode.lower() == "phenopacket":
            self.file.write(" P ")
        elif self.mode.lower() == "manual":
            self.file.write(" R ")

    def write_observed_phenotypic_features(
            self, command_arguments: LiricalManualCommandLineArguments
    ):
        """Write observed HPO ids to command."""
        self.file.write(
            "--observed-phenotypes " + ",".join(command_arguments.observed_phenotypes)
        )

    def write_negated_phenotypic_features(
            self, command_arguments: LiricalManualCommandLineArguments
    ):
        """Write negated HPO ids to command."""
        if command_arguments.negated_phenotypes is not None:
            self.file.write(
                " --negated-phenotypes " + ",".join(command_arguments.negated_phenotypes)
            )

    def write_vcf_file_properties(self,
                                  command_arguments: LiricalManualCommandLineArguments or LiricalPhenopacketCommandLineArguments):
        """Write related VCF arguments to command."""
        self.file.write(
            " --vcf "
            + str(command_arguments.vcf_file_path)
            + " --assembly "
            + command_arguments.assembly
        )

    def write_sample_id(self, command_arguments: LiricalManualCommandLineArguments):
        self.file.write(" --sample-id "
                        + '"'
                        + command_arguments.sample_id
                        + '"')

    def write_lirical_data_dir(self,
                               command_arguments: LiricalManualCommandLineArguments or LiricalPhenopacketCommandLineArguments):
        """Write data directory locations to command."""
        self.file.write(
            " --data "
            + str(command_arguments.lirical_data)
            + " --exomiser "
            + str(command_arguments.exomiser_data)
        )

    # def write_exomiser_data_dir(self, command_arguments: LiricalManualCommandLineArguments or LiricalPhenopacketCommandLineArguments):
    #     if version.parse(config.run.exomiser_configurations.exomiser_version) > version.parse("1.3.4"):
    #         self.file.write(" -e19 " +)

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
        self.write_java_command(command_arguments)
        self.write_observed_phenotypic_features(command_arguments)
        self.write_negated_phenotypic_features(command_arguments)
        self.write_vcf_file_properties(command_arguments)
        self.write_lirical_data_dir(command_arguments)
        self.write_output_parameters(command_arguments)
        self.file.write("\n")

    def close(self) -> None:
        """Close file."""
        try:
            self.file.close()
        except IOError:
            print("Error closing ", self.file)
