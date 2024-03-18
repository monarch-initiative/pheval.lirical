"""LIRICAL Runner"""

from dataclasses import dataclass
from pathlib import Path

from pheval.runners.runner import PhEvalRunner

from pheval_lirical.post_process.post_process import post_process_results_format
from pheval_lirical.run.run import prepare_lirical_commands, run_lirical_local
from pheval_lirical.tool_specific_configuration_parser import LIRICALToolSpecificConfigurations


@dataclass
class LiricalPhEvalRunner(PhEvalRunner):
    """_summary_"""

    input_dir: Path
    testdata_dir: Path
    tmp_dir: Path
    output_dir: Path
    config_file: Path
    version: str

    def prepare(self):
        """prepare"""
        print("preparing")

    def run(self):
        """run"""
        print("running with lirical")
        config = LIRICALToolSpecificConfigurations.parse_obj(
            self.input_dir_config.tool_specific_configuration_options
        )
        prepare_lirical_commands(
            input_dir=self.input_dir,
            testdata_dir=self.testdata_dir,
            raw_results_dir=self.raw_results_dir,
            tool_input_commands_dir=self.tool_input_commands_dir,
            lirical_version=self.version,
            tool_specific_configurations=config,
            gene_analysis=self.input_dir_config.gene_analysis,
            variant_analysis=self.input_dir_config.variant_analysis,
        )
        run_lirical_local(
            testdata_dir=self.testdata_dir, tool_input_commands_dir=self.tool_input_commands_dir
        )

    def post_process(self):
        """post_process"""
        print("post processing")
        config = LIRICALToolSpecificConfigurations.parse_obj(
            self.input_dir_config.tool_specific_configuration_options
        )
        post_process_results_format(
            raw_results_dir=self.raw_results_dir,
            output_dir=self.output_dir,
            config=config,
            disease_analysis=self.input_dir_config.disease_analysis,
            gene_analysis=self.input_dir_config.gene_analysis,
            variant_analysis=self.input_dir_config.variant_analysis,
        )
