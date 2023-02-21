"""Lirical Runner"""
from dataclasses import dataclass
from pathlib import Path

from pheval.runners.runner import PhEvalRunner

from pheval_lirical.config_parser import parse_lirical_config
from pheval_lirical.post_process.post_process import post_process_results_format
from pheval_lirical.run.run import prepare_lirical_commands, run_lirical_local


@dataclass
class LiricalPhEvalRunner(PhEvalRunner):
    """_summary_"""

    input_dir: Path
    testdata_dir: Path
    tmp_dir: Path
    output_dir: Path
    config_file: Path

    def prepare(self):
        """prepare"""
        print("preparing")

    def run(self):
        """run"""
        config = parse_lirical_config(self.config_file)
        print("running with lirical")
        prepare_lirical_commands(
            config=config,
            input_dir=self.input_dir,
            output_dir=self.output_dir,
            testdata_dir=self.testdata_dir,
        )
        run_lirical_local(
            config=config,
            input_dir=self.input_dir,
            testdata_dir=self.testdata_dir,
            output_dir=self.output_dir,
        )

    def post_process(self):
        """post_process"""
        print("post processing")
        config = parse_lirical_config(self.config_file)
        post_process_results_format(
            input_dir=self.input_dir,
            testdata_dir=self.testdata_dir,
            output_dir=self.output_dir,
            config=config,
        )
