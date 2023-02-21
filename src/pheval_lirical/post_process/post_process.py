import os
from pathlib import Path

from pheval_lirical.config_parser import LiricalConfig
from pheval_lirical.post_process.post_process_results_format import create_standardised_results


def post_process_results_format(
    input_dir: Path, testdata_dir: Path, output_dir: Path, config: LiricalConfig
):
    """Create pheval gene and variant result from LIRICAL tsv output."""
    print("...creating pheval results format...")
    run_output_dir = Path(output_dir).joinpath(
        f"lirical_{config.run.version}_{Path(input_dir).name}{os.sep}{os.path.basename(testdata_dir)}_results"
    )
    create_standardised_results(
        results_dir=Path(run_output_dir).joinpath("lirical_results"),
        output_dir=run_output_dir,
        sort_order=config.post_process.sort_order,
    )
    print("done")
