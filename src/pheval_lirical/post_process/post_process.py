from pathlib import Path

from pheval_lirical.post_process.post_process_results_format import create_standardised_results
from pheval_lirical.tool_specific_configuration_parser import LIRICALToolSpecificConfigurations


def post_process_results_format(
    raw_results_dir: Path,
    output_dir: Path,
    phenopacket_dir: Path,
    config: LIRICALToolSpecificConfigurations,
    disease_analysis: bool,
    gene_analysis: bool,
    variant_analysis: bool,
):
    """Create pheval gene and variant result from LIRICAL tsv output."""
    print("...creating pheval results format...")
    create_standardised_results(
        raw_results_dir=raw_results_dir,
        output_dir=output_dir,
        phenopacket_dir=phenopacket_dir,
        sort_order=config.post_process.sort_order,
        disease_analysis=disease_analysis,
        gene_analysis=gene_analysis,
        variant_analysis=variant_analysis,
    )
    print("done")
