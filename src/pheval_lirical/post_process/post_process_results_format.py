from pathlib import Path

import polars as pl
from pheval.post_processing.phenopacket_truth_set import calculate_end_pos
from pheval.post_processing.post_processing import (
    SortOrder,
    generate_disease_result,
    generate_gene_result,
    generate_variant_result,
)
from pheval.utils.file_utils import files_with_suffix
from pheval.utils.phenopacket_utils import (
    GeneIdentifierUpdater,
    create_gene_identifier_map,
)


def read_lirical_result(lirical_result_path: Path) -> pl.DataFrame:
    """Read LIRICAL tsv output and return a dataframe."""
    return pl.read_csv(lirical_result_path, separator="\t", comment_prefix="!")


def extract_disease_results(raw_result: pl.DataFrame) -> pl.DataFrame:
    """
    Extract disease results from LIRICAL results.
    Args:
        raw_result (pl.DataFrame): LIRICAL results dataframe.
    Returns:
        pl.DataFrame: The extracted results.
    """
    return raw_result.select(
        [
            pl.col("diseaseCurie").alias("disease_identifier"),
            pl.when(pl.col("compositeLR").cast(pl.Utf8) == "-∞")
            .then(float("-inf"))
            .otherwise(pl.col("compositeLR").cast(pl.Float64))
            .alias("score"),
        ]
    )


def extract_gene_results(
    raw_result: pl.DataFrame, gene_identifier_updater: GeneIdentifierUpdater
) -> pl.DataFrame:
    """
    Extract gene results from LIRICAL results.
    Args:
        raw_result (pl.DataFrame): LIRICAL results dataframe.
        gene_identifier_updater (GeneIdentifierUpdater): GeneIdentifierUpdater object.
    Returns:
        pl.DataFrame: The extracted results.
    """
    return raw_result.select(
        [
            pl.col("entrezGeneId")
            .str.split(":")
            .list.get(1)
            .map_elements(
                gene_identifier_updater.obtain_gene_symbol_from_identifier, return_dtype=pl.String
            )
            .alias("gene_symbol"),
            pl.when(pl.col("compositeLR").cast(pl.Utf8) == "-∞")
            .then(float("-inf"))
            .otherwise(pl.col("compositeLR").cast(pl.Float64))
            .alias("score"),
        ]
    ).with_columns(
        [
            pl.col("gene_symbol")
            .map_elements(gene_identifier_updater.find_identifier, return_dtype=pl.String)
            .alias("gene_identifier"),
        ]
    )


def extract_variant_results(raw_result: pl.DataFrame) -> pl.DataFrame:
    """
    Extract variant results from LIRICAL results.
    Args:
        raw_result (pl.DataFrame): LIRICAL results dataframe.
    Returns:
        pl.DataFrame: The extracted results.
    """
    return (
        raw_result.with_columns([pl.col("variants").str.split("; ").alias("variant_list")])
        .explode("variant_list")
        .rename({"variant_list": "variant"})
        .drop("variants")
        .with_columns(
            [
                pl.col("variant").str.extract(r"^(\d+):", 1).alias("chrom"),
                pl.col("variant").str.extract(r":(\d+)", 1).cast(pl.Int64).alias("pos"),
                pl.col("variant").str.extract(r"([ACGT]+)>([ACGT]+)", 1).alias("ref"),
                pl.col("variant").str.extract(r"([ACGT]+)>([ACGT]+)", 2).alias("alt"),
                pl.when(pl.col("compositeLR").cast(pl.Utf8) == "-∞")
                .then(float("-inf"))
                .otherwise(pl.col("compositeLR").cast(pl.Float64))
                .alias("score"),
            ]
        )
        .select(
            [
                pl.col("chrom").cast(pl.String),
                pl.col("pos").alias("start").cast(pl.Int64),
                pl.struct("pos", "ref")
                .map_elements(lambda x: calculate_end_pos(x["pos"], x["ref"]))
                .alias("end")
                .cast(pl.String),
                pl.col("ref").cast(pl.String),
                pl.col("alt").cast(pl.String),
                pl.col("score").cast(pl.Float64),
            ]
        )
    )


def create_standardised_results(
    raw_results_dir: Path,
    output_dir: Path,
    phenopacket_dir: Path,
    sort_order: str,
    disease_analysis: bool,
    gene_analysis: bool,
    variant_analysis: bool,
) -> None:
    """Write standardised gene and variant results from LIRICAL tsv output."""
    gene_identifier_updater = GeneIdentifierUpdater(
        gene_identifier="ensembl_id",
        identifier_map=create_gene_identifier_map(),
    )
    sort_order = SortOrder.ASCENDING if sort_order.lower() == "ascending" else SortOrder.DESCENDING
    for result in files_with_suffix(raw_results_dir, ".tsv"):
        lirical_result = read_lirical_result(result)
        if gene_analysis:
            pheval_gene_result = extract_gene_results(lirical_result, gene_identifier_updater)
            generate_gene_result(
                results=pheval_gene_result,
                output_dir=output_dir,
                sort_order=sort_order,
                result_path=result,
                phenopacket_dir=phenopacket_dir,
            )
        if variant_analysis:
            pheval_variant_result = extract_variant_results(lirical_result)
            generate_variant_result(
                results=pheval_variant_result,
                output_dir=output_dir,
                sort_order=sort_order,
                result_path=result,
                phenopacket_dir=phenopacket_dir,
            )
        if disease_analysis:
            pheval_disease_result = extract_disease_results(lirical_result)
            generate_disease_result(
                results=pheval_disease_result,
                output_dir=output_dir,
                sort_order=sort_order,
                result_path=result,
                phenopacket_dir=phenopacket_dir,
            )
