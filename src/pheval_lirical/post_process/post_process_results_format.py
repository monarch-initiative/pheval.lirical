import re
from pathlib import Path

import click
import pandas as pd
from pheval.post_processing.post_processing import (
    PhEvalGeneResult,
    PhEvalVariantResult,
    RankedPhEvalGeneResult,
    RankedPhEvalVariantResult,
    create_pheval_result,
    write_pheval_gene_result,
    write_pheval_variant_result,
)
from pheval.utils.file_utils import files_with_suffix
from pheval.utils.phenopacket_utils import GeneIdentifierUpdater, create_hgnc_dict


def read_lirical_result(lirical_result_path: Path) -> pd.DataFrame:
    """Read LIRICAL tsv output and return a dataframe."""
    return pd.read_csv(lirical_result_path, delimiter="\t", comment="!")


def calculate_end(variant_start: int, variant_ref: str) -> int:
    """Calculate the end position for a variant."""
    # TODO move to PhEval
    return int(variant_start) + int(len(variant_ref)) - 1


def obtain_gene_symbol_from_identifier(
    gene_identifier: str, hgnc_data: dict, identifier: str
) -> str:
    """Obtain gene symbol from a gene identifier."""
    # TODO move to PhEval GeneIdentifierupdator
    for symbol, data in hgnc_data.items():
        if gene_identifier == data[identifier]:
            return symbol


class PhEvalGeneResultFromLirical:
    def __init__(
        self,
        lirical_result: pd.DataFrame,
        hgnc_data: dict,
        gene_identifier_updator: GeneIdentifierUpdater,
    ):
        self.lirical_result = lirical_result
        self.hgnc_data = hgnc_data
        self.gene_identifier_updator = gene_identifier_updator

    @staticmethod
    def obtain_lirical_gene_identifier(result: pd.Series) -> str:
        """Obtain the NCBI gene id from LIRICAL result."""
        return result["entrezGeneId"]

    def obtain_gene_symbol(self, result: pd.Series):
        """Obtain the gene symbol from NCBI gene."""
        gene_identifier = self.obtain_lirical_gene_identifier(result)
        return obtain_gene_symbol_from_identifier(
            gene_identifier.split(":")[1], self.hgnc_data, "entrez_id"
        )

    def obtain_gene_identifier(self, result: pd.Series):
        """Obtain the ensembl gene identifier from gene symbol."""
        gene_symbol = self.obtain_gene_symbol(result)
        return self.gene_identifier_updator.find_identifier(gene_symbol)

    @staticmethod
    def obtain_score(result: pd.Series) -> float:
        """Obtain score from result."""
        return result["compositeLR"]

    def extract_pheval_requirements(self):
        """Extract data required to produce PhEval gene output."""
        simplified_gene_results = []
        for _index, result in self.lirical_result.iterrows():
            simplified_gene_results.append(
                PhEvalGeneResult(
                    gene_symbol=self.obtain_gene_symbol(result),
                    gene_identifier=self.obtain_gene_identifier(result),
                    score=self.obtain_score(result),
                )
            )
        return simplified_gene_results


class PhEvalVariantResultFromLirical:
    def __init__(self, lirical_result: pd.DataFrame):
        self.lirical_result = lirical_result
        self.reg = re.compile(r"(?P<numbers>\d*)(?P<rest>.*)")

    @staticmethod
    def obtain_score(lirical_result_entry: pd.Series) -> float:
        """Obtain score from result."""
        return lirical_result_entry["compositeLR"]

    @staticmethod
    def obtain_variants(lirical_result_entry: pd.Series):
        """Obtain all variants listed in a result."""
        return lirical_result_entry["variants"]

    def split_variants(self, lirical_result_entry):
        """Split variants contained in a single result row."""
        variants = self.obtain_variants(lirical_result_entry)
        return [x.strip() for x in variants.split(";")]

    @staticmethod
    def get_variant_string(variant: str):
        """Obtain variant string."""
        return variant.split(" ")[0]

    def obtain_chromosome(self, variant_str: str):
        """Obtain chromosome from variant string."""
        variant = self.get_variant_string(variant_str)
        return variant.split(":")[0]

    def obtain_start(self, variant_str: str):
        """Obtain start position from variant string."""
        variant = self.get_variant_string(variant_str)
        return int(self.reg.search(variant.split(":")[1]).group("numbers"))

    def obtain_ref(self, variant_str: str):
        """Obtain reference allele from variant string."""
        variant = self.get_variant_string(variant_str)
        return self.reg.search(variant.split(":")[1]).group("rest").split(">")[0]

    def obtain_end(self, variant_str: str):
        """Obtain end position from variant string."""
        variant = self.get_variant_string(variant_str)
        return calculate_end(self.obtain_start(variant), self.obtain_ref(variant))

    def obtain_alt(self, variant_str: str):
        """Obtain alternate allele from variant string."""
        variant = self.get_variant_string(variant_str)
        return self.reg.search(variant.split(":")[1]).group("rest").split(">")[1]

    def extract_pheval_requirements(self):
        """Extract data required to produce PhEval variant output."""
        simplified_variant_results = []
        for _index, result in self.lirical_result.iterrows():
            variant_results = self.split_variants(result)
            for variant_result in variant_results:
                if variant_result != "":
                    variant = self.get_variant_string(variant_result)
                    simplified_variant_results.append(
                        PhEvalVariantResult(
                            chromosome=self.obtain_chromosome(variant),
                            start=self.obtain_start(variant),
                            end=self.obtain_end(variant),
                            ref=self.obtain_ref(variant),
                            alt=self.obtain_alt(variant),
                            score=self.obtain_score(result),
                        )
                    )
        return simplified_variant_results


def create_variant_gene_result_from_lirical(
    lirical_tsv_result: pd.DataFrame,
    sort_order: str,
) -> [RankedPhEvalVariantResult]:
    """Create ranked PhEval variant result from LIRICAL tsv."""
    pheval_variant_result = PhEvalVariantResultFromLirical(
        lirical_tsv_result
    ).extract_pheval_requirements()
    return create_pheval_result(pheval_variant_result, sort_order)


def create_pheval_gene_result_from_lirical(
    lirical_tsv_result: pd.DataFrame,
    gene_identifier_updator: GeneIdentifierUpdater,
    hgnc_data: dict,
    sort_order: str,
) -> [RankedPhEvalGeneResult]:
    """Create ranked PhEval gene result from LIRICAL tsv."""
    pheval_gene_result = PhEvalGeneResultFromLirical(
        lirical_tsv_result, hgnc_data, gene_identifier_updator
    ).extract_pheval_requirements()
    return create_pheval_result(pheval_gene_result, sort_order)


def create_standardised_results(results_dir: Path, output_dir: Path, sort_order: str) -> None:
    """Write standardised gene and variant results from LIRICAL tsv output."""
    output_dir.joinpath("pheval_gene_results/").mkdir(exist_ok=True, parents=True)
    output_dir.joinpath("pheval_variant_results/").mkdir(exist_ok=True, parents=True)
    hgnc_data = create_hgnc_dict()
    gene_identifier_updator = GeneIdentifierUpdater(hgnc_data, gene_identifier="ensembl_id")
    for result in files_with_suffix(results_dir, ".tsv"):
        lirical_result = read_lirical_result(result)
        pheval_gene_result = create_pheval_gene_result_from_lirical(
            lirical_result, gene_identifier_updator, hgnc_data, sort_order
        )
        write_pheval_gene_result(pheval_gene_result, output_dir, result)
        pheval_variant_result = create_variant_gene_result_from_lirical(lirical_result, sort_order)
        write_pheval_variant_result(pheval_variant_result, output_dir, result)


@click.command("post-process")
@click.option(
    "--lirical-file", "-f", required=True, help="Path to Lirical results file.", type=Path
)
@click.option("--output-dir", "-o", required=True, help="Path to output directory.", type=Path)
@click.option(
    "--sort-order",
    "-s",
    required=True,
    help="sort order.",
    type=click.Choice(["ascending", "descending"]),
    default="descending",
    show_default=True,
)
def post_process(lirical_file: Path, output_dir, sort_order):
    """Post-process LIRICAL .tsv results to PhEval gene and variant result format."""
    create_standardised_results(lirical_file, output_dir, sort_order)
