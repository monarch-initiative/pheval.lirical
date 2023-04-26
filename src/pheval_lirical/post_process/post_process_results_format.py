import re
from pathlib import Path

import click
import pandas as pd
from pheval.post_processing.post_processing import (
    PhEvalGeneResult,
    PhEvalVariantResult,
    calculate_end_pos,
    generate_pheval_result,
)
from pheval.utils.file_utils import files_with_suffix
from pheval.utils.phenopacket_utils import (
    GeneIdentifierUpdater,
    create_gene_identifier_map,
    create_hgnc_dict,
)


def read_lirical_result(lirical_result_path: Path) -> pd.DataFrame:
    """Read LIRICAL tsv output and return a dataframe."""
    return pd.read_csv(lirical_result_path, delimiter="\t", comment="!")


class PhEvalGeneResultFromLirical:
    def __init__(
        self,
        lirical_result: pd.DataFrame,
        gene_identifier_updator: GeneIdentifierUpdater,
    ):
        self.lirical_result = lirical_result
        self.gene_identifier_updator = gene_identifier_updator

    @staticmethod
    def obtain_lirical_gene_identifier(result: pd.Series) -> str:
        """Obtain the NCBI gene id from LIRICAL result."""
        return result["entrezGeneId"]

    def obtain_gene_symbol(self, result: pd.Series) -> str:
        """Obtain the gene symbol from NCBI gene."""
        gene_identifier = self.obtain_lirical_gene_identifier(result)
        return self.gene_identifier_updator.obtain_gene_symbol_from_identifier(
            gene_identifier.split(":")[1]
        )

    def obtain_gene_identifier(self, result: pd.Series) -> str:
        """Obtain the ensembl gene identifier from gene symbol."""
        gene_symbol = self.obtain_gene_symbol(result)
        return self.gene_identifier_updator.find_identifier(gene_symbol)

    @staticmethod
    def obtain_score(result: pd.Series) -> float:
        """Obtain score from result."""
        return result["compositeLR"]

    def extract_pheval_requirements(self) -> [PhEvalGeneResult]:
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
    def obtain_variants(lirical_result_entry: pd.Series) -> str:
        """Obtain all variants listed in a result."""
        return lirical_result_entry["variants"]

    def split_variants(self, lirical_result_entry) -> [str]:
        """Split variants contained in a single result row."""
        variants = self.obtain_variants(lirical_result_entry)
        return [x.strip() for x in variants.split(";")]

    @staticmethod
    def get_variant_string(variant: str) -> str:
        """Obtain variant string."""
        return variant.split(" ")[0]

    def obtain_chromosome(self, variant_str: str) -> str:
        """Obtain chromosome from variant string."""
        variant = self.get_variant_string(variant_str)
        return variant.split(":")[0]

    def obtain_start(self, variant_str: str) -> int:
        """Obtain start position from variant string."""
        variant = self.get_variant_string(variant_str)
        return int(self.reg.search(variant.split(":")[1]).group("numbers"))

    def obtain_ref(self, variant_str: str) -> str:
        """Obtain reference allele from variant string."""
        variant = self.get_variant_string(variant_str)
        return self.reg.search(variant.split(":")[1]).group("rest").split(">")[0]

    def obtain_end(self, variant_str: str) -> int:
        """Obtain end position from variant string."""
        variant = self.get_variant_string(variant_str)
        return calculate_end_pos(self.obtain_start(variant), self.obtain_ref(variant))

    def obtain_alt(self, variant_str: str) -> str:
        """Obtain alternate allele from variant string."""
        variant = self.get_variant_string(variant_str)
        return self.reg.search(variant.split(":")[1]).group("rest").split(">")[1]

    def extract_pheval_requirements(self) -> PhEvalVariantResult:
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


def create_standardised_results(raw_results_dir: Path, output_dir: Path, sort_order: str) -> None:
    """Write standardised gene and variant results from LIRICAL tsv output."""
    identifier_map = create_gene_identifier_map()
    hgnc_data = create_hgnc_dict()
    gene_identifier_updator = GeneIdentifierUpdater(
        gene_identifier="ensembl_id", identifier_map=identifier_map, hgnc_data=hgnc_data
    )
    for result in files_with_suffix(raw_results_dir, ".tsv"):
        lirical_result = read_lirical_result(result)
        pheval_gene_result = PhEvalGeneResultFromLirical(
            lirical_result, gene_identifier_updator
        ).extract_pheval_requirements()
        generate_pheval_result(pheval_gene_result, sort_order, output_dir, result)
        pheval_variant_result = PhEvalVariantResultFromLirical(
            lirical_result
        ).extract_pheval_requirements()
        generate_pheval_result(pheval_variant_result, sort_order, output_dir, result)


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
def post_process(lirical_file: Path, output_dir: Path, sort_order: str):
    """Post-process LIRICAL .tsv results to PhEval gene and variant result format."""
    output_dir.joinpath("pheval_gene_results/").mkdir(exist_ok=True, parents=True)
    output_dir.joinpath("pheval_variant_results/").mkdir(exist_ok=True, parents=True)
    create_standardised_results(lirical_file, output_dir, sort_order)
