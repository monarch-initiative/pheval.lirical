from dataclasses import dataclass
from pathlib import Path
from phenopackets import PhenotypicFeature, Phenopacket
from pheval.utils.phenopacket_utils import phenopacket_reader, PhenopacketUtil


@dataclass
class LiricalCommandLineArguments:
    observed_phenotypes: [PhenotypicFeature]
    negated_phenotypes: [PhenotypicFeature] or None
    assembly: str
    vcf_file_path: Path
    lirical_data: Path
    exomiser_data: Path or None


def obtain_negated_phenotypes(phenopacket: Phenopacket):
    negated_phenotypic_features = []
    phenotypes = PhenopacketUtil(phenopacket).phenotypic_features()
    for phenotype in phenotypes:
        if phenotype.excluded:
            negated_phenotypic_features.append(phenotype)
    return negated_phenotypic_features


