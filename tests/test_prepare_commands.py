import unittest
from pathlib import Path

from phenopackets import (
    Diagnosis,
    File,
    GeneDescriptor,
    GenomicInterpretation,
    Individual,
    Interpretation,
    MetaData,
    OntologyClass,
    Phenopacket,
    PhenotypicFeature,
    Resource,
    VariantInterpretation,
    VariationDescriptor,
    VcfRecord,
)

from pheval_lirical.prepare.prepare_commands import (
    LiricalManualCommandLineArguments,
    create_command_line_arguments,
)

interpretations = [
    Interpretation(
        id="test-subject-1-int",
        progress_status="SOLVED",
        diagnosis=Diagnosis(
            genomic_interpretations=[
                GenomicInterpretation(
                    subject_or_biosample_id="test-subject-1",
                    interpretation_status=4,
                    variant_interpretation=VariantInterpretation(
                        acmg_pathogenicity_classification="NOT_PROVIDED",
                        therapeutic_actionability="UNKNOWN_ACTIONABILITY",
                        variation_descriptor=VariationDescriptor(
                            gene_context=GeneDescriptor(
                                value_id="ENSG00000102302",
                                symbol="FGD1",
                                alternate_ids=[
                                    "HGNC:3663",
                                    "ncbigene:2245",
                                    "ensembl:ENSG00000102302",
                                    "symbol:FGD1",
                                ],
                            ),
                            vcf_record=VcfRecord(
                                genome_assembly="GRCh37",
                                chrom="X",
                                pos=54492285,
                                ref="C",
                                alt="T",
                            ),
                            allelic_state=OntologyClass(
                                id="GENO:0000134",
                                label="hemizygous",
                            ),
                        ),
                    ),
                ),
                GenomicInterpretation(
                    subject_or_biosample_id="test-subject-1",
                    interpretation_status=4,
                    variant_interpretation=VariantInterpretation(
                        acmg_pathogenicity_classification="NOT_PROVIDED",
                        therapeutic_actionability="UNKNOWN_ACTIONABILITY",
                        variation_descriptor=VariationDescriptor(
                            gene_context=GeneDescriptor(
                                value_id="ENSG00000176225",
                                symbol="RTTN",
                                alternate_ids=[
                                    "HGNC:18654",
                                    "ncbigene:25914",
                                    "ensembl:ENSG00000176225",
                                    "symbol:RTTN",
                                ],
                            ),
                            vcf_record=VcfRecord(
                                genome_assembly="GRCh37",
                                chrom="18",
                                pos=67691994,
                                ref="G",
                                alt="A",
                            ),
                            allelic_state=OntologyClass(
                                id="GENO:0000402", label="compound heterozygous"
                            ),
                        ),
                    ),
                ),
            ]
        ),
    )
]

phenotypic_features_none_excluded = [
    PhenotypicFeature(type=OntologyClass(id="HP:0000256", label="Macrocephaly")),
    PhenotypicFeature(type=OntologyClass(id="HP:0002059", label="Cerebral atrophy")),
    PhenotypicFeature(type=OntologyClass(id="HP:0100309", label="Subdural hemorrhage")),
    PhenotypicFeature(type=OntologyClass(id="HP:0003150", label="Glutaric aciduria")),
    PhenotypicFeature(type=OntologyClass(id="HP:0001332", label="Dystonia")),
]
phenotypic_features_with_excluded = [
    PhenotypicFeature(type=OntologyClass(id="HP:0000256", label="Macrocephaly")),
    PhenotypicFeature(type=OntologyClass(id="HP:0002059", label="Cerebral atrophy")),
    PhenotypicFeature(type=OntologyClass(id="HP:0100309", label="Subdural hemorrhage")),
    PhenotypicFeature(type=OntologyClass(id="HP:0003150", label="Glutaric aciduria")),
    PhenotypicFeature(type=OntologyClass(id="HP:0001332", label="Dystonia")),
    PhenotypicFeature(
        type=OntologyClass(id="HP:0008494", label="Inferior lens subluxation"), excluded=True
    ),
]
phenopacket_files = [
    File(
        uri="test/path/to/test_1.vcf",
        file_attributes={"fileFormat": "VCF", "genomeAssembly": "GRCh37"},
    ),
    File(
        uri="test_1.ped",
        file_attributes={"fileFormat": "PED", "genomeAssembly": "GRCh37"},
    ),
]
phenopacket_metadata = MetaData(
    created_by="pheval-converter",
    resources=[
        Resource(
            id="hp",
            name="human phenotype ontology",
            url="http://purl.obolibrary.org/obo/hp.owl",
            version="hp/releases/2019-11-08",
            namespace_prefix="HP",
            iri_prefix="http://purl.obolibrary.org/obo/HP_",
        )
    ],
    phenopacket_schema_version="2.0",
)
phenopacket_with_excluded = Phenopacket(
    id="test-subject",
    subject=Individual(id="test-subject-1", sex=1),
    phenotypic_features=phenotypic_features_with_excluded,
    interpretations=interpretations,
    files=phenopacket_files,
    meta_data=phenopacket_metadata,
)

phenopacket_without_excluded = Phenopacket(
    id="test-subject",
    subject=Individual(id="test-subject-1", sex=1),
    phenotypic_features=phenotypic_features_none_excluded,
    interpretations=interpretations,
    files=phenopacket_files,
    meta_data=phenopacket_metadata,
)


class TestCreateCommandLineArguments(unittest.TestCase):
    def test_create_command_line_arguments_none_excluded(self):
        self.assertEqual(
            create_command_line_arguments(
                phenopacket=phenopacket_without_excluded,
                lirical_jar=Path("/path/to/lirical.jar"),
                input_dir=Path("/path/to/lirical/data"),
                exomiser_data_dir=Path("/path/to/exomiser/data"),
                phenopacket_path=Path("/path/to/phenopacket.json"),
                vcf_dir=Path("/path/to/vcf_dir"),
                results_dir=Path("/path/to/results_dir"),
            ),
            LiricalManualCommandLineArguments(
                lirical_jar_file=Path("/path/to/lirical.jar"),
                observed_phenotypes=[
                    "HP:0000256",
                    "HP:0002059",
                    "HP:0100309",
                    "HP:0003150",
                    "HP:0001332",
                ],
                negated_phenotypes=None,
                assembly="GRCh37",
                vcf_file_path="/path/to/vcf_dir/test_1.vcf",
                sample_id="test-subject-1",
                lirical_data=Path("/path/to/lirical/data"),
                exomiser_data=Path("/path/to/exomiser/data"),
                output_dir=Path("/path/to/results_dir"),
                output_prefix="phenopacket",
            ),
        )

    def test_create_command_line_arguments_excluded(self):
        self.assertEqual(
            create_command_line_arguments(
                phenopacket=phenopacket_with_excluded,
                lirical_jar=Path("/path/to/lirical.jar"),
                input_dir=Path("/path/to/lirical/data"),
                exomiser_data_dir=Path("/path/to/exomiser/data"),
                phenopacket_path=Path("/path/to/phenopacket.json"),
                vcf_dir=Path("/path/to/vcf_dir"),
                results_dir=Path("/path/to/results_dir"),
            ),
            LiricalManualCommandLineArguments(
                lirical_jar_file=Path("/path/to/lirical.jar"),
                observed_phenotypes=[
                    "HP:0000256",
                    "HP:0002059",
                    "HP:0100309",
                    "HP:0003150",
                    "HP:0001332",
                ],
                negated_phenotypes=["HP:0008494"],
                assembly="GRCh37",
                vcf_file_path="/path/to/vcf_dir/test_1.vcf",
                sample_id="test-subject-1",
                lirical_data=Path("/path/to/lirical/data"),
                exomiser_data=Path("/path/to/exomiser/data"),
                output_dir=Path("/path/to/results_dir"),
                output_prefix="phenopacket",
            ),
        )
