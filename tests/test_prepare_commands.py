import shutil
import tempfile
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

from pheval_lirical.prepare.prepare_manual_commands import (
    LiricalManualCommandLineArguments,
)
from pheval_lirical.prepare.prepare_commands import CommandCreator, CommandWriter

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


class TestCommandCreator(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.command_creator = CommandCreator(
            phenopacket_path=Path("/path/to/phenopacket.json"),
            phenopacket=phenopacket_with_excluded,
            lirical_jar=Path("/path/to/lirical.jar"),
            input_dir=Path("/path/to/lirical/data"),
            exomiser_data_dir=Path("/path/to/exomiser/data"),
            vcf_dir=Path("/path/to/vcf_dir"),
            results_dir=Path("/path/to/results_dir"),
        )
        cls.command_creator_none_excluded = CommandCreator(
            phenopacket_path=Path("/path/to/phenopacket.json"),
            phenopacket=phenopacket_without_excluded,
            lirical_jar=Path("/path/to/lirical.jar"),
            input_dir=Path("/path/to/lirical/data"),
            exomiser_data_dir=Path("/path/to/exomiser/data"),
            vcf_dir=Path("/path/to/vcf_dir"),
            results_dir=Path("/path/to/results_dir"),
        )

    def test_get_list_negated_phenotypic_features(self):
        self.assertEqual(
            self.command_creator.get_list_negated_phenotypic_features(), ["HP:0008494"]
        )

    def test_get_list_negated_phenotypic_features_none_excluded(self):
        self.assertEqual(
            self.command_creator_none_excluded.get_list_negated_phenotypic_features(), None
        )

    def test_get_list_observed_phenotypic_features(self):
        self.assertEqual(
            self.command_creator.get_list_observed_phenotypic_features(),
            ["HP:0000256", "HP:0002059", "HP:0100309", "HP:0003150", "HP:0001332"],
        )

    def test_get_vcf_path(self):
        self.assertEqual(self.command_creator.get_vcf_path(), "/path/to/vcf_dir/test_1.vcf")

    def test_get_vcf_assembly(self):
        self.assertEqual(self.command_creator.get_vcf_assembly(), "GRCh37")

    def test_add_manual_cli_arguments(self):
        self.assertEqual(
            self.command_creator.add_manual_cli_arguments(),
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

    def test_add_manual_cli_arguments_none_excluded(self):
        self.assertEqual(
            self.command_creator_none_excluded.add_manual_cli_arguments(),
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


class TestCommandWriter(unittest.TestCase):
    def setUp(self) -> None:
        self.test_dir = tempfile.mkdtemp()
        self.command_file_path = Path(self.test_dir).joinpath("test-commands.txt")
        self.command_writer = CommandWriter(self.command_file_path)
        self.command_arguments = LiricalManualCommandLineArguments(
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
        )

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_write_basic_manual_command(self):
        self.command_writer.write_java_command(self.command_arguments)
        self.command_writer.file.close()
        with open(self.command_file_path) as f:
            content = f.readlines()
        f.close()
        self.assertEqual(content, ["java -jar /path/to/lirical.jar R "])

    def test_write_observed_phenotypic_features(self):
        self.command_writer.write_observed_phenotypic_features(self.command_arguments)
        self.command_writer.file.close()
        with open(self.command_file_path) as f:
            content = f.readlines()
        f.close()
        self.assertEqual(
            content,
            ["--observed-phenotypes HP:0000256,HP:0002059,HP:0100309,HP:0003150,HP:0001332"],
        )

    def test_write_negated_phenotypic_features(self):
        self.command_writer.write_negated_phenotypic_features(self.command_arguments)
        self.command_writer.file.close()
        with open(self.command_file_path) as f:
            content = f.readlines()
        f.close()
        self.assertEqual(content, [" --negated-phenotypes HP:0008494"])

    def test_write_vcf(self):
        self.command_writer.write_vcf_file_properties(self.command_arguments)
        self.command_writer.file.close()
        with open(self.command_file_path) as f:
            content = f.readlines()
        f.close()
        self.assertEqual(
            content,
            [' --vcf /path/to/vcf_dir/test_1.vcf --assembly GRCh37 --sample-id "test-subject-1"'],
        )

    def test_write_data_dirs(self):
        self.command_writer.write_lirical_data_dir(self.command_arguments)
        self.command_writer.file.close()
        with open(self.command_file_path) as f:
            content = f.readlines()
        f.close()
        self.assertEqual(
            content, [" --data /path/to/lirical/data --exomiser /path/to/exomiser/data"]
        )

    def test_write_output_parameters(self):
        self.command_writer.write_output_parameters(self.command_arguments)
        self.command_writer.file.close()
        with open(self.command_file_path) as f:
            content = f.readlines()
        f.close()
        self.assertEqual(
            content,
            [" --prefix phenopacket --output-directory /path/to/results_dir --output-format tsv"],
        )

    def test_write_manual_command(self):
        self.command_writer.write_manual_command(self.command_arguments)
        self.command_writer.file.close()
        with open(self.command_file_path) as f:
            content = f.readlines()
        f.close()
        self.assertEqual(
            content,
            [
                "java -jar /path/to/lirical.jar R --observed-phenotypes "
                "HP:0000256,HP:0002059,HP:0100309,HP:0003150,HP:0001332 --negated-phenotypes HP:0008494 "
                '--vcf /path/to/vcf_dir/test_1.vcf --assembly GRCh37 --sample-id "test-subject-1" '
                "--data /path/to/lirical/data --exomiser /path/to/exomiser/data --prefix phenopacket "
                "--output-directory /path/to/results_dir --output-format tsv\n"
            ],
        )

    def test_close(self):
        self.command_writer.close()
        self.assertTrue(self.command_writer.file.closed)
