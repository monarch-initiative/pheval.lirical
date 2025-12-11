import unittest

import polars as pl
from pheval.utils.phenopacket_utils import (
    GeneIdentifierUpdater,
    create_gene_identifier_map,
)

from src.pheval_lirical.post_process.post_process_results_format import (
    extract_disease_results,
    extract_gene_results,
    extract_variant_results,
)

lirical_results = pl.DataFrame(
    [
        {
            "rank": 1,
            "diseaseName": "Glutaric acidemia I",
            "diseaseCurie": "OMIM:231670",
            "pretestprob": 1 / 8371,
            "posttestprob": "65.60%",
            "compositeLR": "4.203",
            "entrezGeneId": "NCBIGene:2639",
            "variants": "19:12998205G>C NM_006563.3:: pathogenicity:0.0 [0/0]; "
            "19:13002033A>G NM_000159.3:c.-35+17A>G:p.(=) pathogenicity:0.0 [1/1]; "
            "19:13002400G>C NM_000159.3:c.127+64G>C:p.(=) pathogenicity:0.0 [0/1]; "
            "19:13002563T>G NM_000159.3:c.128-82T>G:p.(=) pathogenicity:0.0 [0/1]; "
            "19:13007113G>A NM_000159.3:c.730G>A:p.(G244S) pathogenicity:1.0 [1/1]; "
            "19:13008264C>T NM_000159.3:c.1082+22C>T:p.(=) pathogenicity:0.0 [0/1]; "
            "19:13008607G>T NM_000159.3:c.1173G>T:p.(=) pathogenicity:0.0 [0/1]; "
            "19:13010520A>G NM_013976.3:c.1250A>G:p.(Q417R) pathogenicity:0.0 [0/1]; "
            "19:13010643G>T NM_001105578.1:c.612+175C>A:p.(=) pathogenicity:0.0 [0/1]",
        },
        {
            "rank": 2,
            "diseaseName": "Diencephalic-mesencephalic junction dysplasia syndrome 2",
            "diseaseCurie": "OMIM:618646",
            "pretestprob": 1 / 8371,
            "posttestprob": "0.00%",
            "compositeLR": "-1.439",
            "entrezGeneId": "NCBIGene:170825",
            "variants": "4:54950725C>T NM_133267.2:: pathogenicity:0.0 [1/1]; "
            "4:54966667C>T NM_133267.2:c.156C>T:p.(=) pathogenicity:0.0 [1/1]; "
            "4:54966830G>A NM_133267.2:c.319A>A:p.(=) pathogenicity:0.0 [1/1]; "
            "4:54966910G>GCACCAC NM_133267.2:c.402_407dup:p.(H138_H139dup) pathogenicity:0.9 [1/1]; "
            "4:54967709C>A NM_133267.2:c.575-40C>A:p.(=) pathogenicity:0.0 [1/1]; "
            "4:55002331T>C NM_133267.2:: pathogenicity:0.0 [1/1]; "
            "4:55009968C>T NM_133267.2:: pathogenicity:0.0 [1/1]; "
            "4:55010012G>T NM_133267.2:: pathogenicity:0.0 [1/1]; "
            "4:55011769T>C NM_133267.2:: pathogenicity:0.0 [1/1]; "
            "4:55013951T>C NM_133267.2:: pathogenicity:0.0 [1/1]; "
            "4:55026539ACT>A NM_133267.2:: pathogenicity:0.0 [./.]",
        },
    ]
)


class TestPostProcessResultsFormat(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.gene_identifier_updater = GeneIdentifierUpdater(
            "ensembl_id", create_gene_identifier_map()
        )

    def test_extract_disease_results(self):
        self.assertTrue(
            extract_disease_results(lirical_results).equals(
                pl.DataFrame(
                    [
                        {"disease_identifier": "OMIM:231670", "score": 4.203},
                        {"disease_identifier": "OMIM:618646", "score": -1.439},
                    ]
                )
            )
        )

    def test_extract_gene_results(self):
        self.assertTrue(
            extract_gene_results(lirical_results, self.gene_identifier_updater).equals(
                pl.DataFrame(
                    [
                        {
                            "gene_symbol": "GCDH",
                            "score": 4.203,
                            "gene_identifier": "ENSG00000105607",
                        },
                        {
                            "gene_symbol": "GSX2",
                            "score": -1.439,
                            "gene_identifier": "ENSG00000180613",
                        },
                    ]
                )
            )
        )

    def test_extract_variant_results(self):
        print(extract_variant_results(lirical_results))


#
#
# class TestPhEvalVariantResultFromLirical(unittest.TestCase):
#     @classmethod
#     def setUpClass(cls) -> None:
#         cls.variant_result = PhEvalVariantResultFromLirical(lirical_results)
#
#     def test_obtain_score(self):
#         self.assertEqual(self.variant_result.obtain_score(lirical_result), 4.203)
#
#     def test_obtain_variants(self):
#         self.assertEqual(
#             self.variant_result.obtain_variants(lirical_result),
#             "19:12998205G>C NM_006563.3:: pathogenicity:0.0 [0/0]; "
#             "19:13002033A>G NM_000159.3:c.-35+17A>G:p.(=) pathogenicity:0.0 [1/1]; "
#             "19:13002400G>C NM_000159.3:c.127+64G>C:p.(=) pathogenicity:0.0 [0/1]; "
#             "19:13002563T>G NM_000159.3:c.128-82T>G:p.(=) pathogenicity:0.0 [0/1]; "
#             "19:13007113G>A NM_000159.3:c.730G>A:p.(G244S) pathogenicity:1.0 [1/1]; "
#             "19:13008264C>T NM_000159.3:c.1082+22C>T:p.(=) pathogenicity:0.0 [0/1]; "
#             "19:13008607G>T NM_000159.3:c.1173G>T:p.(=) pathogenicity:0.0 [0/1]; "
#             "19:13010520A>G NM_013976.3:c.1250A>G:p.(Q417R) pathogenicity:0.0 [0/1]; "
#             "19:13010643G>T NM_001105578.1:c.612+175C>A:p.(=) pathogenicity:0.0 [0/1]",
#         )
#
#     def test_split_variants(self):
#         self.assertEqual(
#             self.variant_result.split_variants(lirical_result),
#             [
#                 "19:12998205G>C NM_006563.3:: pathogenicity:0.0 [0/0]",
#                 "19:13002033A>G NM_000159.3:c.-35+17A>G:p.(=) pathogenicity:0.0 [1/1]",
#                 "19:13002400G>C NM_000159.3:c.127+64G>C:p.(=) pathogenicity:0.0 [0/1]",
#                 "19:13002563T>G NM_000159.3:c.128-82T>G:p.(=) pathogenicity:0.0 [0/1]",
#                 "19:13007113G>A NM_000159.3:c.730G>A:p.(G244S) pathogenicity:1.0 [1/1]",
#                 "19:13008264C>T NM_000159.3:c.1082+22C>T:p.(=) pathogenicity:0.0 [0/1]",
#                 "19:13008607G>T NM_000159.3:c.1173G>T:p.(=) pathogenicity:0.0 [0/1]",
#                 "19:13010520A>G NM_013976.3:c.1250A>G:p.(Q417R) pathogenicity:0.0 [0/1]",
#                 "19:13010643G>T NM_001105578.1:c.612+175C>A:p.(=) pathogenicity:0.0 [0/1]",
#             ],
#         )
#
#     def test_get_variant_string(self):
#         self.assertEqual(
#             self.variant_result.get_variant_string(
#                 "19:12998205G>C NM_006563.3:: pathogenicity:0.0 [0/0]"
#             ),
#             "19:12998205G>C",
#         )
#
#     def test_obtain_chromosome(self):
#         self.assertEqual(self.variant_result.obtain_chromosome("19:12998205G>C"), "19")
#
#     def test_obtain_start(self):
#         self.assertEqual(self.variant_result.obtain_start("19:12998205G>C"), 12998205)
#
#     def test_obtain_ref(self):
#         self.assertEqual(self.variant_result.obtain_ref("19:12998205G>C"), "G")
#
#     def test_obtain_end(self):
#         self.assertEqual(self.variant_result.obtain_end("19:12998205G>C"), 12998205)
#
#     def test_obtain_end_add_pos(self):
#         self.assertEqual(self.variant_result.obtain_end("19:12998205GGGGGGC>C"), 12998211)
#
#     def test_obtain_alt(self):
#         self.assertEqual(self.variant_result.obtain_alt("19:12998205G>C"), "C")
#
#     def test_extract_pheval_requirements(self):
#         self.assertEqual(
#             self.variant_result.extract_pheval_requirements(),
#             [
#                 PhEvalVariantResult(
#                     chromosome="19", start=12998205, end=12998205, ref="G", alt="C", score=4.203
#                 ),
#                 PhEvalVariantResult(
#                     chromosome="19", start=13002033, end=13002033, ref="A", alt="G", score=4.203
#                 ),
#                 PhEvalVariantResult(
#                     chromosome="19", start=13002400, end=13002400, ref="G", alt="C", score=4.203
#                 ),
#                 PhEvalVariantResult(
#                     chromosome="19", start=13002563, end=13002563, ref="T", alt="G", score=4.203
#                 ),
#                 PhEvalVariantResult(
#                     chromosome="19", start=13007113, end=13007113, ref="G", alt="A", score=4.203
#                 ),
#                 PhEvalVariantResult(
#                     chromosome="19", start=13008264, end=13008264, ref="C", alt="T", score=4.203
#                 ),
#                 PhEvalVariantResult(
#                     chromosome="19", start=13008607, end=13008607, ref="G", alt="T", score=4.203
#                 ),
#                 PhEvalVariantResult(
#                     chromosome="19", start=13010520, end=13010520, ref="A", alt="G", score=4.203
#                 ),
#                 PhEvalVariantResult(
#                     chromosome="19", start=13010643, end=13010643, ref="G", alt="T", score=4.203
#                 ),
#                 PhEvalVariantResult(
#                     chromosome="4", start=54950725, end=54950725, ref="C", alt="T", score=-1.439
#                 ),
#                 PhEvalVariantResult(
#                     chromosome="4", start=54966667, end=54966667, ref="C", alt="T", score=-1.439
#                 ),
#                 PhEvalVariantResult(
#                     chromosome="4", start=54966830, end=54966830, ref="G", alt="A", score=-1.439
#                 ),
#                 PhEvalVariantResult(
#                     chromosome="4",
#                     start=54966910,
#                     end=54966910,
#                     ref="G",
#                     alt="GCACCAC",
#                     score=-1.439,
#                 ),
#                 PhEvalVariantResult(
#                     chromosome="4", start=54967709, end=54967709, ref="C", alt="A", score=-1.439
#                 ),
#                 PhEvalVariantResult(
#                     chromosome="4", start=55002331, end=55002331, ref="T", alt="C", score=-1.439
#                 ),
#                 PhEvalVariantResult(
#                     chromosome="4", start=55009968, end=55009968, ref="C", alt="T", score=-1.439
#                 ),
#                 PhEvalVariantResult(
#                     chromosome="4", start=55010012, end=55010012, ref="G", alt="T", score=-1.439
#                 ),
#                 PhEvalVariantResult(
#                     chromosome="4", start=55011769, end=55011769, ref="T", alt="C", score=-1.439
#                 ),
#                 PhEvalVariantResult(
#                     chromosome="4", start=55013951, end=55013951, ref="T", alt="C", score=-1.439
#                 ),
#                 PhEvalVariantResult(
#                     chromosome="4", start=55026539, end=55026541, ref="ACT", alt="A", score=-1.439
#                 ),
#             ],
#         )
#
#
