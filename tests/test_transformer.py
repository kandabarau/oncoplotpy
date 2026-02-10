"""OncoPlotPy: Genomic waterfall plots in Matplotlib."""
import unittest
import pandas as pd
import io
import numpy as np
from oncoplotpy.transformer import OncoMatrix, SAMPLE_FIELD, GENE_FIELD, VC_FIELD, VC_MULTIHIT_LABEL


class TestOncoMatrix(unittest.TestCase):
    def setUp(self):
        # The "Classic 12" MAF Dataset
        self.toy_maf = (
            f"{SAMPLE_FIELD}\t{GENE_FIELD}\t{VC_FIELD}\tgenes\tsamples\tphenotype\n"
            "SAMPLE-001\t\t\t\t12\tT\n"
            "SAMPLE-002\tDNMT3A\tMissense_Mutation\t4\t7\tT\n"
            "SAMPLE-003\tASXL1\tIn_Frame_Del\t2\t6\tT\n"
            "SAMPLE-003\tASXL1\tMissense_Mutation\t2\t6\tT\n"
            "SAMPLE-004\tASXL1\tMissense_Mutation\t2\t2\tT\n"
            "SAMPLE-004\tRUNX1\tNonsense_Mutation\t3\t2\tT\n"
            "SAMPLE-005\tTP53\tIntron\t1\t8\tT\n"
            "SAMPLE-006\t\t\t\t9\tT\n"
            "SAMPLE-007\tASXL1\tMissense_Mutation\t2\t10\tC\n"
            "SAMPLE-008\tASXL1\tMissense_Mutation\t2\t1\tC\n"
            "SAMPLE-008\tGATA2\tIntron\t5\t1\tC\n"
            "SAMPLE-009\tTP53\tMissense_Mutation\t1\t3\tC\n"
            "SAMPLE-010\tRUNX1\tFrame_Shift_Del\t3\t4\tC\n"
            "SAMPLE-011\t\t\t\t5\tC\n"
            "SAMPLE-012\t\t\t\t11\tT"
        )

    def get_maf_buffer(self):
        return io.StringIO(self.toy_maf)

    def test_logic_with_classic_maf(self):
        """Test basic mutation logic: Multi-hit and Intron filtering."""
        onco = OncoMatrix(maf_path=self.get_maf_buffer())
        onco.prepare()

        # Verify SAMPLE-003 ASXL1 is Multi_Hit
        self.assertEqual(onco.data.loc['ASXL1', 'SAMPLE-003'], VC_MULTIHIT_LABEL)
        # Verify SAMPLE-005 TP53 is NaN (Intron filtered)
        self.assertTrue(pd.isna(onco.data.loc['TP53', 'SAMPLE-005']))

    def test_external_list_ordering(self):
        """Test sorting using external Python lists (CLI/Jupyter style)."""
        onco = OncoMatrix(maf_path=self.get_maf_buffer())
        onco.prepare()

        # Define external lists (decided outside the class)
        forced_gene_order = ['DNMT3A', 'RUNX1', 'ASXL1', 'TP53']
        forced_sample_order = ['SAMPLE-010', 'SAMPLE-002', 'SAMPLE-008']

        # Apply the lists directly
        onco.sort_genes(genes_order=forced_gene_order)
        onco.sort_samples(samples_order=forced_sample_order)

        # Verify Genes match the list order exactly
        self.assertEqual(onco.genes, forced_gene_order)

        # Verify Samples match the list order exactly
        self.assertEqual(onco.samples, forced_sample_order)

        # Verify data integrity: SAMPLE-010 should have RUNX1 mutation
        self.assertEqual(onco.data.loc['RUNX1', 'SAMPLE-010'], 'Frame_Shift_Del')

    def test_phenotype_sorting_waterfall(self):
        """Test grouping by phenotype (Largest Group First) then staircase."""
        onco = OncoMatrix(
            maf_path=self.get_maf_buffer(),
            annotation_fields=['phenotype']
        )
        onco.prepare().sort_genes().sort_samples(sortby=['phenotype'])

        # Phenotype T (n=7) comes before C (n=5)
        self.assertEqual(onco.annotation.loc[onco.samples[0], 'phenotype'], 'T')

        # Within Phenotype C, waterfall logic: 
        # SAMPLE-008/007 (ASXL1) should be before SAMPLE-010 (RUNX1)
        c_samples = [s for s in onco.samples if onco.annotation.loc[s, 'phenotype'] == 'C']
        self.assertLess(c_samples.index('SAMPLE-008'), c_samples.index('SAMPLE-010'))

    def test_annotation_persistence_without_sorting(self):
        """Ensure phenotype is available for plots and matches matrix order even in waterfall sort."""
        onco = OncoMatrix(
            maf_path=self.get_maf_buffer(),
            annotation_fields=['phenotype']
        )
        # Perform a pure waterfall sort (no 'sortby' provided)
        onco.prepare().sort_genes().sort_samples()

        # 1. Check presence
        self.assertIn('phenotype', onco.annotation.columns)

        # 2. Check length
        self.assertEqual(len(onco.annotation), len(onco.samples))

        # 3. Check matching order (Critical Sync Check)
        # The index of the annotation must be identical to the list of samples in the matrix
        self.assertEqual(onco.annotation.index.tolist(), onco.samples)

        # 4. Verify specific mapping
        # In waterfall, SAMPLE-004 (2 mutations) should be early.
        # Let's ensure its phenotype 'T' is still correctly mapped.
        self.assertEqual(onco.annotation.loc['SAMPLE-004', 'phenotype'], 'T')


if __name__ == '__main__':
    unittest.main()
