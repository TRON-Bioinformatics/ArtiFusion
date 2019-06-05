#!/usr/bin/env python

import sys
import os
import unittest

sys.path.append(os.path.join(os.path.dirname( __file__ ), '..'))

from artifuse import ArtiFusion

class TestArtiFusion(unittest.TestCase):
    def setUp(self):
        input_table = "Arti_Fusions_Candidates.csv"
        working_dir = "test"
        reference_genome = "test.fa"
        bed_file = "test.bed"
        gene_symbols_file = "test_gene_symbols.csv"

        self.af = ArtiFusion(input_table, working_dir, reference_genome, bed_file, gene_symbols_file)


    def test_get_longest_transcript(self):
        self.assertEqual(self.af.get_longest_transcript("LARGE1"), "ENST00000476319")

    def test_rev_comp(self):
        self.assertEqual(self.af.rev_comp("ACGT"), "ACGT")

    def test_determine_bp_1(self):
        self.assertEqual(self.af.determine_bp_1([0, 10, 20], [1, 2, 3], "1", 0.5, 6, "+"), ('TAGCT', 5, '1:1:+', [('1', 20, 23), ('1', 10, 12)]))
        self.assertEqual(self.af.determine_bp_1([0, 10, 20], [1, 2, 3], "1", 0.5, 6, "-"), ('TAT', 3, '1:20:-', [('1', 0, 1), ('1', 10, 12)]))

    def test_determine_bp_2(self):
        self.assertEqual(self.af.determine_bp_2([50, 60, 70], [1, 2, 3], "2", 6, "+"), (0, ''))
        self.assertEqual(self.af.determine_bp_2([50, 60, 70], [1, 2, 3], "2", 6, "-"), (0, ''))

    def test_swap(self):
        self.assertEqual(self.af.swap([100, 200, 300], [10, 20, 30], "2", "-", 2, "ACGT", 4), "TATCGACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGGCGCTAGCTACATACGACGATCGATCGATCGATCGATGCGCGATCGATCGATNNNNNNNNNNATCGATCGATCTAGCTAGCGACTAGCGACTGACTGATCGATCGGCGCTACGATCGATCGATCGATCTATCGACTAGCTGACTGATCGACGNNNNNNNNNNNNNNNNNNNNCTAGCTAGCTAGCTAGTTGCGACTAGCATGCTAGCATGCTACTAGCTACGATGCATCTATCTGACGATTACTGACATCGANNNNNNNNNNNNNNNNNNNNNNNNNNACGTCGTAGCTAGCTAGCTAGCTAGCTA")
        self.assertEqual(self.af.swap([100, 200, 300], [10, 20, 30], "2", "+", 2, "ACGT", 4), "TATCGACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGGCGCTAGCTACATACGACGATCGATCGATCGATCGATGCGCGATCGATCGATCGATCGATCGATCGATCGATCTAGCTAGCGACTAGCGACTGACTGATCGATCGGCGCTACGATCGATCGATCGATCTATCGACTAGCTGACTGATCGACGGCATCGATCGGATCAGCTAGCTAGCTAGCTAGCTAGTTGCGACTAGCATGCTAGCATGCTACTAGCTACGATGCATCTATCTGACGATTACTGACATCGAACGTNNNNNNNNNNNNNNNNNNNNNNNNNNCGTAGCTAGCTAGCTAGCTAGCTA")

if __name__ == '__main__':
    unittest.main()
