#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# Created by Roberto Preste and Domenico Simone
import os
from types import SimpleNamespace
import unittest

from ..mtVariantCaller import (
    parse_mismatches_from_cigar_md,
    allele_strand_counter,
    allele_strand_updater,
    mismatch_detection,
)

r = "A00181:108:HLFMYDSXX:2:2205:3495:11303\t161\tk127_149\t569\t40\t146M1D5M\t=\t564\t-157\tCGTAACGGTTTGCTCCGTCTGACACGGCGGTTCCTTATCGAGTTGGTGTTCCCGGGCATCGTGGGCGCCGGGGGAGTTGTGAATGGCGGTACAATACCCGACGAGGAAAAATACCATGATGTTTACGCGCCATTTCATGTTGATGAGATCG\tFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFF\tAS:i:-23\tXN:i:0\tXM:i:3\tXO:i:1\tXG:i:1\tNM:i:4\tMD:Z:32T113^C1G2T0\tYS:i:-28\tYT:Z:DP"

mismatch_dict = {601 : SimpleNamespace(DP=24, POS=601, REF='T', allele_DP=[24], allele_strand_count=[[24, 0]], alleles=['C'])}
test_sam = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "data",
    "test.sam"
)

class TestSNPcalling(unittest.TestCase):

    def setUp(self) -> None:
        self.r = r
        self.mismatch_dict = mismatch_dict

    def test_parse_mismatches_from_cigar_md_default(self):
        # Given
        expected = ([601], [32], ['T'], ['C'], [37], '+')
        
        # When
        result = parse_mismatches_from_cigar_md(self.r)
        
        # Then
        self.assertEqual(expected, result)

    def test_parse_mismatches_from_cigar_md_tail_zero(self):
        # Given
        expected = ([601, 717, 720], [32, 148, 151], ['T', 'G', 'T'], ['C', 'A', 'G'], [37, 37, 37], '+')
        
        # When
        result = parse_mismatches_from_cigar_md(self.r, tail_mismatch=0)
        
        # Then
        self.assertEqual(expected, result)
    
    def test_parse_mismatches_from_cigar_md_depth_positions(self):
        # the read counts towards the depth where a mismatch would be counted:
        # 146M1D5M at 569, so read offsets 5-145 (the first/last 5 bases and
        # the deletion at 715 excluded), minus the two Q11 bases (',')
        result = parse_mismatches_from_cigar_md(self.r, return_depth_positions=True)
        depth_positions = list(result[-1])
        expected = [p for p in range(574, 715) if p not in (605, 705)]
        self.assertEqual(expected, depth_positions)
        # every counted mismatch lies on a position counted in the depth
        self.assertTrue(set(result[0]) <= set(depth_positions))

    def test_allele_strand_counter(self):
        # Given
        expected = [[1, 0], [0, 1]]
        
        # When
        result = [allele_strand_counter("+"), allele_strand_counter("-")]
        
        # Then
        self.assertEqual(expected, result)

    def test_allele_strand_updater(self):
        # Given
        expected = [25, 5]
        # When
        result = allele_strand_updater([0, 1], [25, 4])
        # Then
        self.assertEqual(expected, result)

    def test_mismatch_detection(self):
        # Given
        expected = self.mismatch_dict
        # When
        result = mismatch_detection(sam=test_sam)
        # Then
        self.assertEqual(expected, result)
        
    # 
    # def test_get_mt_genomes(self):
    #     # Given
    #     expected = ["NC_001323.1"]
    # 
    #     # When
    #     result = get_mt_genomes(self.analysis_tab)
    # 
    #     # Then
    #     self.assertEqual(expected, result)

    # def test_fastqc_outputs_raw(self):
    #     # Given
    #     expected = [el.format("results/fastqc_raw")
    #                 for el in FASTQC_OUT]
    # 
    #     # When
    #     result = fastqc_outputs(datasets_tab=self.datasets_tab,
    #                             analysis_tab=self.analysis_tab,
    #                             out="raw")
    # 
    #     # Then
    #     self.assertEqual(expected, result)
    # 
    # def test_fastqc_outputs_filtered(self):
    #     # Given
    #     expected = [el.format("results/fastqc_filtered")
    #                 for el in FASTQC_OUT_FILT]
    # 
    #     # When
    #     result = fastqc_outputs(datasets_tab=self.datasets_tab,
    #                             analysis_tab=self.analysis_tab,
    #                             out="filtered")
    # 
    #     # Then
    #     self.assertEqual(sorted(expected), sorted(result))
    # 
    # def test_fastqc_outputs_error(self):
    #     # Given/When
    #     with self.assertRaises(ValueError):
    #         fastqc_outputs(datasets_tab=self.datasets_tab,
    #                        analysis_tab=self.analysis_tab,
    #                        out="test")
    # 
    # def test_get_genome_vcf_files(self):
    #     # Given
    #     expected = ["results/vcf/NC_001323.1_GCF_000002315.5.vcf"]
    # 
    #     # When
    #     result = get_genome_vcf_files(df=self.analysis_tab)
    # 
    #     # Then
    #     self.assertEqual(expected, result)
    # 
    # def test_get_bed_files(self):
    #     # Given
    #     expected = [
    #         "results/5517_hypo/5517_hypo_NC_001323.1_GCF_000002315.5.bed",
    #         "results/5517_liver/5517_liver_NC_001323.1_GCF_000002315.5.bed"
    #     ]
    # 
    #     # When
    #     result = get_bed_files(df=self.analysis_tab)
    # 
    #     # Then
    #     self.assertEqual(expected, result)
    # 
    # def test_get_fasta_files(self):
    #     # Given
    #     expected = [
    #         "results/5517_hypo/5517_hypo_NC_001323.1_GCF_000002315.5.fasta",
    #         "results/5517_liver/5517_liver_NC_001323.1_GCF_000002315.5.fasta"
    #     ]
    # 
    #     # When
    #     result = get_fasta_files(df=self.analysis_tab)
    # 
    #     # Then
    #     self.assertEqual(expected, result)
    # 
    # def test_get_genome_files(self):
    #     # Given
    #     expected = ["GCF_000002315.5_GRCg6a_genomic_mt.fna"]
    # 
    #     # When
    #     result = get_genome_files(df=self.reference_tab,
    #                               ref_genome_mt="NC_001323.1",
    #                               field="ref_genome_mt_file")
    # 
    #     # Then
    #     self.assertEqual(expected, result)


def test_get_consensus_single_no_variants():
    # a sample identical to the reference has no variants
    from modules.mtVariantCaller import get_consensus_single
    assert get_consensus_single([]) == []


def test_indel_region_end():
    from modules.mtVariantCaller import indel_region_end
    # CA deleted after the T at 4, in the CACACA repeat at 5-10: reads must
    # reach position 11 (G) to show whether the deletion is present
    seq = "GGGTCACACAGTTT"
    assert indel_region_end(seq, 4, "CA", "del") == 11
    # C inserted after the T at 3, in the poly-C at 4-7: first base after is 8
    assert indel_region_end("AATCCCCGT", 3, "C", "ins") == 8
    # no repeat: the base right after the indel
    assert indel_region_end("AATGCAT", 3, "C", "ins") == 4


def test_read_indels():
    from modules.mtVariantCaller import read_indels
    seq = "NNNNN" + "A" * 10 + "C" * 5 + "G" + "T" * 4
    fields = ["r1", "0", "chrM", "100", "60", "5S10M2D5M1I4M", "*", "0", "0", seq, "I" * len(seq)]
    # deletion after 109 (100-109 aligned), insertion of G after 116
    assert read_indels(fields) == [("del", 109, 2), ("ins", 116, "G")]


def test_parse_indels_no_indel_left():
    # all indel reads already discarded for low flanking quality ('delete'):
    # this used to fail inverting an empty (float) numpy array
    import pandas as pd
    from modules.mtVariantCaller import parse_indels
    df = pd.DataFrame([["Del", "r1", "+", 3106, "range(3107, 3108)", ["delete", "delete"]]],
                      columns=["Type", "readName", "strand", "rleft", "genotype", "qs_flanking"])
    assert parse_indels(df, 25, 5, "qs_flanking").empty


def test_fasta_output_two_deletions_same_position(tmp_path):
    # two deletion alleles at the same position, both above hf_max: this
    # used to fail looking for a variant of another type to drop
    from modules.BEDoutput import fasta_output
    ref = "ACGTACGTAC" * 3
    deletions = [[5, ["TA"], 100, ["T"], [95], [["50;45"]], [30], [0.95], [0.9], [1.0], "del"],
                 [5, ["TAC"], 100, ["T"], [92], [["46;46"]], [30], [0.92], [0.9], [1.0], "del"]]
    snp = [12, "G", 100, ["C"], [99], [["50;49"]], "PASS", [0.99], [0.95], [1.0], "mism"]
    out = tmp_path / "c.fasta"
    fasta_output(vcf_dict={"s": deletions + [snp]}, contigs=[((1, 30), ref)], fasta_out=str(out))
    seq = "".join(out.read_text().splitlines()[1:])
    # the first deletion (of position 6) and the substitution at 12 are applied
    assert seq == ref[:5] + ref[6:11] + "C" + ref[12:]
