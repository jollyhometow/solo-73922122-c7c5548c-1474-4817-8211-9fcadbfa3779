"""单倍型生成核心逻辑测试，仅依赖标准库：python -I -m unittest -v"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.bio import CheckError, parse_fasta  # noqa: E402
from app.haplo import (  # noqa: E402
    build_haplotypes,
    format_fasta,
    parse_phased_vcf,
    safe_file_stem,
)

# 1-based: 1A 2C 3G 4T | 5A 6A 7A 8A | 9C 10G 11T 12T
SEQ = "ACGTAAAACGTT"
FASTA = f">chrT\n{SEQ}\n".encode()


def vcf(records, sample="S1"):
    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample,
    ]
    lines.extend(records)
    return ("\n".join(lines) + "\n").encode()


def rec(pos, ref, alt, gt, ps=None, chrom="chrT", fmt=None):
    if ps is not None:
        return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT:PS\t{gt}:{ps}"
    if fmt is not None:
        return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\t{fmt}"
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT\t{gt}"


def parse(data):
    refobj = parse_fasta(FASTA)
    return parse_phased_vcf(data, "test.vcf", refobj)


def build(records):
    calls, _sample, _ps = parse(vcf(records))
    return build_haplotypes(parse_fasta(FASTA), calls)


class PhasedGtTests(unittest.TestCase):
    def test_unphased_gt_rejected(self):
        for gt in ("0/1", "1/1", "0/0"):
            with self.subTest(gt=gt):
                with self.assertRaises(CheckError):
                    parse(vcf([rec(2, "C", "G", gt)]))

    def test_phased_forms_accepted(self):
        for gt in ("0|1", "1|0", "1|1", "0|0"):
            with self.subTest(gt=gt):
                calls, _, _ = parse(vcf([rec(2, "C", "G", gt, ps=7 if "0|1" == gt or "1|0" == gt else None)]))
                if gt == "0|0":
                    self.assertEqual(calls, [])
                else:
                    self.assertEqual(len(calls), 1)

    def test_missing_allele_rejected(self):
        for gt in (".|1", "1|.", ".|."):
            with self.subTest(gt=gt):
                with self.assertRaises(CheckError):
                    parse(vcf([rec(2, "C", "G", gt, ps=1)]))

    def test_wrong_ploidy_rejected(self):
        for gt in ("1", "0|1|1", "0|1|2"):
            with self.subTest(gt=gt):
                with self.assertRaises(CheckError):
                    parse(vcf([rec(2, "C", "G", gt, ps=1)]))

    def test_alt_index_out_of_range_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(2, "C", "G", "0|2", ps=1)]))

    def test_negative_allele_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(2, "C", "G", "-1|0", ps=1)]))

    def test_empty_gt_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf(["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t"]))


class PhaseSetTests(unittest.TestCase):
    def test_het_without_ps_rejected(self):
        with self.assertRaises(CheckError) as ctx:
            parse(vcf([rec(2, "C", "G", "0|1")]))
        self.assertIn("PS", str(ctx.exception))

    def test_het_with_dot_ps_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(2, "C", "G", "0|1", ps=".")]))

    def test_het_with_different_ps_rejected(self):
        data = vcf([
            rec(2, "C", "G", "0|1", ps=100),
            rec(3, "G", "T", "1|0", ps=200),
        ])
        with self.assertRaises(CheckError) as ctx:
            parse(data)
        self.assertIn("PS", str(ctx.exception))

    def test_het_with_same_ps_accepted(self):
        data = vcf([
            rec(2, "C", "G", "0|1", ps=100),
            rec(3, "G", "T", "1|0", ps=100),
        ])
        calls, _, ps = parse(data)
        self.assertEqual(len(calls), 2)
        self.assertEqual(ps, "100")

    def test_hom_only_needs_no_ps(self):
        calls, _, ps = parse(vcf([rec(2, "C", "G", "1|1")]))
        self.assertEqual(len(calls), 1)
        self.assertIsNone(ps)

    def test_hom_ps_ignored_even_if_different(self):
        data = vcf([
            rec(2, "C", "G", "0|1", ps=100),
            rec(3, "G", "T", "1|1", ps=999),
        ])
        _, _, ps = parse(data)
        self.assertEqual(ps, "100")

    def test_ps_not_first_format_key(self):
        calls, _, ps = parse(vcf([
            "chrT\t2\t.\tC\tG\t.\t.\t.\tGT:DP:PS\t0|1:30:100",
        ]))
        self.assertEqual(len(calls), 1)
        self.assertEqual(ps, "100")


class RecordRuleTests(unittest.TestCase):
    def test_multi_alt_rejected(self):
        with self.assertRaises(CheckError) as ctx:
            parse(vcf([rec(2, "C", "G,T", "0|1", ps=1)]))
        self.assertIn("单 ALT", str(ctx.exception))

    def test_mnv_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(1, "AC", "GT", "0|1", ps=1)]))

    def test_insertion_too_long_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(1, "A", "A" + "C" * 21, "0|1", ps=1)]))

    def test_symbolic_alt_rejected(self):
        for alt in (".", "<DEL>", "G*"):
            with self.subTest(alt=alt):
                with self.assertRaises(CheckError):
                    parse(vcf([rec(2, "C", alt, "0|1", ps=1)]))

    def test_ref_mismatch_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(2, "G", "T", "0|1", ps=1)]))

    def test_ref_mismatch_rejected_even_for_ref_hom(self):
        # 0|0 不产生编辑，但 REF 校验仍然生效
        with self.assertRaises(CheckError):
            parse(vcf([rec(2, "G", "T", "0|0")]))

    def test_ref_hom_skips_variant_level_checks(self):
        # 与核对流程一致：0|0 记录的 ALT 不做形态/长度校验，也不产生调用
        calls, _, _ = parse(vcf([rec(1, "AC", "GT", "0|0")]))
        self.assertEqual(calls, [])

    def test_chrom_mismatch_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(2, "C", "G", "0|1", ps=1, chrom="chrX")]))

    def test_pos_out_of_range_rejected(self):
        with self.assertRaises(CheckError):
            parse(vcf([rec(99, "C", "G", "0|1", ps=1)]))

    def test_multi_sample_header_rejected(self):
        data = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
            "chrT\t2\t.\tC\tG\t.\t.\t.\tGT:PS\t0|1:1\t0|1:1\n"
        ).encode()
        with self.assertRaises(CheckError):
            parse(data)

    def test_sample_name_returned(self):
        _, sample, _ = parse(vcf([rec(2, "C", "G", "1|1")], sample="HG002"))
        self.assertEqual(sample, "HG002")


class BuildHaplotypeTests(unittest.TestCase):
    def test_snv_strand_assignment(self):
        hap1, hap2 = build([
            rec(2, "C", "G", "1|0", ps=1),
            rec(3, "G", "T", "0|1", ps=1),
        ])
        self.assertEqual(hap1, "AGGTAAAACGTT")
        self.assertEqual(hap2, "ACTTAAAACGTT")

    def test_homozygous_applies_to_both(self):
        hap1, hap2 = build([rec(2, "C", "G", "1|1")])
        self.assertEqual(hap1, "AGGTAAAACGTT")
        self.assertEqual(hap2, "AGGTAAAACGTT")

    def test_ref_hom_no_edit(self):
        hap1, hap2 = build([rec(2, "C", "G", "0|0")])
        self.assertEqual(hap1, SEQ)
        self.assertEqual(hap2, SEQ)

    def test_insertion_applied(self):
        hap1, hap2 = build([rec(4, "T", "TA", "1|0", ps=1)])
        self.assertEqual(hap1, "ACGTAAAAACGTT")
        self.assertEqual(hap2, SEQ)

    def test_insertion_normalized_from_alternative_writing(self):
        # 重复区里的另一种写法（pos5 A>AA）归一后与 pos4 T>TA 相同
        hap1, _ = build([rec(5, "A", "AA", "1|0", ps=1)])
        self.assertEqual(hap1, "ACGTAAAAACGTT")

    def test_deletion_applied(self):
        hap1, hap2 = build([rec(4, "TA", "T", "0|1", ps=1)])
        self.assertEqual(hap1, SEQ)
        self.assertEqual(hap2, "ACGTAAACGTT")

    def test_coordinates_stay_on_original_reference(self):
        # 同一拷贝上 4 号位插入 + 9 号位 SNV：SNV 仍按原参考坐标 9 定位
        hap1, hap2 = build([
            rec(4, "T", "TA", "1|0", ps=1),
            rec(9, "C", "G", "1|0", ps=1),
        ])
        self.assertEqual(hap1, "ACGTAAAAAGGTT")
        self.assertEqual(hap2, SEQ)

    def test_same_copy_snv_inside_deletion_span_fails(self):
        with self.assertRaises(CheckError) as ctx:
            build([
                rec(2, "C", "G", "1|0", ps=1),
                rec(2, "CG", "C", "1|0", ps=1),
            ])
        self.assertIn("交叠", str(ctx.exception))

    def test_same_copy_insertion_anchor_overlaps_snv_fails(self):
        # 插入跨度含锚定碱基 pos4，与同位置 SNV 交叠
        with self.assertRaises(CheckError):
            build([
                rec(4, "T", "TA", "1|0", ps=1),
                rec(4, "T", "C", "1|0", ps=1),
            ])

    def test_same_copy_insertion_then_next_base_snv_ok(self):
        hap1, _ = build([
            rec(4, "T", "TA", "1|0", ps=1),
            rec(5, "A", "T", "1|0", ps=1),
        ])
        self.assertEqual(hap1, "ACGTATAAACGTT")

    def test_same_copy_deletion_then_flanking_snv_ok(self):
        hap1, _ = build([
            rec(4, "TA", "T", "1|0", ps=1),   # 跨度 [4,5]
            rec(6, "A", "T", "1|0", ps=1),    # 跨度 [6,6]
        ])
        self.assertEqual(hap1, "ACGTTAACGTT")

    def test_same_copy_deletion_span_covers_snv_fails(self):
        with self.assertRaises(CheckError):
            build([
                rec(4, "TA", "T", "1|0", ps=1),   # 跨度 [4,5]
                rec(5, "A", "T", "1|0", ps=1),    # 跨度 [5,5]
            ])

    def test_two_insertions_same_anchor_same_copy_fail(self):
        with self.assertRaises(CheckError):
            build([
                rec(4, "T", "TA", "1|0", ps=1),
                rec(4, "T", "TC", "1|0", ps=1),
            ])

    def test_same_variant_two_writings_same_copy_fail(self):
        # 同一插入的两种写法落在同一拷贝：归一后跨度相同，视为冲突
        with self.assertRaises(CheckError):
            build([
                rec(4, "T", "TA", "1|0", ps=1),
                rec(5, "A", "AA", "1|0", ps=1),
            ])

    def test_same_variant_on_different_copies_ok(self):
        hap1, hap2 = build([
            rec(4, "T", "TA", "1|0", ps=1),
            rec(5, "A", "AA", "0|1", ps=1),
        ])
        self.assertEqual(hap1, "ACGTAAAAACGTT")
        self.assertEqual(hap2, "ACGTAAAAACGTT")

    def test_overlapping_calls_on_different_copies_ok(self):
        # 拷贝 1 的缺失跨度 [4,5] 与拷贝 2 的 SNV [5,5] 交叠：允许生成
        hap1, hap2 = build([
            rec(4, "TA", "T", "1|0", ps=1),
            rec(5, "A", "C", "0|1", ps=1),
        ])
        self.assertEqual(hap1, "ACGTAAACGTT")
        self.assertEqual(hap2, "ACGTCAAACGTT")

    def test_hom_plus_het_same_variant_conflicts(self):
        # 1|1 与 0|1 使拷贝 2 收到同一插入两次
        with self.assertRaises(CheckError):
            build([
                rec(4, "T", "TA", "1|1"),
                rec(5, "A", "AA", "0|1", ps=1),
            ])

    def test_empty_vcf_gives_two_reference_copies(self):
        hap1, hap2 = build([])
        self.assertEqual(hap1, SEQ)
        self.assertEqual(hap2, SEQ)


class OutputTests(unittest.TestCase):
    def test_format_fasta_wraps_at_60(self):
        seq = "A" * 130
        text = format_fasta("x_hap1", seq)
        lines = text.splitlines()
        self.assertEqual(lines[0], ">x_hap1")
        self.assertEqual([len(ln) for ln in lines[1:]], [60, 60, 10])
        self.assertTrue(text.endswith("\n"))

    def test_safe_file_stem(self):
        self.assertEqual(safe_file_stem("HG002"), "HG002")
        self.assertEqual(safe_file_stem("HG 002/β"), "HG_002")
        self.assertEqual(safe_file_stem(".."), "sample")
        self.assertEqual(safe_file_stem(""), "sample")


if __name__ == "__main__":
    unittest.main()
