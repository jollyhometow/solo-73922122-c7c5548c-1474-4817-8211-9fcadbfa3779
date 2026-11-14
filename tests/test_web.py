"""Web 层端到端测试，需要 fastapi/httpx 依赖：
PYTHONPATH=/tmp/pylibs python -I -m unittest -v tests.test_web
"""

import io
import re
import sys
import unittest
from pathlib import Path

# 本地无 pip 环境时，依赖安装在隔离目录；容器内该目录不存在则跳过
for _cand in (Path("/tmp/pylibs"), Path("/workspace/.pylibs")):
    if _cand.is_dir() and str(_cand) not in sys.path:
        sys.path.insert(0, str(_cand))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from fastapi.testclient import TestClient  # noqa: E402

from app.main import app  # noqa: E402
from app.bio import MAX_FASTA_BYTES, MAX_VCF_BYTES  # noqa: E402

client = TestClient(app)

FASTA = b">chrT\nACGTAAAACGTT\n"
HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"


def vcf(rows):
    return (HEADER + "".join(r + "\n" for r in rows)).encode()


def upload(fasta=FASTA, a=None, b=None):
    return client.post(
        "/compare",
        files={
            "fasta": ("ref.fa", fasta, "text/plain"),
            "vcf_a": ("a.vcf", a if a is not None else b"", "text/plain"),
            "vcf_b": ("b.vcf", b if b is not None else b"", "text/plain"),
        },
    )


class WebTests(unittest.TestCase):
    def test_index_renders(self):
        resp = client.get("/")
        self.assertEqual(resp.status_code, 200)
        self.assertIn("双 VCF 变异核对".encode(), resp.content)

    def test_happy_path_groups_and_line_numbers(self):
        a = vcf([
            "chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1",       # SNV 一致
            "chrT\t4\t.\tT\tTA\t.\t.\t.\tGT\t0/1",      # 插入（写法一）
            "chrT\t10\t.\tG\tT\t.\t.\t.\tGT\t1/1",      # 仅 A 有
        ])
        b = vcf([
            "chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1",
            "chrT\t5\t.\tA\tAA\t.\t.\t.\tGT\t0/1",      # 同一插入的另一种写法
            "chrT\t3\t.\tG\tC\t.\t.\t.\tGT\t0/1",       # 仅 B 有
        ])
        resp = upload(a=a, b=b)
        self.assertEqual(resp.status_code, 200)
        body = resp.text
        self.assertIn("两侧一致", body)
        self.assertIn("4:TinsTA", body)          # 归一变异
        self.assertIn("原始行号：4", body)        # 插入记录的物理行号
        self.assertIn("5:A&gt;AA", body)         # B 侧保留原始写法
        self.assertIn("仅单侧出现", body)
        self.assertIn("剂量不符", body)

    def test_ref_mismatch_shows_failure_page_without_partial_results(self):
        a = vcf(["chrT\t2\t.\tG\tT\t.\t.\t.\tGT\t0/1"])  # pos2 参考是 C
        b = vcf(["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1"])
        resp = upload(a=a, b=b)
        self.assertEqual(resp.status_code, 422)
        self.assertIn("核对失败", resp.text)
        self.assertNotIn("核对概览", resp.text)

    def test_size_limit_enforced(self):
        big = b">chrT\n" + b"A" * (MAX_FASTA_BYTES + 10)
        resp = upload(fasta=big, a=vcf([]), b=vcf([]))
        self.assertEqual(resp.status_code, 422)
        self.assertIn("超过大小上限", resp.text)

        big_vcf = HEADER.encode() + b"#\n" * (MAX_VCF_BYTES // 2)
        resp = upload(a=big_vcf, b=vcf([]))
        self.assertEqual(resp.status_code, 422)

    def test_empty_file_rejected(self):
        resp = upload(a=b"", b=vcf([]))
        self.assertEqual(resp.status_code, 422)
        self.assertIn("为空", resp.text)

    def test_api_json_happy_path(self):
        a = vcf([
            "chrT\t2\t.\tC\tG,T\t.\t.\t.\tGT\t1/2",
            "chrT\t4\t.\tT\tTA\t.\t.\t.\tGT\t1/1",
        ])
        b = vcf([
            "chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1",  # G 剂量不符（1 vs 2）
            "chrT\t4\t.\tT\tTA\t.\t.\t.\tGT\t1/1",
        ])
        resp = client.post(
            "/api/compare",
            files={
                "fasta": ("ref.fa", FASTA, "text/plain"),
                "vcf_a": ("a.vcf", a, "text/plain"),
                "vcf_b": ("b.vcf", b, "text/plain"),
            },
        )
        self.assertEqual(resp.status_code, 200)
        data = resp.json()
        self.assertTrue(data["ok"])
        labels = {e["variant"]["label"] for e in data["result"]["consistent"]}
        self.assertIn("4:TinsTA", labels)
        # 多 ALT 记录 1/2：G、T 各剂量 1；B 侧只有 G 剂量 1 → T 仅 A 有，G 一致
        mismatch = data["result"]["dosage_mismatch"]
        self.assertEqual(mismatch, [])
        left_only = {e["variant"]["label"]: e for e in data["result"]["left_only"]}
        self.assertIn("2:C>T", left_only)
        self.assertEqual(left_only["2:C>T"]["left"]["dosage"], 1)

    def test_api_json_dosage_mismatch(self):
        a = vcf(["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t1/1"])
        b = vcf(["chrT\t2\t.\tC\tG\t.\t.\t.\tGT\t0/1"])
        resp = client.post(
            "/api/compare",
            files={
                "fasta": ("ref.fa", FASTA, "text/plain"),
                "vcf_a": ("a.vcf", a, "text/plain"),
                "vcf_b": ("b.vcf", b, "text/plain"),
            },
        )
        self.assertEqual(resp.status_code, 200)
        mismatch = resp.json()["result"]["dosage_mismatch"]
        self.assertEqual(len(mismatch), 1)
        self.assertEqual(mismatch[0]["left"]["dosage"], 2)
        self.assertEqual(mismatch[0]["right"]["dosage"], 1)

    def test_api_json_failure(self):
        resp = client.post(
            "/api/compare",
            files={
                "fasta": ("ref.fa", FASTA, "text/plain"),
                "vcf_a": ("a.vcf", vcf(["chrT\t2\t.\tG\tT\t.\t.\t.\tGT\t0/1"]), "text/plain"),
                "vcf_b": ("b.vcf", vcf([]), "text/plain"),
            },
        )
        self.assertEqual(resp.status_code, 422)
        self.assertFalse(resp.json()["ok"])


# --------------------------------------------------------------------------- #
# 单倍型 FASTA 生成
# --------------------------------------------------------------------------- #
PHASED_HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDEMO\n"
)


def pvcf(rows):
    return (PHASED_HEADER + "".join(r + "\n" for r in rows)).encode()


def upload_haplo(fasta=FASTA, vcf_bytes=None):
    return client.post(
        "/haplotype",
        files={
            "fasta": ("ref.fa", fasta, "text/plain"),
            "vcf": ("p.vcf", vcf_bytes if vcf_bytes is not None else b"", "text/plain"),
        },
    )


class HaplotypeWebTests(unittest.TestCase):
    HAPPY_ROWS = [
        "chrT\t2\t.\tC\tG\t.\t.\t.\tGT:PS\t1|0:100",   # SNV -> 拷贝 1
        "chrT\t3\t.\tG\tT\t.\t.\t.\tGT:PS\t0|1:100",   # SNV -> 拷贝 2
        "chrT\t9\t.\tC\tG\t.\t.\t.\tGT\t1|1",          # 纯合 -> 两条拷贝
        "chrT\t10\t.\tG\tA\t.\t.\t.\tGT\t0|0",         # 不编辑，仅校验
    ]

    def test_index_offers_both_entries(self):
        resp = client.get("/")
        self.assertEqual(resp.status_code, 200)
        self.assertIn('action="/compare"', resp.text)
        self.assertIn('action="/haplotype"', resp.text)

    def test_happy_path_two_downloads_with_sample_name(self):
        resp = upload_haplo(vcf_bytes=pvcf(self.HAPPY_ROWS))
        self.assertEqual(resp.status_code, 200)
        self.assertIn("DEMO", resp.text)
        self.assertIn("PS：100", resp.text)

        links = re.findall(r'/haplotype/download/([A-Za-z0-9_-]+)/([12])', resp.text)
        self.assertEqual(len(links), 2)
        self.assertEqual(links[0][0], links[1][0])  # 同一令牌
        token = links[0][0]

        r1 = client.get(f"/haplotype/download/{token}/1")
        r2 = client.get(f"/haplotype/download/{token}/2")
        self.assertEqual(r1.status_code, 200)
        self.assertEqual(r2.status_code, 200)
        self.assertIn('filename="DEMO_hap1.fa"', r1.headers["content-disposition"])
        self.assertIn('filename="DEMO_hap2.fa"', r2.headers["content-disposition"])
        # 拷贝 1：pos2 C>G + pos9 C>G；拷贝 2：pos3 G>T + pos9 C>G
        self.assertEqual(r1.text, ">DEMO_hap1 reference=chrT\nAGGTAAAAGGTT\n")
        self.assertEqual(r2.text, ">DEMO_hap2 reference=chrT\nACTTAAAAGGTT\n")

    def test_failure_delivers_no_partial_sequences(self):
        bad = pvcf([
            "chrT\t2\t.\tC\tG\t.\t.\t.\tGT:PS\t1|0:100",
            "chrT\t3\t.\tG\tT\t.\t.\t.\tGT\t0/1",       # 未分相：整批失败
        ])
        resp = upload_haplo(vcf_bytes=bad)
        self.assertEqual(resp.status_code, 422)
        self.assertIn("单倍型生成失败", resp.text)
        self.assertNotIn("/haplotype/download/", resp.text)

    def test_same_copy_overlap_fails(self):
        bad = pvcf([
            "chrT\t4\t.\tT\tTA\t.\t.\t.\tGT:PS\t1|0:7",  # 插入，跨度含锚点 [4,4]
            "chrT\t4\t.\tT\tC\t.\t.\t.\tGT:PS\t1|0:7",   # 同拷贝同位 SNV
        ])
        resp = upload_haplo(vcf_bytes=bad)
        self.assertEqual(resp.status_code, 422)
        self.assertIn("交叠", resp.text)
        self.assertNotIn("/haplotype/download/", resp.text)

    def test_cross_copy_overlap_still_generates(self):
        ok = pvcf([
            "chrT\t4\t.\tTA\tT\t.\t.\t.\tGT:PS\t1|0:7",  # 拷贝 1 缺失，跨度 [4,5]
            "chrT\t5\t.\tA\tC\t.\t.\t.\tGT:PS\t0|1:7",   # 拷贝 2 SNV，跨度 [5,5]
        ])
        resp = upload_haplo(vcf_bytes=ok)
        self.assertEqual(resp.status_code, 200)
        self.assertIn("/haplotype/download/", resp.text)

    def test_download_with_unknown_token_404(self):
        resp = client.get("/haplotype/download/not-a-token/1")
        self.assertEqual(resp.status_code, 404)

    def test_api_haplotype_happy_path(self):
        resp = client.post(
            "/api/haplotype",
            files={
                "fasta": ("ref.fa", FASTA, "text/plain"),
                "vcf": ("p.vcf", pvcf(self.HAPPY_ROWS), "text/plain"),
            },
        )
        self.assertEqual(resp.status_code, 200)
        data = resp.json()
        self.assertTrue(data["ok"])
        self.assertEqual(data["sample"], "DEMO")
        self.assertEqual(data["ps"], "100")
        self.assertEqual(len(data["variants"]), 3)  # 0|0 不产生编辑
        self.assertIn("AGGTAAAAGGTT", data["haplotypes"]["hap1"]["fasta"])
        self.assertIn("ACTTAAAAGGTT", data["haplotypes"]["hap2"]["fasta"])

    def test_api_haplotype_failure(self):
        resp = client.post(
            "/api/haplotype",
            files={
                "fasta": ("ref.fa", FASTA, "text/plain"),
                "vcf": ("p.vcf", pvcf(["chrT\t2\t.\tC\tG,T\t.\t.\t.\tGT:PS\t0|1:1"]), "text/plain"),
            },
        )
        self.assertEqual(resp.status_code, 422)
        self.assertFalse(resp.json()["ok"])


if __name__ == "__main__":
    unittest.main()
