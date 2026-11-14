"""FastAPI 入口：

1. 双 VCF 核对：上传参考 FASTA + 两份单样本 VCF，展示归一化核对结果；
2. 单倍型生成：上传参考 FASTA + 一份分相单样本 VCF，生成两条单倍型 FASTA 下载。
"""

import secrets
import threading
import time
from pathlib import Path

from fastapi import FastAPI, File, Request, UploadFile
from fastapi.responses import HTMLResponse, JSONResponse, Response
from fastapi.templating import Jinja2Templates

from .bio import (
    MAX_FASTA_BYTES,
    MAX_VCF_BYTES,
    CheckError,
    compare,
    parse_fasta,
    parse_vcf,
)
from .haplo import (
    build_haplotypes,
    format_fasta,
    parse_phased_vcf,
    safe_file_stem,
)

app = FastAPI(
    title="变异核对与单倍型生成",
    description=(
        "上传单条染色体参考 FASTA 与单样本 VCF："
        "核对两份 VCF 的归一化差异，或由分相变异生成两条单倍型 FASTA。"
    ),
)

templates = Jinja2Templates(directory=str(Path(__file__).parent / "templates"))

# 单倍型下载链接的有效期（秒）与最大在存份数
DOWNLOAD_TTL_SECONDS = 30 * 60
DOWNLOAD_MAX_ENTRIES = 64


class _DownloadStore:
    """内存下载暂存：生成成功的两条单倍型 FASTA 以一次性令牌对外提供下载。

    只在生成全部成功后才写入；失败时不产生任何可下载的部分结果。
    """

    def __init__(self, ttl_seconds: int = DOWNLOAD_TTL_SECONDS, max_entries: int = DOWNLOAD_MAX_ENTRIES):
        self._ttl = ttl_seconds
        self._max = max_entries
        self._items: dict[str, tuple[float, dict]] = {}
        self._lock = threading.Lock()

    def put(self, payload: dict) -> str:
        with self._lock:
            now = time.time()
            for token in [t for t, (ts, _) in self._items.items() if now - ts > self._ttl]:
                del self._items[token]
            while len(self._items) >= self._max:
                oldest = min(self._items, key=lambda t: self._items[t][0])
                del self._items[oldest]
            token = secrets.token_urlsafe(16)
            self._items[token] = (now, payload)
            return token

    def get(self, token: str) -> dict | None:
        with self._lock:
            item = self._items.get(token)
            if item is None:
                return None
            ts, payload = item
            if time.time() - ts > self._ttl:
                del self._items[token]
                return None
            return payload


_downloads = _DownloadStore()


def _read_limited(file: UploadFile, limit: int, label: str) -> bytes:
    data = file.file.read(limit + 1)
    if len(data) > limit:
        raise CheckError(
            f"{label} 超过大小上限 {limit} 字节"
            f"（FASTA 上限 {MAX_FASTA_BYTES // 1024} KiB，VCF 上限 {MAX_VCF_BYTES // 1024} KiB）。"
        )
    if not data:
        raise CheckError(f"{label} 为空文件。")
    return data


def _index_context(**extra):
    ctx = {
        "fasta_limit_kib": MAX_FASTA_BYTES // 1024,
        "vcf_limit_kib": MAX_VCF_BYTES // 1024,
    }
    ctx.update(extra)
    return ctx


# --------------------------------------------------------------------------- #
# 双 VCF 核对
# --------------------------------------------------------------------------- #
def _run_compare(fasta_bytes: bytes, vcf_a_bytes: bytes, vcf_b_bytes: bytes):
    refobj = parse_fasta(fasta_bytes)
    calls_a, sample_a = parse_vcf(vcf_a_bytes, "VCF A", refobj)
    calls_b, sample_b = parse_vcf(vcf_b_bytes, "VCF B", refobj)
    groups = compare(refobj, calls_a, calls_b)
    return refobj, sample_a, sample_b, groups


@app.get("/", response_class=HTMLResponse)
def index(request: Request):
    return templates.TemplateResponse(request, "index.html", _index_context())


@app.post("/compare", response_class=HTMLResponse)
def compare_page(
    request: Request,
    fasta: UploadFile = File(...),
    vcf_a: UploadFile = File(...),
    vcf_b: UploadFile = File(...),
):
    try:
        fasta_bytes = _read_limited(fasta, MAX_FASTA_BYTES, "参考 FASTA")
        vcf_a_bytes = _read_limited(vcf_a, MAX_VCF_BYTES, "VCF A")
        vcf_b_bytes = _read_limited(vcf_b, MAX_VCF_BYTES, "VCF B")
        refobj, sample_a, sample_b, groups = _run_compare(
            fasta_bytes, vcf_a_bytes, vcf_b_bytes
        )
    except CheckError as exc:
        # 整次核对失败：只展示错误，不展示任何局部结论
        return templates.TemplateResponse(
            request,
            "index.html",
            _index_context(error=str(exc), error_title="核对失败"),
            status_code=422,
        )

    return templates.TemplateResponse(
        request,
        "result.html",
        {
            "ref_name": refobj.name,
            "ref_length": len(refobj.sequence),
            "sample_a": sample_a,
            "sample_b": sample_b,
            "groups": groups,
        },
    )


@app.post("/api/compare")
def compare_api(
    fasta: UploadFile = File(...),
    vcf_a: UploadFile = File(...),
    vcf_b: UploadFile = File(...),
):
    """JSON 接口，便于脚本化调用；失败返回 422 与错误信息。"""
    try:
        fasta_bytes = _read_limited(fasta, MAX_FASTA_BYTES, "参考 FASTA")
        vcf_a_bytes = _read_limited(vcf_a, MAX_VCF_BYTES, "VCF A")
        vcf_b_bytes = _read_limited(vcf_b, MAX_VCF_BYTES, "VCF B")
        refobj, sample_a, sample_b, groups = _run_compare(
            fasta_bytes, vcf_a_bytes, vcf_b_bytes
        )
    except CheckError as exc:
        return JSONResponse(status_code=422, content={"ok": False, "error": str(exc)})

    def dump_entry(entry):
        def dump_call(call):
            if call is None:
                return None
            return {
                "dosage": call.dosage,
                "line_numbers": call.line_numbers,
                "raw": call.raw_descriptions,
            }

        v = entry["variant"]
        return {
            "variant": {
                "pos": v.pos,
                "ref": v.ref,
                "alt": v.alt,
                "kind": v.kind,
                "label": v.label,
            },
            "left": dump_call(entry["left"]),
            "right": dump_call(entry["right"]),
        }

    return {
        "ok": True,
        "reference": {"name": refobj.name, "length": len(refobj.sequence)},
        "samples": {"left": sample_a, "right": sample_b},
        "result": {key: [dump_entry(e) for e in rows] for key, rows in groups.items()},
    }


# --------------------------------------------------------------------------- #
# 单倍型 FASTA 生成
# --------------------------------------------------------------------------- #
def _run_haplotype(fasta_bytes: bytes, vcf_bytes: bytes):
    """解析 -> 构建 -> 格式化全部完成后才返回；任何一步失败都抛 CheckError。"""
    refobj = parse_fasta(fasta_bytes)
    calls, sample_name, ps = parse_phased_vcf(vcf_bytes, "VCF", refobj)
    hap1, hap2 = build_haplotypes(refobj, calls)
    stem = safe_file_stem(sample_name)
    fasta1 = format_fasta(f"{stem}_hap1 reference={refobj.name}", hap1)
    fasta2 = format_fasta(f"{stem}_hap2 reference={refobj.name}", hap2)
    return refobj, calls, sample_name, stem, ps, hap1, hap2, fasta1, fasta2


@app.post("/haplotype", response_class=HTMLResponse)
def haplotype_page(
    request: Request,
    fasta: UploadFile = File(...),
    vcf: UploadFile = File(...),
):
    try:
        fasta_bytes = _read_limited(fasta, MAX_FASTA_BYTES, "参考 FASTA")
        vcf_bytes = _read_limited(vcf, MAX_VCF_BYTES, "VCF")
        refobj, calls, sample_name, stem, ps, hap1, hap2, fasta1, fasta2 = _run_haplotype(
            fasta_bytes, vcf_bytes
        )
    except CheckError as exc:
        # 整批失败：只展示错误，不交付任何部分序列
        return templates.TemplateResponse(
            request,
            "index.html",
            _index_context(error=str(exc), error_title="单倍型生成失败"),
            status_code=422,
        )

    token = _downloads.put({"stem": stem, "fasta1": fasta1, "fasta2": fasta2})
    return templates.TemplateResponse(
        request,
        "haplotype_result.html",
        {
            "ref_name": refobj.name,
            "ref_length": len(refobj.sequence),
            "sample": sample_name,
            "ps": ps,
            "calls": calls,
            "n_copy1": sum(1 for c in calls if c.on_copy1),
            "n_copy2": sum(1 for c in calls if c.on_copy2),
            "hap1_len": len(hap1),
            "hap2_len": len(hap2),
            "token": token,
            "ttl_minutes": DOWNLOAD_TTL_SECONDS // 60,
        },
    )


@app.get("/haplotype/download/{token}/{copy}")
def haplotype_download(token: str, copy: int):
    """按令牌下载一条单倍型 FASTA；令牌过期或不存在返回 404。"""
    item = _downloads.get(token)
    if item is None or copy not in (1, 2):
        return HTMLResponse(
            "<p>下载链接不存在或已过期，请重新生成。</p><p><a href='/'>返回首页</a></p>",
            status_code=404,
        )
    content = item["fasta1"] if copy == 1 else item["fasta2"]
    filename = f"{item['stem']}_hap{copy}.fa"
    return Response(
        content,
        media_type="text/plain; charset=utf-8",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


@app.post("/api/haplotype")
def haplotype_api(
    fasta: UploadFile = File(...),
    vcf: UploadFile = File(...),
):
    """JSON 接口：直接返回两条单倍型 FASTA 文本；失败返回 422 与错误信息。"""
    try:
        fasta_bytes = _read_limited(fasta, MAX_FASTA_BYTES, "参考 FASTA")
        vcf_bytes = _read_limited(vcf, MAX_VCF_BYTES, "VCF")
        refobj, calls, sample_name, stem, ps, hap1, hap2, fasta1, fasta2 = _run_haplotype(
            fasta_bytes, vcf_bytes
        )
    except CheckError as exc:
        return JSONResponse(status_code=422, content={"ok": False, "error": str(exc)})

    return {
        "ok": True,
        "reference": {"name": refobj.name, "length": len(refobj.sequence)},
        "sample": sample_name,
        "ps": ps,
        "variants": [
            {
                "line_no": c.line_no,
                "label": c.variant.label,
                "kind": c.variant.kind,
                "raw": c.raw,
                "on_copy1": c.on_copy1,
                "on_copy2": c.on_copy2,
            }
            for c in calls
        ],
        "haplotypes": {
            "hap1": {"length": len(hap1), "fasta": fasta1},
            "hap2": {"length": len(hap2), "fasta": fasta2},
        },
    }
