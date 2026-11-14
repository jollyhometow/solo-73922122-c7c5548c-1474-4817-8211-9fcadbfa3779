"""分相变异 -> 两条单倍型样本序列 的核心逻辑（仅依赖标准库）。

处理流程：
1. 解析单样本 VCF，仅接受：
   - 单 ALT 的 A/C/G/T 单碱基替换，或事件长度 <= 20 bp 的纯插入 / 纯缺失；
   - 二倍体分相 GT（形如 0|1、1|0、1|1、0|0），拒绝未分相 '/' 写法；
   - 杂合调用（0|1 / 1|0）必须携带 PS，且全部杂合调用的 PS 相同。
2. REF 校验与插缺左对齐归一化沿用 bio 模块的同一实现；
   所有变异一律按原参考坐标定位。
3. 按 GT 把变异分配到两条参考拷贝：1|0 -> 拷贝 1，0|1 -> 拷贝 2，1|1 -> 两条。
4. 同一拷贝上变异的参考跨度（含插缺锚定碱基）两两不得交叠，否则整批失败；
   分属不同拷贝的交叠不受影响。
5. 任何一步失败都抛 CheckError，不交付任何部分序列。
"""

import re
from dataclasses import dataclass

from .bio import (
    CheckError,
    Reference,
    Variant,
    check_alt_bases,
    classify,
    decode_text,
    normalize,
    parse_vcf_header,
    split_format_sample,
    validate_locus,
)

# 输出 FASTA 每行碱基数
FASTA_LINE_WIDTH = 60

# 下载文件名 / FASTA 头中允许出现的字符
_STEM_UNSAFE = re.compile(r"[^A-Za-z0-9._-]+")


@dataclass(frozen=True)
class PhasedCall:
    """一条待应用的分相变异调用（已归一化）。"""

    variant: Variant
    on_copy1: bool
    on_copy2: bool
    line_no: int
    raw: str  # 原始写法 "POS:REF>ALT"


# --------------------------------------------------------------------------- #
# GT / PS
# --------------------------------------------------------------------------- #
def parse_phased_gt(gt: str, line_no: int, n_alt: int) -> tuple[int, int]:
    """解析二倍体分相 GT（形如 0|1），返回两个等位的 ALT 索引（0 表示 REF）。"""
    if not gt:
        raise CheckError(f"第 {line_no} 行：GT 为空。")
    if "/" in gt:
        raise CheckError(
            f"第 {line_no} 行：单倍型生成仅接受分相 GT（'|'），发现未分相标记 '/'：{gt}"
        )

    parts = gt.split("|")
    if len(parts) != 2:
        raise CheckError(f"第 {line_no} 行：仅接受二倍体分相 GT（形如 0|1），收到 {gt!r}。")

    alleles: list[int] = []
    for part in parts:
        if part in ("", "."):
            raise CheckError(f"第 {line_no} 行：GT 含缺失等位，本工具不接受缺失 GT：{gt}")
        try:
            idx = int(part)
        except ValueError:
            raise CheckError(f"第 {line_no} 行：GT 等位 {part!r} 不是整数索引。")
        if idx < 0:
            raise CheckError(f"第 {line_no} 行：GT 等位索引不能为负：{gt}")
        if idx > n_alt:
            raise CheckError(
                f"第 {line_no} 行：GT 等位索引 {idx} 超出 ALT 数量（{n_alt} 个 ALT）。"
            )
        alleles.append(idx)
    return alleles[0], alleles[1]


# --------------------------------------------------------------------------- #
# VCF
# --------------------------------------------------------------------------- #
def _parse_phased_record(
    fields: list[str], line_no: int, refobj: Reference
) -> tuple[PhasedCall | None, str | None]:
    """解析一条数据行，返回 (PhasedCall 或 None, 杂合调用的 PS 或 None)。

    0|0 记录完成全部结构校验后不产生任何拷贝编辑（返回的调用为 None）。
    """
    pos, ref = validate_locus(fields, line_no, refobj)

    alt_s = fields[4]
    if "," in alt_s:
        raise CheckError(
            f"第 {line_no} 行：ALT={alt_s!r} 含多个等位，单倍型生成仅接受单 ALT 记录。"
        )
    check_alt_bases(alt_s, line_no)

    fmt_keys, sample_values = split_format_sample(fields, line_no)
    gt = sample_values[0] if sample_values else ""
    a, b = parse_phased_gt(gt, line_no, 1)

    ps = None
    if "PS" in fmt_keys:
        ps_idx = fmt_keys.index("PS")
        if ps_idx < len(sample_values) and sample_values[ps_idx] not in ("", "."):
            ps = sample_values[ps_idx]

    if a != b:
        # 杂合调用：必须有 PS，由调用方汇总校验全部杂合 PS 一致
        if ps is None:
            raise CheckError(
                f"第 {line_no} 行：杂合分相调用 {gt} 缺少 PS（相位集）标记，"
                "无法确定与其他杂合调用的连锁关系。"
            )
    if a == 0 and b == 0:
        return None, None  # 双拷贝均为 REF：只校验，不编辑

    npos, nref, nalt = normalize(pos, ref, alt_s, refobj.sequence)
    kind = classify(nref, nalt)  # 同时承担形态与长度校验
    variant = Variant(npos, nref, nalt, kind)
    call = PhasedCall(
        variant=variant,
        on_copy1=(a == 1),
        on_copy2=(b == 1),
        line_no=line_no,
        raw=f"{pos}:{ref}>{alt_s}",
    )
    return call, ps if a != b else None


def parse_phased_vcf(
    data: bytes, filename: str, refobj: Reference
) -> tuple[list[PhasedCall], str, str | None]:
    """解析单样本分相 VCF，返回 (调用列表, 样本名, 杂合调用的共同 PS)。

    任何硬性错误（含杂合 PS 缺失或不一致）都抛出 CheckError，整批失败。
    """
    text = decode_text(data, filename)
    lines = text.splitlines()
    header_idx, sample_name = parse_vcf_header(lines, filename)

    calls: list[PhasedCall] = []
    het_ps_values: set[str] = set()
    for offset, raw in enumerate(lines[header_idx + 1:], start=header_idx + 2):
        line = raw.rstrip("\r")
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        call, het_ps = _parse_phased_record(fields, offset, refobj)
        if het_ps is not None:
            het_ps_values.add(het_ps)
        if call is not None:
            calls.append(call)

    if len(het_ps_values) > 1:
        raise CheckError(
            f"{filename}：杂合调用出现多个不同的 PS 相位集"
            f"（{', '.join(sorted(het_ps_values))}），无法确定跨相位集的连锁关系，整批失败。"
        )
    ps = next(iter(het_ps_values), None)
    return calls, sample_name, ps


# --------------------------------------------------------------------------- #
# 单倍型构建
# --------------------------------------------------------------------------- #
def _span(variant: Variant) -> tuple[int, int]:
    """变异在参考上的跨度（1-based 闭区间；插缺含锚定碱基）。"""
    return variant.pos, variant.pos + len(variant.ref) - 1


def _check_copy_overlaps(ordered: list[PhasedCall], copy_idx: int) -> None:
    """同一拷贝上参考跨度两两不得交叠（含插缺锚点），否则整批失败。"""
    prev: PhasedCall | None = None
    for call in ordered:
        if prev is not None:
            _, prev_end = _span(prev.variant)
            cur_start, _ = _span(call.variant)
            if cur_start <= prev_end:
                raise CheckError(
                    f"拷贝 {copy_idx} 上变异 {prev.variant.label}（第 {prev.line_no} 行）"
                    f"与 {call.variant.label}（第 {call.line_no} 行）的参考跨度交叠"
                    "（插缺跨度含锚定碱基），无法确定编辑次序，整批失败。"
                )
        prev = call


def _apply_calls(sequence: str, ordered: list[PhasedCall]) -> str:
    """按原参考坐标把同一拷贝上的变异依次应用到参考序列。"""
    parts: list[str] = []
    cursor = 0  # 0-based，已复制到参考的右端（不含）
    for call in ordered:
        v = call.variant
        start = v.pos - 1
        parts.append(sequence[cursor:start])
        parts.append(v.alt)
        cursor = start + len(v.ref)
    parts.append(sequence[cursor:])
    return "".join(parts)


def build_haplotypes(
    refobj: Reference, calls: list[PhasedCall]
) -> tuple[str, str]:
    """按分相调用分别编辑两条参考拷贝，返回 (拷贝 1 序列, 拷贝 2 序列)。"""
    per_copy: tuple[list[PhasedCall], list[PhasedCall]] = ([], [])
    for call in calls:
        if call.on_copy1:
            per_copy[0].append(call)
        if call.on_copy2:
            per_copy[1].append(call)

    sequences: list[str] = []
    for copy_idx, copy_calls in enumerate(per_copy, start=1):
        ordered = sorted(
            copy_calls, key=lambda c: (c.variant.pos, c.variant.ref, c.variant.alt)
        )
        _check_copy_overlaps(ordered, copy_idx)
        sequences.append(_apply_calls(refobj.sequence, ordered))
    return sequences[0], sequences[1]


# --------------------------------------------------------------------------- #
# 输出
# --------------------------------------------------------------------------- #
def safe_file_stem(name: str) -> str:
    """把样本名整理成安全的文件名片段（仅保留字母数字与 . _ -）。"""
    stem = _STEM_UNSAFE.sub("_", name).strip("._")
    return stem or "sample"


def format_fasta(header: str, sequence: str) -> str:
    """生成换行宽度固定的 FASTA 文本。"""
    lines = [f">{header}"]
    lines.extend(
        sequence[i:i + FASTA_LINE_WIDTH]
        for i in range(0, len(sequence), FASTA_LINE_WIDTH)
    )
    return "\n".join(lines) + "\n"
