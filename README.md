# 变异核对与单倍型生成（FastAPI）

同一参考、同一样本下的两个测序数据处理入口：

1. **双 VCF 变异核对**：核对两套变异识别程序输出的单样本 VCF，按**归一化后的等位变异**
   列出两侧一致、仅单侧出现与剂量不符的结果，并保留每条结论的 VCF 原始行号，供人工复查。
2. **单倍型 FASTA 生成**：读取分相（phased）单样本 VCF，按 GT 分别编辑两条参考拷贝，
   输出两条可直接比对的单倍型样本序列（带样本名的两份 FASTA 下载）。

## 功能与判定规则

### 共用输入限制

- 单条染色体参考 FASTA，纯 A/C/G/T 序列，不超过 **1 MiB**；
- 单样本 VCF（恰好 10 列、1 个样本），各不超过 **256 KiB**；
- 所有记录做 **REF 校验**（含坐标越界、CHROM 一致性），indel 对照参考做**左对齐归一化**
  （与 `bcftools norm` / `vt normalize` 等价）。

### 双 VCF 核对

- **接受的变异**：A/C/G/T 单碱基替换（SNV）、事件长度不超过 **20 bp** 的纯插入 / 纯缺失；
  GT 必须是**非缺失、未分相的二倍体**（`0/1`、`1/1`，多 ALT 可为 `1/2`）。
  其他形态（MNV、符号化 ALT、`*`、分相 `|`、缺失 `.`、非二倍体等）一律判失败。
- **剂量计算**：多 ALT 记录按 GT 分别计算每个 ALT 的 0/1/2 剂量，剂量为 0 的变异忽略；
  同一文件内同一归一变异的多条记录剂量累加，超过二倍体上限 2 即判定为**互相矛盾**。
- **整次失败**（返回 422，只显示错误、不展示任何局部结论）：
  任一记录 REF 与参考序列不符；同一文件中同一归一变异剂量互相矛盾；任何不满足输入规则的情形。
- **结果分组**（均按归一变异排序，并标注两侧原始行号与原始 `POS:REF>ALT` 写法）：
  两侧一致 / 仅 VCF A 或 B 出现 / 剂量不符。

### 单倍型 FASTA 生成

- **接受的变异**：仅**单 ALT** 记录——A/C/G/T 单碱基替换，或事件长度不超过 **20 bp** 的
  纯插入 / 纯缺失；REF 校验与插缺归一化与核对流程共用同一实现。
- **GT 必须是二倍体分相写法**（`0|1`、`1|0`、`1|1`、`0|0`），拒绝未分相 `/` 写法；
  **杂合调用（`0|1` / `1|0`）必须携带 PS，且全部杂合调用的 PS 相同**；纯合调用不要求 PS。
- **分拷贝编辑**：`1|0` 编辑拷贝 1，`0|1` 编辑拷贝 2，`1|1` 两条都编辑，`0|0` 只校验不编辑；
  所有变异一律**按原参考坐标定位**。
- **跨度交叠检查**：同一拷贝上变异的参考跨度（**含插缺锚定碱基**）两两不得交叠，
  否则整批失败；分属不同拷贝的交叠调用不受影响，仍可正常生成。
- **整批失败**（返回 422，不交付任何部分序列）：REF 不符、多 ALT、未分相 GT、
  杂合 PS 缺失或不一致、同拷贝跨度交叠，以及任何不满足输入规则的情形。
- **交付**：成功后在结果页提供 `{样本名}_hap1.fa` / `{样本名}_hap2.fa` 两份下载
  （链接 30 分钟内有效），并列出每条变异分配到哪条拷贝。

## 目录结构

```
app/            FastAPI 应用与核心逻辑（均只用标准库）
  main.py       路由、上传校验、下载令牌暂存
  bio.py        FASTA/VCF 解析、REF 校验、左对齐归一化、剂量与分组
  haplo.py      分相 GT/PS 解析、跨度交叠检查、单倍型序列构建
  templates/    上传页、核对结果页、单倍型结果页（Jinja2）
tests/          单元 / 端到端测试（test_bio + test_haplo 仅标准库，test_web 需依赖）
samples/        可直接试用的演示数据（demo.fa + 两份核对 VCF + 分相 VCF）
scripts/        演示 VCF 生成器（坐标由参考自动计算）
Dockerfile, docker-compose.yml, requirements.txt
```

## 用 Docker Compose 启动（推荐）

```bash
docker compose up --build
```

浏览器入口：

- 上传页面（核对 + 单倍型生成两个表单）：<http://localhost:8000/>
- 交互式 API 文档（Swagger UI）：<http://localhost:8000/docs>

停止：`docker compose down`。容器以只读根文件系统、`no-new-privileges` 运行，
内置 `/` 健康检查。

## 本地运行（无 Docker 时）

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## 试用

打开 <http://localhost:8000/>：

- **双 VCF 核对**：参考 `samples/demo.fa`，VCF A `samples/caller_a.vcf`，VCF B `samples/caller_b.vcf`。
  演示数据覆盖：SNV 一致、重复区内插入 / 缺失的两种等价写法归一、多 ALT 按 GT 拆分、
  `0/0` 剂量 0 被忽略、剂量不符（1/1 vs 0/1）、仅单侧出现。
- **单倍型 FASTA 生成**：参考 `samples/demo.fa`，VCF `samples/phased_demo.vcf`。
  演示数据覆盖：`1|0` / `0|1` 杂合（同一 PS=100）、重复区插入与缺失的归一化、
  分属不同拷贝的跨度交叠、`1|1` 纯合双侧编辑、`0|0` 仅校验不编辑。

也可以走 JSON 接口：

```bash
# 双 VCF 核对
curl -F 'fasta=@samples/demo.fa' \
     -F 'vcf_a=@samples/caller_a.vcf' \
     -F 'vcf_b=@samples/caller_b.vcf' \
     http://localhost:8000/api/compare

# 单倍型生成（响应内直接携带两条 FASTA 文本）
curl -F 'fasta=@samples/demo.fa' \
     -F 'vcf=@samples/phased_demo.vcf' \
     http://localhost:8000/api/haplotype
```

## 测试

```bash
# 核心逻辑（仅标准库）
python3 -I -m unittest discover -s tests -t . -p 'test_bio.py' -v
python3 -I -m unittest discover -s tests -t . -p 'test_haplo.py' -v

# 全部测试（Web 层需要先安装 requirements.txt）
python3 -I -m unittest discover -s tests -t . -v
```
