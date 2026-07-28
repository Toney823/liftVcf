# swap_ref — VCF 参考基因组互换工具

**一句话**：把 VCF 中任意一个样品变成参考基因组，原始参考基因组变成一列样品。

> [English](README.md)

## 快速开始

```bash
# 1. 看一眼有哪些样品
python3 -c "import gzip; f=gzip.open('toy.vcf.gz','rt');
[print(l.split()[9]) for l in f if l.startswith('#CHROM')]; break"

# 2. 选一个样品，run（不需要任何额外参数）
python3 swap_ref.py toy.vcf.gz --sample 1_H3

# 3. 验证结果正确
python3 verify_swap.py toy.vcf.gz --swap --sample 1_H3
```

## 它做了什么？

假设你选了样品 `1_H3`，对每个变异位点：

| 1_H3 的基因型 | 含义 | 结果 |
|---|---|---|
| `0/0` | 和参考一样 | 不动 |
| `1/1` | 纯合突变 | **交换**：突变碱基变成新 REF，原 REF 变成新 ALT |
| `0/1` | 杂合（两种碱基都有） | **不动**：样品自己都不确定，没法定义"参考" |
| `./.` | 缺失 | 跳过 |

交换后，VCF 多出一列 `ORIGINAL_REF`，代表原始参考基因组的基因型。

## 选参数？不用

**绝大多数情况不需要任何额外参数**，默认行为就是最正确的：

```bash
python3 swap_ref.py input.vcf.gz --sample 样品名
```

## `--iupac` 和 `--force` 是什么？（可选，非必须）

这两个参数只影响**杂合位点**（样品基因型是 `0/1` 的那些位点）。

默认情况下，杂合位点不动——因为样品自己携带两种碱基，没法说哪个才是"参考"。但如果你非要在这些位点也改 REF，就面临一个问题：**选哪个碱基做新 REF？**

| 参数 | 策略 | 举例 |
|------|------|------|
| （默认，不加） | 杂合位点不动 | REF=A, ALT=G, 样品 0/1 → 保持 REF=A |
| `--iupac` | 用 IUPAC 歧义编码同时表示两个碱基 | 样品 0/1 (A/G) → REF 写成 `R`（=A 或 G） |
| `--force` | 看测序深度，reads 多的那个碱基做 REF | 样品 0/1, AD=3,15 → ALT 深度更高 → REF 变 G |

**什么时候用？**
- `--iupac`：你想让 REF 如实反映"这个样品在此处有多态性"。注意：下游工具（bwa、GATK 等）不一定支持 IUPAC 碱基，慎用。
- `--force`：你不在乎杂合，只想尽可能多地把 REF 改成样品的等位基因。
- 两个都不加：默认，最安全，适用于绝大多数场景。

## 输出

- 文件名自动生成：`toy_ref-1_H3.vcf.gz`（`<输入>_ref-<样品名>.vcf.gz`）
- 比原 VCF 多一列 `ORIGINAL_REF`（原始参考基因组的基因型）
- Header 中记录了用哪个样品做的参考：`##swapRef_Sample=<ID=1_H3,...>`

## 验证原理

等位基因重编码是双射（一一映射），样品间的等位基因共享关系不变，因此：

```
原 VCF 的 IBS 距离矩阵  ==  交换后 VCF 的 IBS 距离矩阵
         ↓                          ↓
      聚类树                      聚类树
         ↓                          ↓
    拓扑结构        ===        拓扑结构
```

距离矩阵不变 → 树拓扑不变。`verify_swap.py` 自动完成以上检查。

## 全部参数

```
python3 swap_ref.py input.vcf.gz \
    --sample SAMPLE        # 必填：作为新参考的样品名
    -o output.vcf.gz       # 输出文件（默认自动命名）
    --iupac                # 可选：杂合位点用 IUPAC 歧义碱基
    --force                # 可选：杂合位点按 AD 选主等位基因
    --new-sample-name NAME # 原始参考列名（默认 ORIGINAL_REF）
    --keep-pl              # 保留 PL 字段
    -q                     # 安静模式
```
