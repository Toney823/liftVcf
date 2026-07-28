# swap_ref — VCF 参考基因组互换工具

**一句话**：把 VCF 中任意一个样品变成参考基因组，原始参考基因组变成一列样品。

## 快速开始

```bash
# 1. 看一眼有哪些样品
python3 -c "import gzip; f=gzip.open('toy.vcf.gz','rt');
[print(l.split()[9]) for l in f if l.startswith('#CHROM')]; break"

# 2. 选一个样品，run
python3 swap_ref.py toy.vcf.gz --sample 1_H3

# 3. 验证结果正确
python3 verify_swap.py toy.vcf.gz --swap --sample 1_H3
```

## 三种模式

| 模式 | 命令 | 行为 |
|------|------|------|
| **strict**（默认·推荐） | `--sample 1_H3` | 仅纯合非 REF 位点交换，杂合位点不动 |
| **iupac** | `--sample 1_H3 --iupac` | 杂合位点 REF 写成 IUPAC 编码（A/C→M, G/T→K...） |
| **force** | `--sample 1_H3 --force` | 杂合位点按测序深度选主等位基因做 REF |

## 输出

- 文件名自动生成：`toy_ref-1_H3.vcf.gz`（`<输入>_ref-<样品名>.vcf.gz`）
- 比原 VCF 多一列 `ORIGINAL_REF`（原始参考基因组的基因型）
- Header 中记录了用哪个样品做的参考：`##swapRef_Sample=<ID=1_H3,...>`

## 验证原理

```
原 VCF 的 IBS 距离矩阵  ==  交换后 VCF 的 IBS 距离矩阵
         ↓                          ↓
      聚类树                      聚类树
         ↓                          ↓
    拓扑结构        ===        拓扑结构
```

等位基因重编码是双射（一一映射），所以样品间的等位基因共享关系不变 → 距离矩阵不变 → 树拓扑不变。

## 全部参数

```
python3 swap_ref.py input.vcf.gz \
    --sample SAMPLE        # 必填：作为新参考的样品名
    -o output.vcf.gz       # 输出文件（默认自动命名）
    --iupac                # 杂合位点 IUPAC 编码
    --force                # 杂合位点按 AD 选主等位基因
    --new-sample-name NAME # 原始参考列名（默认 ORIGINAL_REF）
    --keep-pl              # 保留 PL 字段
    -q                     # 安静模式
```
