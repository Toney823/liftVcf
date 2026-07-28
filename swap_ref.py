#!/usr/bin/env python3
"""
swap_ref.py — Reorient a multi-sample VCF to a selected sample's genotype
=========================================================================

以指定样品的基因型为参照，重新定向多样品 VCF 的 REF/ALT 等位基因：
  - 该样品纯合非 REF 的位点：该样品的等位基因 → 新 REF，旧 REF → 新 ALT
  - 杂合位点：默认保持不动（strict 模式）
  - 新增 "ORIGINAL_REF" 列：原始参考基因组作为一个额外样品保留
  - 所有样品的等位基因共享关系保持不变（IBS 距离矩阵不变）

用法:
  python swap_ref.py toy.vcf.gz --sample 1_A1 -o swapped.vcf.gz

选项:
  --sample NAME       目标样品名（必填）
  -o, --output FILE   输出文件（默认自动命名）
  --iupac             [EXPERIMENTAL] 杂合位点用 IUPAC 歧义碱基编码 REF
  --force             杂合位点按 AD 深度选主等位基因做新 REF（启发式，不推断相位）
  --new-sample-name   原始参考列的名称（默认 ORIGINAL_REF）
  -q, --quiet         不输出统计信息到 stderr

注意: PL 字段始终丢弃（等位基因重排后需重新计算，建议用 bcftools 补算）

推荐默认用法:
  python swap_ref.py input.vcf.gz --sample SAMPLE
  （默认 strict 模式：仅纯合非 REF 位点交换，最安全）
"""

import sys
import gzip
import argparse

# ============================================================
# IUPAC 碱基编码表（双碱基 → 歧义编码 & 反向查找）
# ============================================================
IUPAC_PAIRS = {
    frozenset('AG'): 'R', frozenset('CT'): 'Y',
    frozenset('CG'): 'S', frozenset('AT'): 'W',
    frozenset('GT'): 'K', frozenset('AC'): 'M',
}
IUPAC_TO_BASES = {
    'R': ('A', 'G'), 'Y': ('C', 'T'), 'S': ('C', 'G'),
    'W': ('A', 'T'), 'K': ('G', 'T'), 'M': ('A', 'C'),
}


def iupac_encode(b1, b2):
    """返回两个碱基的 IUPAC 歧义编码，不适用时返回 None。"""
    key = frozenset([b1.upper(), b2.upper()])
    return IUPAC_PAIRS.get(key)


def iupac_decode(code):
    """解码 IUPAC 歧义编码，返回两个碱基的元组。"""
    return IUPAC_TO_BASES.get(code.upper())


def is_snp(allele):
    """判断等位基因是否为单核苷酸。"""
    return len(allele) == 1 and allele.upper() in 'ACGT'


# ============================================================
# GT 解析与重编码
# ============================================================

def parse_gt(gt_str):
    """
    解析 GT 字符串，返回 (alleles: tuple, sep: str)。
    - "0/1" → ((0, 1), '/')
    - "1|2" → ((1, 2), '|')
    - "./." → ((None, None), '/')
    - "0"   → ((0,), '/')   (haploid)
    """
    if not gt_str or gt_str == '.':
        return ((None,), '/')
    if '|' in gt_str:
        sep = '|'
    elif '/' in gt_str:
        sep = '/'
    else:
        # 单倍体或未知格式
        return ((None,), '/') if gt_str == '.' else ((int(gt_str),), '/')

    parts = gt_str.split(sep)
    alleles = tuple(int(p) if p != '.' else None for p in parts)
    return (alleles, sep)


def format_gt(alleles, sep):
    """将等位基因元组和分隔符编码为 GT 字符串。
    unphased '/' 按 VCF 惯例小索引在前（0/2 而非 2/0）；
    phased '|' 保留原始顺序（相位有意义）。
    """
    non_missing = [a for a in alleles if a is not None]
    missing_count = sum(1 for a in alleles if a is None)
    if sep == '/':
        non_missing.sort()
    parts = [str(a) for a in non_missing] + ['.'] * missing_count
    return sep.join(parts)


def recode_gt(gt_str, old_to_new):
    """
    用等位基因映射表重新编码 GT 字符串。
    保持分隔符（/ 或 |），缺失值保持不变。
    """
    if not gt_str or gt_str == '.':
        return gt_str
    alleles, sep = parse_gt(gt_str)
    new_alleles = tuple(
        old_to_new.get(a, a) if a is not None else None
        for a in alleles
    )
    return format_gt(new_alleles, sep)


def gt_is_homozygous(gt_str):
    """判断 GT 是否为纯合（两条等位基因相同且非缺失）。"""
    alleles, _ = parse_gt(gt_str)
    if None in alleles:
        return False
    return len(set(alleles)) == 1


def gt_is_het(gt_str):
    """判断 GT 是否为杂合。"""
    alleles, _ = parse_gt(gt_str)
    if None in alleles or len(alleles) < 2:
        return False
    return len(set(alleles)) > 1


def gt_is_missing(gt_str):
    """判断 GT 是否缺失。"""
    alleles, _ = parse_gt(gt_str)
    return all(a is None for a in alleles)


def gt_alleles(gt_str):
    """返回 GT 中的等位基因编号集合（不含缺失）。"""
    alleles, _ = parse_gt(gt_str)
    return {a for a in alleles if a is not None}


# ============================================================
# AD 解析
# ============================================================

def parse_ad(ad_str):
    """
    解析 AD 字符串，返回深度列表。
    "10,5,3" → [10, 5, 3]
    """
    if not ad_str or ad_str == '.':
        return []
    try:
        return [int(x) for x in ad_str.split(',')]
    except (ValueError, TypeError):
        return []


def reorder_list_by_map(lst, old_to_new):
    """
    按映射表重排序列表（用于 AD 重排）。
    new_list[old_to_new[i]] = lst[i]
    """
    n = len(lst)
    new_lst = [0] * n
    for old_idx, new_idx in old_to_new.items():
        if old_idx < n and new_idx < n:
            new_lst[new_idx] = lst[old_idx]
    return new_lst


# ============================================================
# 核心算法
# ============================================================

def build_old_alleles(ref, alt):
    """构建旧等位基因列表：[REF, ALT1, ALT2, ...]"""
    alts = alt.split(',') if alt and alt != '.' else []
    return [ref] + alts


def determine_new_ref_index(sample_gt_str, old_alleles, ad_values, strategy):
    """
    确定新 REF 在 old_alleles 中的索引。

    参数:
        sample_gt_str: 目标样品的 GT 字符串
        old_alleles: [REF, ALT1, ALT2, ...]
        ad_values: AD 深度列表（可能为空）
        strategy: 'strict' | 'iupac' | 'force'

    返回:
        new_ref_idx: 新 REF 在 old_alleles 中的索引
        action: 'swap', 'unchanged', 'iupac', 'skipped'
        iupac_code_or_none: IUPAC 编码（仅 iupac 策略时可能非 None）
    """
    if gt_is_missing(sample_gt_str):
        return 0, 'skipped', None  # 缺失，保持原样

    alleles_set = gt_alleles(sample_gt_str)
    if not alleles_set:
        return 0, 'skipped', None

    # 如果样品携带 REF (0)，且是纯合的
    if alleles_set == {0}:
        return 0, 'unchanged', None  # 纯合 REF，无需交换

    # 纯合非 REF —— 干净交换
    if gt_is_homozygous(sample_gt_str):
        new_ref_idx = list(alleles_set)[0]
        if new_ref_idx == 0:
            return 0, 'unchanged', None
        return new_ref_idx, 'swap', None

    # ===== 杂合位点 =====
    if strategy == 'strict':
        # strict 模式：杂合位点保持不动
        return 0, 'unchanged', None

    elif strategy == 'iupac':
        # iupac 模式：尝试用 IUPAC 编码
        # 首先收集样品携带的碱基
        sample_bases = []
        for idx in sorted(alleles_set):
            if idx < len(old_alleles) and is_snp(old_alleles[idx]):
                sample_bases.append(old_alleles[idx].upper())

        if len(sample_bases) == 2:
            code = iupac_encode(sample_bases[0], sample_bases[1])
            if code:
                # iupac 模式下，REF 变成 IUPAC 编码，但等位基因映射保持不变
                # 返回特殊标记，由上层处理
                return 0, 'iupac', code

        # IUPAC 不适用时（非 SNP 或无法编码），降级为 strict
        return 0, 'unchanged', None

    elif strategy == 'force':
        # force 模式：用 AD 深度选主等位基因
        non_ref_alleles = alleles_set - {0}
        if ad_values and non_ref_alleles:
            # 找出非 REF 等位基因中 AD 最大的
            best_allele = max(non_ref_alleles, key=lambda a: (
                ad_values[a] if a < len(ad_values) else 0
            ))
            # 如果最大 AD 等位基因的 AD > 0，选它
            if best_allele < len(ad_values) and ad_values[best_allele] > 0:
                return best_allele, 'swap', None

        # 降级：选最小的非 REF 等位基因编号
        if non_ref_alleles:
            return min(non_ref_alleles), 'swap', None

        return 0, 'unchanged', None

    return 0, 'unchanged', None


def build_allele_map(old_alleles, new_ref_idx):
    """
    构建 old_idx → new_idx 的映射。

    新等位基因顺序：new_REF 排第一，其余按原顺序（去重）。
    """
    n = len(old_alleles)
    new_alleles = [old_alleles[new_ref_idx]]
    for i, a in enumerate(old_alleles):
        if i != new_ref_idx and a not in new_alleles:
            new_alleles.append(a)

    old_to_new = {}
    for old_idx, allele in enumerate(old_alleles):
        if allele in new_alleles:
            old_to_new[old_idx] = new_alleles.index(allele)
        else:
            old_to_new[old_idx] = old_idx  # fallback，不应发生

    return old_to_new, new_alleles


def new_ref_alt_from_alleles(new_alleles):
    """从新等位基因列表提取 REF 和 ALT。"""
    ref = new_alleles[0]
    alts = new_alleles[1:] if len(new_alleles) > 1 else ['.']
    alt_str = ','.join(alts)
    return ref, alt_str


def generate_original_ref_gt(old_to_new, sep='/'):
    """
    生成 ORIGINAL_REF 样品的 GT。
    旧 REF (index 0) 在新编码中的纯合基因型。
    """
    new_idx = old_to_new.get(0, 0)
    return format_gt((new_idx, new_idx), sep)


# ============================================================
# FORMAT 字段处理
# ============================================================

def process_format_fields(fmt_str, sample_str, old_to_new):
    """
    处理单个样品的 FORMAT 字段。
    返回新的 sample 字符串。
    """
    fmt_keys = fmt_str.split(':')
    values = sample_str.split(':')

    if len(values) != len(fmt_keys):
        # 字段数不匹配，直接返回
        return sample_str

    # 需要 GT 重编码的字段（基因型格式的字段）
    GT_LIKE_FIELDS = {'GT', 'PGT'}

    new_values = []
    for key, val in zip(fmt_keys, values):
        if key in GT_LIKE_FIELDS:
            new_val = recode_gt(val, old_to_new)
        elif key == 'AD':
            # 重排 AD 顺序
            ad = parse_ad(val)
            if ad:
                ad_reordered = reorder_list_by_map(ad, old_to_new)
                new_val = ','.join(str(x) for x in ad_reordered)
            else:
                new_val = val
        elif key == 'PL':
            # PL 丢弃：等位基因重排后 PL 的 genotype combination ordering 不再有效
            # 建议用户用 bcftools +fill-tags 重新计算
            new_val = '.'
        else:
            new_val = val
        new_values.append(new_val)

    return ':'.join(new_values)


# ============================================================
# 主处理流程
# ============================================================

def process_vcf(input_path, output_path, sample_name, strategy,
                new_sample_name, quiet):
    """处理 VCF 文件的主函数。"""
    # 统计计数器
    stats = {'swap': 0, 'unchanged': 0, 'iupac': 0, 'skipped': 0, 'total': 0}

    # 打开输入文件
    if input_path.endswith('.gz') or input_path.endswith('.bgz'):
        infile = gzip.open(input_path, 'rt')
    else:
        infile = open(input_path, 'r')

    # 确定输出
    if output_path:
        if output_path.endswith('.gz'):
            outfile = gzip.open(output_path, 'wt')
        else:
            outfile = open(output_path, 'w')
    else:
        outfile = sys.stdout

    sample_col_idx = None  # 目标样品在 samples 列表中的索引
    header_lines = []
    raw_header_line = None  # #CHROM 行

    try:
        # ========== 第一遍：读取 header ==========
        for line in infile:
            if line.startswith('##'):
                header_lines.append(line.rstrip())
            elif line.startswith('#CHROM'):
                raw_header_line = line.rstrip()
                break
            else:
                # 没有 header 的数据行？不应发生
                break

        if raw_header_line is None:
            print("ERROR: No #CHROM header line found", file=sys.stderr)
            sys.exit(1)

        # 解析 #CHROM 行
        columns = raw_header_line.split('\t')
        if len(columns) < 10:
            print("ERROR: VCF has fewer than 10 columns", file=sys.stderr)
            sys.exit(1)

        fixed_cols = columns[:9]
        sample_names = columns[9:]

        # 查找目标样品
        try:
            sample_col_idx = sample_names.index(sample_name)
        except ValueError:
            print(f"ERROR: Sample '{sample_name}' not found in VCF.", file=sys.stderr)
            print(f"Available samples: {', '.join(sample_names[:20])}...", file=sys.stderr)
            sys.exit(1)

        # ========== 写入输出 header ==========
        for hl in header_lines:
            outfile.write(hl + '\n')

        # 添加 swapRef 命令行记录
        import datetime
        outfile.write(
            f'##swapRef_Version=1.0.0\n'
        )
        outfile.write(
            f'##swapRef_Sample=<ID={sample_name},'
            f'Description="Sample used as new reference genome">\n'
        )
        outfile.write(
            f'##swapRef_Command=<CommandLine="swap_ref.py --sample {sample_name}'
            f' --strategy {strategy}",'
            f'Date="{datetime.datetime.now().strftime("%B %d, %Y at %I:%M:%S %p %Z")}">\n'
        )

        # PL 始终丢弃，更新 FORMAT header
        for hl in header_lines:
            if hl.startswith('##FORMAT=<ID=PL,'):
                outfile.write(
                    '##FORMAT=<ID=PL,Number=.,Type=String,'
                    'Description="PL field dropped by swap_ref — '
                    'no longer valid after allele reordering. '
                    'Regenerate with bcftools +fill-tags">\n'
                )
                break

        # 写入新 #CHROM 行（添加新样品列）
        new_sample_names = list(sample_names) + [new_sample_name]
        outfile.write('\t'.join(fixed_cols + new_sample_names) + '\n')

        # ========== 第二遍：处理数据行 ==========
        for line in infile:
            if line.startswith('#'):
                continue

            parts = line.rstrip('\n\r').split('\t')
            if len(parts) < 10:
                continue

            chrom, pos, vid, ref, alt, qual, filt, info, fmt = parts[:9]
            sample_values = parts[9:]

            # 确保 sample 数量正确（补齐缺失的）
            if len(sample_values) < len(sample_names):
                sample_values += ['.'] * (len(sample_names) - len(sample_values))
            sample_values = sample_values[:len(sample_names)]

            stats['total'] += 1

            # 获取目标样品的 GT & AD
            target_sample_str = sample_values[sample_col_idx]
            target_fmt_values = dict(zip(fmt.split(':'), target_sample_str.split(':')))
            target_gt = target_fmt_values.get('GT', './.')
            target_ad = parse_ad(target_fmt_values.get('AD', ''))

            # 构建旧等位基因列表
            old_alleles = build_old_alleles(ref, alt)

            # 确定新 REF
            new_ref_idx, action, iupac_code = determine_new_ref_index(
                target_gt, old_alleles, target_ad, strategy
            )
            stats[action] += 1

            # 构建等位基因映射
            if action == 'swap' and new_ref_idx != 0:
                old_to_new, new_alleles = build_allele_map(old_alleles, new_ref_idx)
                new_ref, new_alt = new_ref_alt_from_alleles(new_alleles)
            elif action == 'iupac' and iupac_code:
                # IUPAC 模式：REF 变为 IUPAC 编码，ALT 不变（加了 iupac_code）
                # 但严格来说这改变了位点的语义...
                # 更实际的做法：REF 改成 IUPAC 编码，ALT 加上原始两个碱基
                old_to_new = {i: i for i in range(len(old_alleles))}
                new_alleles = list(old_alleles)
                new_ref = iupac_code
                # ALT 保持原样（包含了样品携带的两个碱基作为 ALT）
                # 实际上这里需要更复杂的处理，暂时简化
                new_alt = alt
            else:
                # 不变
                old_to_new = {i: i for i in range(len(old_alleles))}
                new_alleles = list(old_alleles)
                new_ref = ref
                new_alt = alt

            # 处理每个样品
            new_sample_values = []
            for sv in sample_values:
                new_sv = process_format_fields(fmt, sv, old_to_new)
                new_sample_values.append(new_sv)

            # 生成 ORIGINAL_REF 样品
            orig_ref_gt = generate_original_ref_gt(old_to_new, '/')
            # 用原始 FORMAT 字段模板，只填 GT，其余为缺失
            fmt_keys = fmt.split(':')
            orig_values = []
            for k in fmt_keys:
                if k == 'GT':
                    orig_values.append(orig_ref_gt)
                else:
                    orig_values.append('.')
            new_sample_values.append(':'.join(orig_values))

            # 写入输出
            out_parts = [chrom, pos, vid, new_ref, new_alt, qual, filt, info, fmt]
            out_parts.extend(new_sample_values)
            outfile.write('\t'.join(out_parts) + '\n')

    finally:
        infile.close()
        if output_path:
            outfile.close()

    # ========== 输出统计 ==========
    if not quiet:
        print(f"\n{'='*60}", file=sys.stderr)
        print(f"  swap_ref — VCF 参考基因组互换", file=sys.stderr)
        print(f"  新参考样品:      {sample_name}", file=sys.stderr)
        print(f"  策略:            {strategy}", file=sys.stderr)
        print(f"  {'-'*50}", file=sys.stderr)
        print(f"  总位点数:        {stats['total']}", file=sys.stderr)
        print(f"  已翻转 (swap):   {stats['swap']}", file=sys.stderr)
        print(f"  保持不变:        {stats['unchanged']}", file=sys.stderr)
        if stats['iupac'] > 0:
            print(f"  IUPAC 编码:      {stats['iupac']}", file=sys.stderr)
        if stats['skipped'] > 0:
            print(f"  跳过 (缺失):     {stats['skipped']}", file=sys.stderr)
        print(f"  输出样品数:      {len(sample_names) + 1} (含 {new_sample_name})", file=sys.stderr)
        if output_path:
            print(f"  输出文件:        {output_path}", file=sys.stderr)
        print(f"{'='*60}", file=sys.stderr)


# ============================================================
# 命令行接口
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='Reorient 多样品 VCF 的 REF/ALT 到指定样品的等位基因状态',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python swap_ref.py toy.vcf.gz --sample 1_A1
    → 输出: toy_ref-1_A1.vcf.gz (自动命名)
  python swap_ref.py toy.vcf.gz --sample 2_B3 -o custom.vcf.gz

默认行为:
  - strict 模式：仅交换纯合非 REF 位点（最安全）
  - 输出文件名自动生成为 <input>_ref-<SAMPLE>.vcf.gz
  - PL 字段始终丢弃（等位基因重排后需重新计算，可用 bcftools +fill-tags 恢复）
  - --iupac 和 --force 是可选的特殊模式，绝大多数情况不需要
        """
    )

    parser.add_argument('input', help='输入 VCF 文件（支持 .gz 压缩）')
    parser.add_argument('--sample', '-s', required=True,
                        help='作为新参考基因组的样品名（必填）')
    parser.add_argument('--output', '-o', default=None,
                        help='输出 VCF 文件（默认 stdout；.gz 后缀自动压缩）')
    parser.add_argument('--new-sample-name', default='ORIGINAL_REF',
                        help='原始参考列的名称（默认 ORIGINAL_REF）')

    strategy_group = parser.add_mutually_exclusive_group()
    strategy_group.add_argument('--iupac', action='store_true',
                                help='[EXPERIMENTAL] 杂合位点用 IUPAC 歧义碱基编码 REF（仅 SNP；'
                                     '产出非标准 VCF，下游工具可能不兼容）')
    strategy_group.add_argument('--force', action='store_true',
                                help='杂合位点按 AD 深度选主等位基因做新 REF'
                                     '（启发式方法，不推断真实相位）')

    parser.add_argument('--quiet', '-q', action='store_true',
                        help='不输出统计信息到 stderr')

    args = parser.parse_args()

    # 确定策略
    if args.iupac:
        strategy = 'iupac'
    elif args.force:
        strategy = 'force'
    else:
        strategy = 'strict'

    # 如果没有指定输出文件，自动生成：input_ref-SAMPLE.vcf.gz
    if args.output is None:
        import os
        base = os.path.basename(args.input)
        # 去掉 .gz 后缀
        if base.endswith('.gz'):
            base = base[:-3]
        if base.endswith('.vcf'):
            base = base[:-4]
        args.output = f'{base}_ref-{args.sample}.vcf.gz'

    process_vcf(
        input_path=args.input,
        output_path=args.output,
        sample_name=args.sample,
        strategy=strategy,
        new_sample_name=args.new_sample_name,
        quiet=args.quiet,
    )


# ============================================================
# 自测
# ============================================================

def self_test():
    """内置单元测试，验证核心逻辑。"""
    import copy

    passed = 0
    failed = 0

    def check(name, actual, expected):
        nonlocal passed, failed
        if actual == expected:
            passed += 1
            print(f"  [PASS] {name}")
        else:
            failed += 1
            print(f"  [FAIL] {name}: got {actual!r}, expected {expected!r}")

    print("=== swap_ref self-test ===\n")

    # --- parse_gt ---
    print("[parse_gt]")
    check("biallelic ref", parse_gt("0/0"), ((0, 0), '/'))
    check("biallelic het", parse_gt("0/1"), ((0, 1), '/'))
    check("biallelic alt", parse_gt("1/1"), ((1, 1), '/'))
    check("phased", parse_gt("0|1"), ((0, 1), '|'))
    check("missing", parse_gt("./."), ((None, None), '/'))
    check("multi-allelic", parse_gt("1/2"), ((1, 2), '/'))
    check("triploid", parse_gt("0/0/1"), ((0, 0, 1), '/'))
    check("empty", parse_gt("."), ((None,), '/'))
    check("empty_str", parse_gt(""), ((None,), '/'))

    # --- format_gt ---
    print("\n[format_gt]")
    check("simple", format_gt((0, 1), '/'), "0/1")
    check("phased", format_gt((1, 0), '|'), "1|0")
    check("missing", format_gt((None, None), '/'), "./.")

    # --- recode_gt ---
    print("\n[recode_gt]")
    # 0→1, 1→0 (simple swap)
    swap_map = {0: 1, 1: 0}
    check("0/0→1/1", recode_gt("0/0", swap_map), "1/1")
    check("0/1→0/1 (sorted)", recode_gt("0/1", swap_map), "0/1")
    check("1/1→0/0", recode_gt("1/1", swap_map), "0/0")
    check("./. unchanged", recode_gt("./.", swap_map), "./.")
    check("0|1→1|0", recode_gt("0|1", swap_map), "1|0")

    # 0→2, 1→0, 2→1 (cyclic)
    cyclic_map = {0: 2, 1: 0, 2: 1}
    check("0/0→2/2", recode_gt("0/0", cyclic_map), "2/2")
    check("1/2→0/1", recode_gt("1/2", cyclic_map), "0/1")

    # --- gt is_xxx ---
    print("\n[gt is_xxx]")
    check("0/0 homo", gt_is_homozygous("0/0"), True)
    check("0/1 not homo", gt_is_homozygous("0/1"), False)
    check("./. not homo", gt_is_homozygous("./."), False)
    check("0/1 het", gt_is_het("0/1"), True)
    check("0/0 not het", gt_is_het("0/0"), False)
    check("./. missing", gt_is_missing("./."), True)
    check("0/0 not missing", gt_is_missing("0/0"), False)

    # --- determine_new_ref_index ---
    print("\n[determine_new_ref_index — strict mode]")
    old_alleles_snp = ['A', 'G']
    old_alleles_multi = ['G', 'T', 'A']

    # 纯合 REF
    idx, action, code = determine_new_ref_index("0/0", old_alleles_snp, [], 'strict')
    check("00→unchanged", (idx, action), (0, 'unchanged'))

    # 纯合 ALT
    idx, action, code = determine_new_ref_index("1/1", old_alleles_snp, [], 'strict')
    check("11→swap idx=1", (idx, action), (1, 'swap'))

    # 杂合 — strict 模式不动
    idx, action, code = determine_new_ref_index("0/1", old_alleles_snp, [10, 3], 'strict')
    check("01 strict→unchanged", (idx, action), (0, 'unchanged'))

    # 缺失
    idx, action, code = determine_new_ref_index("./.", old_alleles_snp, [], 'strict')
    check("./. → skipped", (idx, action), (0, 'skipped'))

    # 多等位基因位点 2/2
    idx, action, code = determine_new_ref_index("2/2", old_alleles_multi, [], 'strict')
    check("22→swap idx=2", (idx, action), (2, 'swap'))

    # --- determine_new_ref_index — force mode ---
    print("\n[determine_new_ref_index — force mode]")
    # 杂合 0/1, AD: REF=3, ALT=10 → 应选 ALT
    idx, action, code = determine_new_ref_index("0/1", old_alleles_snp, [3, 10], 'force')
    check("01 AD[3,10]→swap idx=1", (idx, action), (1, 'swap'))

    # 杂合 0/1, AD: REF=10, ALT=3 → 非 REF 只有 1，选它
    idx, action, code = determine_new_ref_index("0/1", old_alleles_snp, [10, 3], 'force')
    check("01 AD[10,3]→swap idx=1", (idx, action), (1, 'swap'))

    # 杂合 1/2, AD: T=8, A=2 → 选 1
    idx, action, code = determine_new_ref_index("1/2", old_alleles_multi, [5, 8, 2], 'force')
    check("12 AD[5,8,2]→swap idx=1", (idx, action), (1, 'swap'))

    # --- determine_new_ref_index — iupac mode ---
    print("\n[determine_new_ref_index — iupac mode]")
    # A/G 杂合 → IUPAC R
    idx, action, code = determine_new_ref_index("0/1", ['A', 'G'], [], 'iupac')
    check("01 A/G→iupac R", (idx, action, code), (0, 'iupac', 'R'))

    # C/T 杂合 → IUPAC Y
    idx, action, code = determine_new_ref_index("0/1", ['C', 'T'], [], 'iupac')
    check("01 C/T→iupac Y", (idx, action, code), (0, 'iupac', 'Y'))

    # 非 SNP（INDEL）→ 降级为 strict
    idx, action, code = determine_new_ref_index("0/1", ['C', 'CT'], [], 'iupac')
    check("01 indel→unchanged", (idx, action, code), (0, 'unchanged', None))

    # --- build_allele_map ---
    print("\n[build_allele_map]")
    old_map, new_alleles = build_allele_map(['A', 'G'], 1)  # ALT 变 REF
    check("swap A/G→G,A", (old_map, new_alleles), ({0: 1, 1: 0}, ['G', 'A']))

    old_map, new_alleles = build_allele_map(['G', 'T', 'A'], 2)  # 第二 ALT 变 REF
    check("multi: 2→REF", new_alleles, ['A', 'G', 'T'])
    check("multi: map", old_map, {0: 1, 1: 2, 2: 0})

    old_map, new_alleles = build_allele_map(['A', 'G'], 0)  # REF 不变
    check("no swap", (old_map, new_alleles), ({0: 0, 1: 1}, ['A', 'G']))

    # --- generate_original_ref_gt ---
    print("\n[generate_original_ref_gt]")
    # 发生了交换 (0→1)，旧 REF 在新编码中为 1
    check("old ref→1/1", generate_original_ref_gt({0: 1, 1: 0}), "1/1")
    # 没交换
    check("old ref→0/0", generate_original_ref_gt({0: 0, 1: 1}), "0/0")
    # 多等位基因：旧 REF→2
    check("old ref→2/2", generate_original_ref_gt({0: 2, 1: 0, 2: 1}), "2/2")

    # --- iupac ---
    print("\n[iupac]")
    check("A/G→R", iupac_encode('A', 'G'), 'R')
    check("C/T→Y", iupac_encode('C', 'T'), 'Y')
    check("C/G→S", iupac_encode('C', 'G'), 'S')
    check("A/T→W", iupac_encode('A', 'T'), 'W')
    check("G/T→K", iupac_encode('G', 'T'), 'K')
    check("A/C→M", iupac_encode('A', 'C'), 'M')
    check("same base→None", iupac_encode('A', 'A'), None)

    # --- parse_ad ---
    print("\n[parse_ad]")
    check("normal", parse_ad("10,5,3"), [10, 5, 3])
    check("empty", parse_ad(""), [])
    check("dot", parse_ad("."), [])
    check("single", parse_ad("42"), [42])

    # --- reorder_list_by_map ---
    print("\n[reorder_list_by_map]")
    check("swap [10,5]→[5,10]",
          reorder_list_by_map([10, 5], {0: 1, 1: 0}), [5, 10])
    check("cyclic [5,8,2]→[2,5,8]",
          reorder_list_by_map([5, 8, 2], {0: 1, 1: 2, 2: 0}), [2, 5, 8])

    # --- 结果 ---
    total = passed + failed
    print(f"\n{'='*40}")
    print(f"Results: {passed}/{total} passed", end='')
    if failed > 0:
        print(f", {failed} FAILED")
    else:
        print(" -- ALL PASS")
    print(f"{'='*40}")

    return failed == 0


if __name__ == '__main__':
    if '--test' in sys.argv or '--self-test' in sys.argv:
        ok = self_test()
        sys.exit(0 if ok else 1)
    else:
        main()
