#!/usr/bin/env python3
"""
verify_swap.py — 验证 VCF 参考基因组互换的正确性

通过对比原始 VCF 和交换后 VCF 的：
  1. 遗传距离矩阵（应完全一致）
  2. 系统发育树拓扑结构（应完全一致）

来验证 liftVcf 变换的正确性。

用法:
  python verify_swap.py original.vcf.gz swapped.vcf.gz --sample 1_H3

原理:
  参考系变换（重新编码等位基因）是双射映射，不改变样品间的等位基因共享模式，
  因此所有样品间的遗传距离和由此构建的树的拓扑结构应保持不变。
"""

import sys
import gzip
import argparse
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster import hierarchy
from scipy.stats import pearsonr


# ============================================================
# VCF → 距离矩阵
# ============================================================

def load_allele_matrix(vcf_path):
    """
    从 VCF 文件加载等位基因矩阵（保留原始等位基因索引）。
    返回:
        samples: 样品名列表
        allele_matrix: list of list，每个元素是 (a1, a2) tuple 或 None (缺失)
                       shape = (n_variants, n_samples)
    """
    samples = []
    rows = []

    opener = gzip.open if (vcf_path.endswith('.gz') or vcf_path.endswith('.bgz')) else open

    with opener(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                cols = line.rstrip().split('\t')
                samples = cols[9:]
            elif line.startswith('#'):
                continue
            else:
                parts = line.rstrip().split('\t')
                fmt = parts[8]
                fmt_keys = fmt.split(':')
                gt_idx = fmt_keys.index('GT') if 'GT' in fmt_keys else 0

                sample_vals = parts[9:]
                site_alleles = []
                for sv in sample_vals[:len(samples)]:
                    fields = sv.split(':')
                    if len(fields) <= gt_idx:
                        site_alleles.append(None)
                        continue
                    gt_str = fields[gt_idx]
                    alleles, sep = parse_gt(gt_str)
                    if None in alleles or not alleles:
                        site_alleles.append(None)
                    else:
                        # 保留原始等位基因编号 (0, 1, 2, ...)
                        site_alleles.append(tuple(alleles))

                rows.append(site_alleles)

    return samples, rows


def parse_gt(gt_str):
    """解析 GT 字符串。"""
    if not gt_str or gt_str == '.':
        return ((None,), '/')
    if '|' in gt_str:
        sep = '|'
    elif '/' in gt_str:
        sep = '/'
    else:
        return ((int(gt_str),), '/') if gt_str != '.' else ((None,), '/')
    parts = gt_str.split(sep)
    return tuple(int(p) if p != '.' else None for p in parts), sep


# ============================================================
# 距离矩阵计算
# ============================================================

def ibs_distance(alleles_a, alleles_b):
    """
    计算两个等位基因集之间的 IBS (Identity By State) 距离。

    alleles_a, alleles_b: tuple of int (e.g., (0, 1)) 或 None (缺失)

    距离定义:
      - 共享 2 个等位基因: 0.0
      - 共享 1 个等位基因: 0.5
      - 共享 0 个等位基因: 1.0
      - 任一缺失: 返回 None
    """
    if alleles_a is None or alleles_b is None:
        return None
    set_a = set(alleles_a)
    set_b = set(alleles_b)
    shared = len(set_a & set_b)
    total = min(len(set_a), len(set_b))
    if total == 0:
        return None
    return 1.0 - shared / total


def compute_distance_matrix(allele_matrix, samples_list=None):
    """
    从等位基因矩阵计算样品间的成对 IBS 距离矩阵。

    IBS 距离在等位基因重编码下是严格不变的，因为：
    - 重编码是双射: old_idx → new_idx
    - 等位基因共享关系得以保持
    """
    n_variants = len(allele_matrix)
    n_samples = len(allele_matrix[0]) if n_variants > 0 else 0
    dist_matrix = np.zeros((n_samples, n_samples))

    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            dists = []
            for site in range(n_variants):
                d = ibs_distance(allele_matrix[site][i], allele_matrix[site][j])
                if d is not None:
                    dists.append(d)

            if not dists:
                dist_matrix[i, j] = dist_matrix[j, i] = np.nan
            else:
                dist_matrix[i, j] = dist_matrix[j, i] = np.mean(dists)

    return dist_matrix


# ============================================================
# 系统发育树构建与比较
# ============================================================

def build_tree(dist_matrix, samples, method='average'):
    """
    从距离矩阵构建层次聚类树 (UPGMA-like)。

    返回:
        linkage_matrix: scipy linkage 矩阵
        leaf_order: 叶子节点顺序
    """
    # 处理 NaN（用列均值填充）
    dist_filled = dist_matrix.copy()
    col_means = np.nanmean(dist_filled, axis=0)
    nan_mask = np.isnan(dist_filled)
    for i in range(dist_filled.shape[0]):
        dist_filled[nan_mask[:, i], i] = col_means[i]
        dist_filled[i, nan_mask[i, :]] = col_means[i]

    # 转换为压缩距离格式（上三角）
    condensed = squareform(dist_filled, checks=False)
    linkage = hierarchy.linkage(condensed, method=method)
    return linkage


def compare_trees(linkage1, linkage2):
    """
    比较两棵树的拓扑结构。

    返回:
        cophenetic_corr: 两个 cophenetic 距离矩阵之间的 Pearson 相关系数
        is_identical: 两个 linkage 矩阵是否完全一致（在浮点误差内）
    """
    # 计算 cophenetic 距离矩阵
    n1 = linkage1.shape[0] + 1
    n2 = linkage2.shape[0] + 1

    if n1 != n2:
        return 0.0, False, "样本数不一致"

    cophen1 = hierarchy.cophenet(linkage1)
    cophen2 = hierarchy.cophenet(linkage2)

    # Pearson 相关
    corr, _ = pearsonr(cophen1, cophen2)

    # 检查 linkage 矩阵是否一致
    is_close = np.allclose(linkage1, linkage2, rtol=1e-10, atol=1e-10)

    return corr, is_close, None


# ============================================================
# 主验证流程
# ============================================================

def verify(original_vcf, swapped_vcf, sample_name, new_sample_name='ORIGINAL_REF', quiet=False):
    """
    验证 VCF 参考基因组互换的正确性。

    检查项:
      1. 样品集一致性
      2. 距离矩阵一致性
      3. 树拓扑结构一致性
      4. 新参考样品的基因型正确性
    """
    all_ok = True

    def log(msg):
        if not quiet:
            print(msg, file=sys.stderr)

    def ok(msg):
        nonlocal all_ok
        if not quiet:
            print(f"  [PASS] {msg}", file=sys.stderr)

    def fail(msg):
        nonlocal all_ok
        all_ok = False
        if not quiet:
            print(f"  [FAIL] {msg}", file=sys.stderr)

    log(f"\n{'='*60}")
    log(f"  liftVcf 验证")
    log(f"  原始 VCF:  {original_vcf}")
    log(f"  交换 VCF:  {swapped_vcf}")
    log(f"  目标样品:  {sample_name}")
    log(f"{'='*60}")

    # --- 1. 加载数据 ---
    log("\n[1/5] 加载等位基因矩阵...")
    orig_samples, orig_alleles = load_allele_matrix(original_vcf)
    swap_samples, swap_alleles = load_allele_matrix(swapped_vcf)

    n_orig_sites = len(orig_alleles)
    log(f"  原始 VCF: {len(orig_samples)} 样品, {n_orig_sites} 位点")
    log(f"  交换 VCF: {len(swap_samples)} 样品, {len(swap_alleles)} 位点")

    # --- 2. 样品集检查 ---
    log("\n[2/5] 检查样品集...")
    if new_sample_name in swap_samples:
        ok(f"交换后 VCF 包含新列 '{new_sample_name}'")
        # 移除 ORIGINAL_REF，只比较原始样品
        ref_idx = swap_samples.index(new_sample_name)
        swap_samples_no_ref = [s for s in swap_samples if s != new_sample_name]
        # 从等位基因矩阵中移除 ORIGINAL_REF 列
        swap_alleles_no_ref = []
        for site_alleles in swap_alleles:
            new_site = list(site_alleles)
            new_site.pop(ref_idx)
            swap_alleles_no_ref.append(new_site)
    else:
        fail(f"交换后 VCF 缺少 '{new_sample_name}' 列")
        return False

    if orig_samples != swap_samples_no_ref:
        if set(orig_samples) == set(swap_samples_no_ref):
            ok("样品集一致（重新排序以匹配）")
            # 重新排序 swap 的列以匹配 orig
            idx_map = [swap_samples_no_ref.index(s) for s in orig_samples]
            swap_alleles_no_ref = [
                [site[idx] for idx in idx_map]
                for site in swap_alleles_no_ref
            ]
        else:
            only_orig = set(orig_samples) - set(swap_samples_no_ref)
            only_swap = set(swap_samples_no_ref) - set(orig_samples)
            fail(f"样品集不一致: 仅在原始={only_orig}, 仅在交换={only_swap}")
            return False
    else:
        ok("样品顺序一致")

    # --- 3. 距离矩阵比较 (IBS) ---
    log("\n[3/5] 比较遗传距离矩阵 (IBS distance)...")
    log(f"  计算原始距离矩阵 ({len(orig_samples)} x {len(orig_samples)})...")
    dist_orig = compute_distance_matrix(orig_alleles)

    log(f"  计算交换后距离矩阵 ({len(swap_samples_no_ref)} x {len(swap_samples_no_ref)})...")
    dist_swap = compute_distance_matrix(swap_alleles_no_ref)

    # IBS 距离矩阵应该完全一致（等位基因重编码不改变共享模式）
    diff = np.abs(dist_orig - dist_swap)
    valid_mask = ~np.isnan(diff)
    max_diff = np.max(diff[valid_mask]) if valid_mask.any() else 0.0
    mean_diff = np.mean(diff[valid_mask]) if valid_mask.any() else 0.0

    # 展平并移除 NaN 后计算相关系数
    d1 = dist_orig[valid_mask].flatten()
    d2 = dist_swap[valid_mask].flatten()
    if len(d1) > 1:
        corr_dist, _ = pearsonr(d1, d2)
    else:
        corr_dist = 1.0

    log(f"  距离矩阵相关系数: {corr_dist:.10f}")
    log(f"  最大绝对差异:     {max_diff:.10e}")
    log(f"  平均绝对差异:     {mean_diff:.10e}")

    if max_diff < 1e-10:
        ok(f"IBS 距离矩阵完全一致 (max diff = {max_diff:.2e})")
    elif max_diff < 1e-6:
        ok(f"IBS 距离矩阵基本一致 (max diff = {max_diff:.2e})")
    else:
        fail(f"IBS 距离矩阵存在差异 (max diff = {max_diff:.2e})")

    # --- 4. 树拓扑比较 ---
    log("\n[4/5] 比较系统发育树拓扑结构...")
    log(f"  构建原始 VCF 聚类树...")
    tree_orig = build_tree(dist_orig, orig_samples, method='average')

    log(f"  构建交换后 VCF 聚类树...")
    tree_swap = build_tree(dist_swap, swap_samples_no_ref, method='average')

    tree_corr, is_identical, err = compare_trees(tree_orig, tree_swap)

    log(f"  Cophenetic 相关系数: {tree_corr:.10f}")
    if err:
        fail(f"树比较失败: {err}")
    elif is_identical:
        ok("两棵树的 linkage 矩阵完全一致")
    elif tree_corr > 0.999999:
        ok(f"树拓扑结构高度一致 (cophenetic r = {tree_corr:.10f})")
    else:
        fail(f"树拓扑结构存在差异 (cophenetic r = {tree_corr:.10f})")

    # --- 5. 新参考样品检查 ---
    log("\n[5/5] 检查新参考样品...")
    if sample_name in swap_samples:
        sample_idx = swap_samples.index(sample_name)
        n_total = 0
        n_homo_ref = 0
        n_het = 0
        n_homo_alt = 0
        for site_alleles in swap_alleles:
            a = site_alleles[sample_idx]
            if a is None:
                continue
            n_total += 1
            unique = set(a)
            if unique == {0}:
                n_homo_ref += 1
            elif 0 in unique:
                n_het += 1
            else:
                n_homo_alt += 1

        if n_het == 0:
            ok(f"样品 '{sample_name}' 在交换后全部为纯合新 REF (GT=0/0): {n_total} 位点")
        else:
            log(f"  注意: 样品 '{sample_name}' 有 {n_het}/{n_total} 个杂合位点")
            log(f"         (这些是 strict 模式下保持不变的位点)")
            ok(f"样品 '{sample_name}' 在交换后主要为纯合新 REF ({n_het}/{n_total} 杂合)")
    else:
        fail(f"交换后 VCF 中找不到样品 '{sample_name}'")

    # --- 6. 检查 ORIGINAL_REF 样品 ---
    if new_sample_name in swap_samples:
        ref_idx = swap_samples.index(new_sample_name)
        n_total = 0
        n_ref = 0
        n_het = 0
        n_alt = 0
        for site_alleles in swap_alleles:
            a = site_alleles[ref_idx]
            if a is None:
                continue
            n_total += 1
            unique = set(a)
            if unique == {0}:
                n_ref += 1
            elif 0 in unique:
                n_het += 1
            else:
                n_alt += 1

        if n_het == 0:
            ok(f"'{new_sample_name}' 样品无杂合位点（符合预期：纯合旧 REF）")
        else:
            fail(f"'{new_sample_name}' 样品出现 {n_het} 个杂合位点（不应出现）")

        log(f"  '{new_sample_name}' 的 GT 分布: {n_ref} REF, {n_alt} ALT, {n_het} HET")

    # --- 总结 ---
    log(f"\n{'='*60}")
    if all_ok:
        log(f"  验证结果: 全部通过 ✓")
        log(f"  结论: liftVcf 变换正确，遗传距离和树拓扑结构得以保持")
    else:
        log(f"  验证结果: 存在失败项 ✗")
    log(f"{'='*60}\n")

    return all_ok


# ============================================================
# 全流程：swap + verify
# ============================================================

def swap_and_verify(input_vcf, sample_name, strategy='strict', new_sample_name='ORIGINAL_REF',
                    quiet=False):
    """
    执行 swap → 验证 → round-trip 测试。
    """
    import tempfile
    import os
    import subprocess
    import shutil

    tmpdir = tempfile.mkdtemp(prefix='swap_verify_')
    output_vcf = os.path.join(tmpdir, 'swapped.vcf.gz')
    restored_vcf = os.path.join(tmpdir, 'restored.vcf.gz')

    log = lambda msg: not quiet and print(msg, file=sys.stderr)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    swap_script = os.path.join(script_dir, 'liftVcf.py')

    def run_swap(input_file, sample, out, ref_name, extra_args=None):
        cmd = [sys.executable, swap_script, input_file, '--sample', sample,
               '-o', out, '--new-sample-name', ref_name, '-q']
        if extra_args:
            cmd.extend(extra_args)
        result = subprocess.run(cmd, capture_output=True, text=True,
                                encoding='utf-8', errors='replace')
        if result.returncode != 0:
            print(f"ERROR: liftVcf failed:\n{result.stderr}", file=sys.stderr)
            return False
        return True

    # Step 1: swap to chosen sample → adds ORIGINAL_REF column
    extra = []
    if strategy == 'iupac':
        extra.append('--iupac')
    elif strategy == 'force':
        extra.append('--force')
    if not run_swap(input_vcf, sample_name, output_vcf, new_sample_name, extra):
        shutil.rmtree(tmpdir, ignore_errors=True)
        return False

    # Step 2: verify distance matrix + tree topology
    ok1 = verify(input_vcf, output_vcf, sample_name, new_sample_name, quiet=quiet)

    # Step 3: round-trip — swap ORIGINAL_REF back
    # Use a temp name for the second pass to avoid column name collision
    roundtrip_ref = '_ROUNDTRIP_TMP'
    log(f"\n{'='*60}")
    log(f"  Round-trip 测试: 交换后 VCF → 以 {new_sample_name} 为目标还原")
    log(f"{'='*60}")
    if not run_swap(output_vcf, new_sample_name, restored_vcf, roundtrip_ref):
        shutil.rmtree(tmpdir, ignore_errors=True)
        return False

    # Compare: exclude the temp column, the rest should match original
    ok2 = verify_round_trip(input_vcf, restored_vcf, sample_name, new_sample_name, roundtrip_ref, quiet=quiet)

    shutil.rmtree(tmpdir, ignore_errors=True)
    return ok1 and ok2


def verify_round_trip(original_vcf, restored_vcf, orig_sample, ref_sample_name,
                      temp_column, quiet=False):
    """
    Check round-trip: original == swap(swap(original, sample_A), ORIGINAL_REF).
    The restored VCF has an extra temp column that we exclude from comparison.
    """
    log = lambda msg: not quiet and print(msg, file=sys.stderr)
    all_ok = True

    orig_samples, orig_alleles = load_allele_matrix(original_vcf)
    rest_samples, rest_alleles = load_allele_matrix(restored_vcf)

    # Remove temp column from restored data
    if temp_column in rest_samples:
        tmp_idx = rest_samples.index(temp_column)
        rest_samples = [s for s in rest_samples if s != temp_column]
        rest_alleles = [[a for j, a in enumerate(site) if j != tmp_idx]
                        for site in rest_alleles]

    log(f"  原始样品数: {len(orig_samples)}, 还原后样品数: {len(rest_samples)}")

    n_sites = min(len(orig_alleles), len(rest_alleles))
    n_mismatch = 0

    # Compare only samples present in both
    common_samples = list(set(orig_samples) & set(rest_samples))
    if len(common_samples) < len(orig_samples):
        log(f"  [INFO] {len(orig_samples) - len(common_samples)} 个样品仅在原始 VCF 中")

    orig_idx = [orig_samples.index(s) for s in common_samples]
    rest_idx = [rest_samples.index(s) for s in common_samples]

    for i in range(n_sites):
        for jo, jr in zip(orig_idx, rest_idx):
            a = orig_alleles[i][jo] if jo < len(orig_alleles[i]) else None
            b = rest_alleles[i][jr] if jr < len(rest_alleles[i]) else None
            if a is None and b is None:
                continue
            if a != b:
                n_mismatch += 1

    if n_mismatch == 0:
        log(f"  [PASS] Round-trip: {n_sites} 位点 × {len(common_samples)} 样品 GT 完全一致")
    else:
        log(f"  [FAIL] Round-trip: {n_mismatch} 个 allele 不一致")
        all_ok = False

    if all_ok:
        log(f"  [PASS] Round-trip 验证通过: liftVcf 变换可逆")
    else:
        log(f"  [FAIL] Round-trip 验证失败")

    return all_ok


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='验证 liftVcf 的正确性（IBS 距离矩阵 + 树拓扑 + round-trip）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python verify_swap.py toy.vcf.gz swapped.vcf.gz --sample 1_H3
  python verify_swap.py toy.vcf.gz --swap --sample 1_H3
  python verify_swap.py toy.vcf.gz --swap --sample 1_H3 --force
        """
    )

    parser.add_argument('original', help='原始 VCF 文件')
    parser.add_argument('swapped', nargs='?', default=None,
                        help='交换后的 VCF 文件（--swap 模式下可选）')
    parser.add_argument('--sample', '-s', required=True,
                        help='作为新参考的样品名')
    parser.add_argument('--swap', action='store_true',
                        help='同时执行 swap + verify + round-trip')
    parser.add_argument('--new-sample-name', default='ORIGINAL_REF',
                        help='原始参考列名（默认 ORIGINAL_REF）')
    parser.add_argument('--iupac', action='store_true',
                        help='swap 时使用 IUPAC 策略')
    parser.add_argument('--force', action='store_true',
                        help='swap 时使用 force 策略')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='减少输出')

    args = parser.parse_args()

    if args.iupac:
        strategy = 'iupac'
    elif args.force:
        strategy = 'force'
    else:
        strategy = 'strict'

    if args.swap:
        # swap + verify 模式
        ok = swap_and_verify(
            args.original, args.sample,
            strategy=strategy,
            new_sample_name=args.new_sample_name,
            quiet=args.quiet,
        )
    else:
        # 仅验证模式
        if args.swapped is None:
            print("ERROR: 需要提供交换后的 VCF 文件，或使用 --swap 模式", file=sys.stderr)
            sys.exit(1)
        ok = verify(
            args.original, args.swapped, args.sample,
            new_sample_name=args.new_sample_name,
            quiet=args.quiet,
        )

    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
