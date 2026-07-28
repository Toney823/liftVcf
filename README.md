# swap_ref — VCF Reference Swapping Tool

**One-liner**: pick any sample in a VCF as the new reference genome, and turn the original reference into a sample column.

> [中文版](README_CN.md)

## Quick Start

```bash
# 1. List sample names
python3 -c "import gzip; f=gzip.open('toy.vcf.gz','rt');
[print(l.split()[9]) for l in f if l.startswith('#CHROM')]; break"

# 2. Pick a sample and run (no extra flags needed)
python3 swap_ref.py toy.vcf.gz --sample 1_H3

# 3. Verify correctness
python3 verify_swap.py toy.vcf.gz --swap --sample 1_H3
```

## What It Does

Say you pick sample `1_H3`. At each variant site:

| 1_H3 genotype | Meaning | Result |
|---|---|---|
| `0/0` | Same as reference | Unchanged |
| `1/1` | Homozygous ALT | **Swap**: ALT becomes new REF, old REF becomes ALT |
| `0/1` | Heterozygous (both alleles) | **Unchanged**: sample is uncertain, can't define a single "reference" |
| `./.` | Missing | Skipped |

After swapping, the VCF gains one extra column: `ORIGINAL_REF` — the original reference genome's genotypes.

## No Flags Needed

For most use cases, the defaults are all you need:

```bash
python3 swap_ref.py input.vcf.gz --sample SAMPLE
```

## What Are `--iupac` and `--force`? (Optional)

These only affect **heterozygous sites** (where the sample's genotype is `0/1`).

By default, heterozygous sites are left alone — the sample carries both alleles, so there's no single "reference" to define. If you insist on changing the reference even at these sites, you must decide: **which allele should become the new REF?**

| Flag | Strategy | Example |
|------|----------|---------|
| (default, no flag) | Leave heterozygous sites unchanged | REF=A, ALT=G, sample 0/1 → REF stays A |
| `--iupac` | Use IUPAC ambiguity code for both alleles | sample 0/1 (A/G) → REF becomes `R` (=A or G) |
| `--force` | Use allele depth: the allele with more reads becomes REF | sample 0/1, AD=3,15 → ALT has higher depth → REF becomes G |

**When to use them?**
- `--iupac`: you want the REF column to honestly reflect "this sample is polymorphic here". Note: most downstream tools (bwa, GATK) don't handle IUPAC bases — use with caution.
- `--force`: you don't care about heterozygosity and just want to maximize the number of swapped sites.
- Neither: safe default, works for virtually all use cases.

## Output

- Auto-named: `toy_ref-1_H3.vcf.gz` (`<input>_ref-<sample>.vcf.gz`)
- One extra column `ORIGINAL_REF` — the original reference genome's genotypes
- VCF header records the sample used: `##swapRef_Sample=<ID=1_H3,...>`

## Verification

Allele re-encoding is a bijection — allele-sharing between samples is preserved:

```
IBS distance matrix (original)  ==  IBS distance matrix (swapped)
              ↓                              ↓
         clustering tree               clustering tree
              ↓                              ↓
          topology          ===          topology
```

Same distances → same tree. `verify_swap.py` runs this check automatically.

## All Options

```
python3 swap_ref.py input.vcf.gz \
    --sample SAMPLE        # required: sample to use as new reference
    -o output.vcf.gz       # output file (auto-named if omitted)
    --iupac                # optional: IUPAC ambiguity codes at heterozygous sites
    --force                # optional: use AD to pick major allele at heterozygous sites
    --new-sample-name NAME # name for original-reference column (default: ORIGINAL_REF)
    --keep-pl              # preserve PL fields
    -q                     # quiet mode
```
