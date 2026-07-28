# swap_ref — VCF Reference Swapping Tool

**One-liner**: pick any sample in a VCF as the new reference genome, and turn the original reference into a sample column.

> [中文版](README_CN.md)

## Quick Start

```bash
# 1. List sample names
python3 -c "import gzip; f=gzip.open('toy.vcf.gz','rt');
[print(l.split()[9]) for l in f if l.startswith('#CHROM')]; break"

# 2. Pick a sample and run
python3 swap_ref.py toy.vcf.gz --sample 1_H3

# 3. Verify correctness
python3 verify_swap.py toy.vcf.gz --swap --sample 1_H3
```

## Three Modes

| Mode | Command | Behavior |
|------|---------|----------|
| **strict** (default·recommended) | `--sample 1_H3` | Only swap homozygous non-REF sites; hets unchanged |
| **iupac** | `--sample 1_H3 --iupac` | Hets encoded as IUPAC ambiguity codes (A/C→M, G/T→K...) |
| **force** | `--sample 1_H3 --force` | Hets use allele depth to pick the major allele as new REF |

## Output

- Auto-named: `toy_ref-1_H3.vcf.gz` (`<input>_ref-<sample>.vcf.gz`)
- One extra column `ORIGINAL_REF` — the original reference genome's genotypes
- VCF header records the sample used: `##swapRef_Sample=<ID=1_H3,...>`

## Verification

The allele re-encoding is a bijection, so allele-sharing between samples is preserved:

```
IBS distance matrix (original)  ==  IBS distance matrix (swapped)
              ↓                              ↓
         clustering tree               clustering tree
              ↓                              ↓
          topology          ===          topology
```

## All Options

```
python3 swap_ref.py input.vcf.gz \
    --sample SAMPLE        # required: sample to use as new reference
    -o output.vcf.gz       # output file (auto-named if omitted)
    --iupac                # IUPAC encoding at heterozygous sites (SNP only)
    --force                # use AD to pick major allele at heterozygous sites
    --new-sample-name NAME # name for original-reference column (default: ORIGINAL_REF)
    --keep-pl              # preserve PL fields
    -q                     # quiet mode
```
