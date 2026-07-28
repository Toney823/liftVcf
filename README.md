# swap_ref — Reorient a Multi-Sample VCF to a Selected Sample

**Reorient VCF REF/ALT alleles using a selected sample's genotype as the new reference allele state.** The original reference genome is preserved as an extra sample column.

> [中文版](README_CN.md)

## Quick Start

```bash
# 1. List sample names
python3 -c "import gzip; f=gzip.open('toy.vcf.gz','rt');
[print(l.split()[9]) for l in f if l.startswith('#CHROM')]; break"

# 2. Pick a sample and run (no extra flags needed)
python3 swap_ref.py toy.vcf.gz --sample 1_H3

# 3. Verify correctness (distance matrix + tree topology + round-trip)
python3 verify_swap.py toy.vcf.gz --swap --sample 1_H3
```

## What It Does

Say you pick sample `1_H3`. At each variant site:

| 1_H3 genotype | Meaning | Result |
|---|---|---|
| `0/0` | Same as reference | Unchanged |
| `1/1` | Homozygous ALT | **Swap**: ALT → new REF, old REF → new ALT |
| `0/1` | Heterozygous | Unchanged (sample carries both alleles) |
| `./.` | Missing | Skipped |

After reorientation, the VCF gains one extra column: `ORIGINAL_REF` — the original reference genome as a pseudo-sample.

## No Flags Needed

```bash
python3 swap_ref.py input.vcf.gz --sample SAMPLE
```

## Optional Flags for Heterozygous Sites

By default, heterozygous sites (`0/1`) are left unchanged — the sample carries both alleles, so there's no single allele to define as reference. These flags override that behavior:

| Flag | Strategy | When to use |
|------|----------|-------------|
| (none) | Leave het sites unchanged | **Default. Works for virtually all cases.** |
| `--force` | Pick the allele with higher sequencing depth | You want maximum swaps; heuristic, does not infer phase |
| `--iupac` | Encode both alleles as IUPAC code (A/G→R) | Experimental. Produces non-standard VCF; most tools won't handle it |

## Output

- Auto-named: `toy_ref-1_H3.vcf.gz` (`<input>_ref-<sample>.vcf.gz`)
- VCF header records the sample used: `##swapRef_Sample=<ID=1_H3,...>`
- PL fields are dropped (no longer valid after allele reordering). Regenerate with `bcftools +fill-tags`.

## Verification

Allele re-encoding is a bijection — allele-sharing between samples is preserved, so:

- IBS distance matrix is **identical** (correlation = 1.0)
- Tree topology is **identical** (cophenetic correlation = 1.0)
- **Round-trip**: swap → inverse swap restores the original VCF exactly

`verify_swap.py` runs all three checks automatically.

## Important: This is NOT a Genome Liftover

This tool reorients **VCF allele definitions**, not genomic coordinates. It does not:

- Change CHROM, POS, or coordinate systems
- Generate new FASTA sequences
- Remap reads or call variants against a new assembly

For coordinate-based liftover between assemblies, use CrossMap, Picard LiftoverVcf, or bcftools + chain files.

## All Options

```
python3 swap_ref.py input.vcf.gz \
    --sample SAMPLE        # required
    -o output.vcf.gz       # output file (auto-named if omitted)
    --force                # optional: use AD at heterozygous sites (heuristic)
    --iupac                # optional: [EXPERIMENTAL] IUPAC codes at het sites
    --new-sample-name NAME # column name for original reference (default: ORIGINAL_REF)
    -q                     # quiet mode
```
