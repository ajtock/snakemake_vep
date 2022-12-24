# Snakemake workflow for annotating variants with pathogenicity predictions using ensembl-vep

* * *

NOTE: docker installation of ensemblorg/ensembl-vep is required, including genome FASTA, cache, and (optionally) plugins.

For docker installation instructions, see https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html

For use with custom files (e.g., ClinVar VCF, including TBI index), see examples at https://www.ensembl.org/info/docs/tools/vep/script/vep_custom.html

* * *

### Usage:

Run the following commands in the workflow base directory:

```
open -a Docker
conda activate snakemake
snakemake -p --cores
conda deactivate
```

### Useful Snakemake parameters

- `--cores` specifies the maximum number of threads
- `-n` performs a dry run
- `-p` prints commands
- `--use-conda`
- `--conda-prefix ~/.myconda`
- `--forcerun vcf_vep_anno` forces rerun of a given rule (e.g., `vcf_vep_anno`)
