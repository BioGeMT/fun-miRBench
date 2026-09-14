# Common Resources

This directory contains external biological resources shared by multiple standardized predictor
pipelines and the benchmark.

The downloaded resource files are not tracked by Git. A pipeline or benchmark command downloads a
resource automatically when it is first required, then subsequent runs and other pipelines reuse
the cached copy.

Expected directory structure:

```text
data/common_resources/
├── README.md
├── mirbase/
│   └── mature.fa
└── ensembl/
    ├── Homo_sapiens.GRCh38.115.gtf.gz
    └── protein_coding_gene_ids.txt
```

- `mirbase/mature.fa` is the miRBase 22.1 mature-miRNA FASTA used by all standardized predictor
  pipelines.
- `ensembl/Homo_sapiens.GRCh38.115.gtf.gz` is the Ensembl release 115 annotation shared by
  microT-CNN, miRAW, TargetScan, and protein-coding benchmark filtering.
- `ensembl/protein_coding_gene_ids.txt` is derived from the shared GTF and cached by the benchmark.

Deleting any downloaded file is safe; it will be downloaded or regenerated the next time it is
needed.
