# Single-cell experiments

This file lists the single-cell experiments of the benchmark and how each was processed. A comparison is one treated group against one control group. The final tables hold 15,311 comparisons from 100 datasets, and 70 of them are selected for evaluation.

The `single-cell-processing` branch holds the code that builds the per-gene tables, as `v2-processing` does for bulk. The `single-cell-evaluation` branch scores them, as `v2-evaluation` does for bulk.

| Group | Comparisons | Datasets | What changed in the cells | Use |
|---|---:|---:|---|---|
| Main | 5 | 4 | One miRNA raised or lowered | Scored against the 21 bulk experiments |
| Machinery | 62 | 11 | A gene that makes miRNAs or carries out repression was knocked down or out | 25 miRNA-making comparisons are plotted; all 62 are in the tables |
| Supplement | 3 | 3 | Two miRNAs raised together, or one miRNA knocked out in a clone | Tables only |

## Main experiments

| Comparison | Accession | miRNA | Cells and time | Treated | Control | Treated cells | Control cells | Replicates per side | Design |
|---|---|---|---|---|---|---|---|---|---|
| `miR-874-3p_mimic_MCF7_36h` | GSE214508 | miR-874-3p | MCF7 breast cancer cells, 36 h | precursor mimic | control miRNA | 8,057 | 13,114 | 3 pseudo-replicates | `~ condition` |
| `mimic-365_vs_antimiR-365_healthy_iPSC-CM_48h` | GSE185688 | miR-365a-3p | iPSC heart muscle cells, healthy line, 25 nM, 48 h | mimic | antimiR-365 | 2,935 | 2,389 | 3 pseudo-replicates | `~ condition` |
| `miR-98-5p_mimic_H9_d02_d03` | GSE280065 | miR-98-5p | H9 stem cells during resetting to naive pluripotency, days 2 and 3 | mimic | control mimic | 5,459 | 5,906 | 2 real replicates | `~ batch + condition` |
| `miR-192-5p_LNP_mimic_PBMC_24h` | GSE292981 | miR-192-5p | blood immune cells (PBMC), 24 h | mimic in lipid nanoparticles | saline | 7,260 | 5,872 | 3 pseudo-replicates | `~ condition` |
| `antimiR-365_vs_antimiR-Ctrl_LQT1_iPSC-CM_48h` | GSE185688 | miR-365a-3p | iPSC heart muscle cells, LQT1 patient line, 25 nM, 48 h | antimiR | control antimiR | 2,039 | 651 | 3 pseudo-replicates | `~ condition` |

Outcome of each comparison. States are at d = 0.5, on the 19,226 listed genes. The seed effect is how much more genes with a 7mer or 8mer site moved in the expected direction than matched genes with no site, in log2FC. It is reported and never used to keep or drop an experiment.

| Comparison | Responsive | No response | Unresolved | Seed effect and rank of 2,022 families |
|---|---|---|---|---|
| `miR-874-3p_mimic_MCF7_36h` | 2,384 | 8,932 | 7,910 | miR-874-3p: 0.300, rank 1 |
| `mimic-365_vs_antimiR-365_healthy_iPSC-CM_48h` | 87 | 15,623 | 3,516 | miR-365a-3p: 0.103, rank 1 |
| `miR-98-5p_mimic_H9_d02_d03` | 9 | 17,990 | 1,227 | miR-98-5p: 0.054, rank 1 |
| `miR-192-5p_LNP_mimic_PBMC_24h` | 570 | 11,490 | 7,166 | miR-192-5p: 0.014, rank 502 |
| `antimiR-365_vs_antimiR-Ctrl_LQT1_iPSC-CM_48h` | 113 | 14,515 | 4,598 | miR-365a-3p: 0.014, rank 315 |

Sources and cell filters:

| Comparison | Paper | Cell filter |
|---|---|---|
| `miR-874-3p_mimic_MCF7_36h` | 10.1038/s41598-022-23205-w | at most 10% mitochondrial, in the deposit's Cell Ranger filter (authors) |
| `mimic-365_vs_antimiR-365_healthy_iPSC-CM_48h` | 10.1038/s41467-021-27856-7 | 300 to 6,500 genes, at most 20% mitochondrial (authors) |
| `miR-98-5p_mimic_H9_d02_d03` | preprint 10.1101/2025.10.14.682465 | 2,000 or more genes, 5,000 to 100,000 counts, at most 15% mitochondrial (authors) |
| `miR-192-5p_LNP_mimic_PBMC_24h` | none recorded | barcodes with a hashtag call in the authors' table; no thresholds stated |
| `antimiR-365_vs_antimiR-Ctrl_LQT1_iPSC-CM_48h` | 10.1038/s41467-021-27856-7 | 300 to 6,500 genes, at most 20% mitochondrial (authors) |

Caveats:

- `miR-874-3p_mimic_MCF7_36h`: one library per arm: three technical pseudo-replicates, labels mean changed in this culture; miR-874 drives apoptosis by 36 h, dying cells may be included.
- `mimic-365_vs_antimiR-365_healthy_iPSC-CM_48h`: one library per arm: three technical pseudo-replicates, labels mean changed in this culture; no neutral control: raising against lowering the same miRNA, so effects are not one-sided.
- `miR-98-5p_mimic_H9_d02_d03`: days 2 and 3 stand as real replicates, day is the batch; let-7 family.
- `miR-192-5p_LNP_mimic_PBMC_24h`: one library per arm: three technical pseudo-replicates, labels mean changed in this culture; hashed in one lane; control is saline, not a control-mimic nanoparticle, so the carrier response is in the effect; GEO says saline 24 h, the paper says 4 h; fractional STARsolo counts rounded per cell for the cell store.
- `antimiR-365_vs_antimiR-Ctrl_LQT1_iPSC-CM_48h`: one library per arm: three technical pseudo-replicates, labels mean changed in this culture; hashed in one lane; arms unbalanced 3:1.

### What the replicates mean

Only miR-98 has real replicates: day 2 and day 3 each hold one mimic and one control library, and the day is the batch. The other four have one sample per group. Their cells are split at random into three pseudo-replicates, which measure cell sampling and not differences between cultures, so their labels mean changed in this culture.

The treated and control cells of miR-874 and of the miR-365 mimic comparison are separate libraries. The control check below shows that this alone gives hundreds to thousands of genes called responsive, so their AUROC values are unreliable. The cells of miR-192 and of the miR-365 antimiR comparison are hashed into one lane and share a library. Each group is still one culture, and that difference was not tested.

## Machinery experiments

When a gene needed to make miRNAs or to carry out repression is lowered, targets of the miRNAs the cell expresses should rise. No single miRNA is changed, so these experiments are scored on the miRNAs each cell line expresses.

In the comparison IDs, `_KO` means the gene was targeted. All screens are CRISPRi knockdowns except the A549 screen, which is a Cas9 knockout.

| Screen | Source | Cell line | Method | Replicate units | Median guides per cell | miRNA-making comparisons | Repression comparisons | Counts |
|---|---|---|---|---|---|---|---|---|
| `Replogle2022_K562_GW` | figshare+ 20029387 | K562 | CRISPRi | 33 sets of 8 lanes (267 lanes) | 1 | 3 | 5 | authors' raw counts; only the cells and the 8,248 genes the authors kept |
| `Replogle2022_K562_essential` | figshare+ 20029387 | K562 | CRISPRi | 6 sets of 8 lanes (48 lanes) | 1 | 2 | 0 | authors' raw counts; only the cells and the 8,563 genes the authors kept |
| `VIPerturb_K562` | Zenodo 18460279 | K562 | CRISPRi with KRAB-dCas9-MeCP2 | 10x lanes of one culture | 1 | 2 | 4 | raw counts of the three full bins; shallow depth per gene; the response-filtered object is not used |
| `XAtlasOrion_HEK293T` | figshare+ 29190726, through HKU 10.25442/hku.34002315.v1 | HEK293T | CRISPRi | 27 sets of 8 lanes (223 lanes) | 1 | 3 | 6 | counts reconstructed from the HKU log-normalised file, which leaves out perturbations below 20 cells |
| `XAtlasOrion_HCT116` | figshare+ 29190726, through HKU 10.25442/hku.34002315.v1 | HCT116 | CRISPRi | 13 sets of 8 lanes (109 lanes) | 1 | 3 | 6 | counts reconstructed from the HKU log-normalised file, which leaves out perturbations below 20 cells |
| `GSE218566_HEK293` | GSE218566 | HEK293 | CRISPRi | library dates (sub-libraries of one pooled culture) | 1 | 3 | 1 | authors' whole-transcriptome counts; culture replicates are not labelled per cell |
| `Nadig2025_HepG2` | GSE264667 | HepG2 | CRISPRi | 7 sets of 8 lanes (56 lanes) | 1 | 1 | 1 | authors' raw counts; only the cells and the 9,624 genes the authors kept |
| `Nadig2025_Jurkat` | GSE264667 | Jurkat | CRISPRi | 6 sets of 8 lanes (55 lanes) | 1 | 2 | 1 | authors' raw counts; only the cells and the 8,882 genes the authors kept |
| `GSE345058_A549` | GSE345058 | A549 | inducible Cas9 knockout | the authors' 4 batches of 8 lanes | 1 | 1 | 1 | authors' raw counts, 1,000 targeted genes |
| `GSE249595_Jurkat_activated` | GSE249595 | Jurkat, activated | CRISPRi at high guide load | 2 sets of 8 lanes (16 lanes) | 11 | 4 | 7 | raw barcodes; the authors' cell, hash and guide filters were reapplied; their mitochondrial outlier filter gives no threshold and is not reapplied |
| `GSE249595_Jurkat_untreated` | GSE249595 | Jurkat, untreated | CRISPRi at high guide load | 2 sets of 8 lanes (16 lanes) | 12 | 1 | 5 | raw barcodes; the authors' cell, hash and guide filters were reapplied; their mitochondrial outlier filter gives no threshold and is not reapplied |

The 62 comparisons are the ones that met two rules, out of 201 comparisons of the 11 genes in the final tables: a measured miRNA profile for the cell line (110 have none) and at least 50 treated cells (29 have fewer). The profiles are small RNA sequencing for K562, HEK293T, HepG2 and A549, and the CCLE panel for HCT116 and Jurkat. CCLE names without an arm go to the arm with more reads in the sequencing profiles, because miRBase releases old enough to hold those names could not be downloaded.

Experiments of one screen share their control cells, so the machinery plot shows one point per screen, the median of its experiments. Repression genes are left off the plot. That choice was made after seeing their weak signal, and these genes have paralogues that may compensate. They stay in the tables.

The two `GSE249595` screens differ from the rest in two ways. Cells carry a median of 11 or 12 guides, so the control is the cells of the same lane set that carry no guide of the target or of the other listed pathway genes. Each side has only 2 replicate units. Their 5 miRNA-making comparisons score below 0 on decoy z for all predictors, about -2 at the screen median, and none has a responsive gene. The shrunken log2FC of the targeted gene is close to 0 even where the unshrunken log2FC is clear (`DICER1`, -1.65 and -2.66), so shrinkage may be flattening these comparisons. This has not been tested.

| Comparison | Gene type | Treated cells | Control cells | Replicates per side | Responsive | Target log2FC | Target shrunken log2FC |
|---|---|---|---|---|---|---|---|
| `GSE218566_HEK293__DICER1_KO` | miRNA-making | 375 | 2,758 | 11 | 12 | -2.72 | -2.67 |
| `GSE218566_HEK293__DROSHA_KO` | miRNA-making | 574 | 2,758 | 11 | 18 | -1.67 | -1.60 |
| `GSE218566_HEK293__XPO5_KO` | miRNA-making | 1,583 | 2,758 | 11 | 0 | -1.11 | -1.07 |
| `GSE249595_Jurkat_activated__DGCR8_KO` | miRNA-making | 138 | 622,130 | 2 | 0 | -0.31 | -0.02 |
| `GSE249595_Jurkat_activated__DICER1_KO` | miRNA-making | 716 | 622,130 | 2 | 0 | -1.65 | -0.00 |
| `GSE249595_Jurkat_activated__DROSHA_KO` | miRNA-making | 308 | 622,130 | 2 | 0 | -0.13 | -0.00 |
| `GSE249595_Jurkat_activated__XPO5_KO` | miRNA-making | 425 | 622,130 | 2 | 0 | 0.18 | 0.00 |
| `GSE249595_Jurkat_untreated__DICER1_KO` | miRNA-making | 103 | 76,195 | 2 | 0 | -2.66 | -0.00 |
| `GSE345058_A549__DGCR8_KO` | miRNA-making | 569 | 45,320 | 4 | 48 | -0.24 | -0.04 |
| `Nadig2025_HepG2__XPO5_KO` | miRNA-making | 56 | 4,311 | 6 | 0 | -2.97 | -2.53 |
| `Nadig2025_Jurkat__DICER1_KO` | miRNA-making | 94 | 12,013 | 6 | 21 | -3.61 | -3.36 |
| `Nadig2025_Jurkat__XPO5_KO` | miRNA-making | 134 | 12,013 | 6 | 8 | -3.59 | -0.02 |
| `Replogle2022_K562_GW__DGCR8_KO` | miRNA-making | 122 | 61,590 | 27 | 87 | -1.64 | -1.52 |
| `Replogle2022_K562_GW__DICER1_KO` | miRNA-making | 342 | 75,328 | 33 | 275 | -3.06 | -3.04 |
| `Replogle2022_K562_GW__XPO5_KO` | miRNA-making | 476 | 75,328 | 33 | 11 | -2.38 | -2.32 |
| `Replogle2022_K562_essential__DICER1_KO` | miRNA-making | 158 | 10,691 | 6 | 11 | -4.95 | -0.01 |
| `Replogle2022_K562_essential__XPO5_KO` | miRNA-making | 221 | 10,691 | 6 | 2 | -2.61 | -2.53 |
| `VIPerturb_K562__DROSHA_KO` | miRNA-making | 68 | 7,949 | 2 | 29 | -1.44 | -1.26 |
| `VIPerturb_K562__XPO5_KO` | miRNA-making | 69 | 7,949 | 2 | 0 | -1.24 | -1.19 |
| `XAtlasOrion_HCT116__DGCR8_KO` | miRNA-making | 126 | 165,777 | 13 | 4 | -1.51 | -1.31 |
| `XAtlasOrion_HCT116__DICER1_KO` | miRNA-making | 82 | 152,668 | 12 | 0 | -1.50 | -1.31 |
| `XAtlasOrion_HCT116__DROSHA_KO` | miRNA-making | 214 | 165,777 | 13 | 32 | -1.80 | -1.71 |
| `XAtlasOrion_HEK293T__DGCR8_KO` | miRNA-making | 71 | 111,420 | 13 | 3 | -0.99 | -0.75 |
| `XAtlasOrion_HEK293T__DICER1_KO` | miRNA-making | 87 | 157,826 | 18 | 10 | -0.94 | -0.91 |
| `XAtlasOrion_HEK293T__DROSHA_KO` | miRNA-making | 189 | 213,505 | 26 | 47 | -0.88 | -0.86 |
| `GSE218566_HEK293__AGO2_KO` | repression | 409 | 2,758 | 11 | 1 | -2.15 | -2.08 |
| `GSE249595_Jurkat_activated__AGO1_KO` | repression | 353 | 622,130 | 2 | 0 | -0.87 | -0.00 |
| `GSE249595_Jurkat_activated__AGO2_KO` | repression | 695 | 622,130 | 2 | 1 | -0.02 | -0.00 |
| `GSE249595_Jurkat_activated__AGO3_KO` | repression | 330 | 622,130 | 2 | 0 | -0.28 | -0.00 |
| `GSE249595_Jurkat_activated__AGO4_KO` | repression | 803 | 622,130 | 2 | 0 | -0.75 | -0.00 |
| `GSE249595_Jurkat_activated__TNRC6A_KO` | repression | 863 | 622,130 | 2 | 0 | -0.35 | -0.00 |
| `GSE249595_Jurkat_activated__TNRC6B_KO` | repression | 905 | 622,130 | 2 | 0 | -0.39 | -0.00 |
| `GSE249595_Jurkat_activated__TNRC6C_KO` | repression | 974 | 622,130 | 2 | 0 | 0.14 | 0.00 |
| `GSE249595_Jurkat_untreated__AGO2_KO` | repression | 100 | 76,195 | 2 | 0 | 0.59 | 0.00 |
| `GSE249595_Jurkat_untreated__AGO4_KO` | repression | 94 | 76,195 | 2 | 0 | 0.80 | 0.00 |
| `GSE249595_Jurkat_untreated__TNRC6A_KO` | repression | 126 | 76,195 | 2 | 0 | -0.44 | -0.00 |
| `GSE249595_Jurkat_untreated__TNRC6B_KO` | repression | 109 | 76,195 | 2 | 0 | -0.26 | -0.00 |
| `GSE249595_Jurkat_untreated__TNRC6C_KO` | repression | 120 | 76,195 | 2 | 0 | 0.04 | 0.00 |
| `GSE345058_A549__TNRC6A_KO` | repression | 501 | 45,320 | 4 | 20 | -0.31 | -0.28 |
| `Nadig2025_HepG2__TNRC6A_KO` | repression | 89 | 4,976 | 7 | 9 | -3.09 | -2.90 |
| `Nadig2025_Jurkat__TNRC6A_KO` | repression | 168 | 12,013 | 6 | 23 | -2.47 | -2.35 |
| `Replogle2022_K562_GW__AGO1_KO` | repression | 160 | 69,253 | 30 | 6 | -0.46 | -0.29 |
| `Replogle2022_K562_GW__AGO3_KO` | repression | 183 | 68,392 | 30 | 0 | -2.35 | -2.25 |
| `Replogle2022_K562_GW__TNRC6A_KO` | repression | 379 | 75,328 | 33 | 3 | -2.13 | -2.06 |
| `Replogle2022_K562_GW__TNRC6B_KO` | repression | 59 | 30,372 | 13 | 0 | -0.44 | -0.00 |
| `Replogle2022_K562_GW__TNRC6C_KO` | repression | 237 | 70,998 | 31 | 0 | -0.06 | -0.00 |
| `VIPerturb_K562__AGO3_KO` | repression | 102 | 7,949 | 2 | 0 | -1.42 | -0.00 |
| `VIPerturb_K562__AGO4_KO` | repression | 52 | 7,949 | 2 | 2 | -1.15 | -0.00 |
| `VIPerturb_K562__TNRC6A_KO` | repression | 65 | 7,949 | 2 | 1 | -1.28 | -0.00 |
| `VIPerturb_K562__TNRC6C_KO` | repression | 65 | 7,949 | 2 | 2 | -0.07 | -0.00 |
| `XAtlasOrion_HCT116__AGO1_KO` | repression | 159 | 165,777 | 13 | 2 | -1.66 | -0.02 |
| `XAtlasOrion_HCT116__AGO2_KO` | repression | 118 | 165,777 | 13 | 1 | -1.63 | -1.55 |
| `XAtlasOrion_HCT116__AGO3_KO` | repression | 123 | 165,777 | 13 | 0 | -1.33 | -1.21 |
| `XAtlasOrion_HCT116__AGO4_KO` | repression | 353 | 165,777 | 13 | 2 | -0.65 | -0.00 |
| `XAtlasOrion_HCT116__TNRC6B_KO` | repression | 242 | 165,777 | 13 | 0 | -0.06 | -0.00 |
| `XAtlasOrion_HCT116__TNRC6C_KO` | repression | 164 | 165,777 | 13 | 2 | -0.16 | -0.00 |
| `XAtlasOrion_HEK293T__AGO1_KO` | repression | 217 | 205,455 | 25 | 21 | -0.94 | -0.91 |
| `XAtlasOrion_HEK293T__AGO2_KO` | repression | 223 | 210,364 | 26 | 10 | -0.67 | -0.66 |
| `XAtlasOrion_HEK293T__AGO3_KO` | repression | 173 | 199,899 | 24 | 15 | -0.57 | -0.56 |
| `XAtlasOrion_HEK293T__AGO4_KO` | repression | 501 | 218,838 | 27 | 2 | -0.43 | -0.33 |
| `XAtlasOrion_HEK293T__TNRC6B_KO` | repression | 273 | 218,838 | 27 | 7 | -0.25 | -0.25 |
| `XAtlasOrion_HEK293T__TNRC6C_KO` | repression | 199 | 206,182 | 25 | 11 | 0.16 | 0.16 |

Target log2FC is information only, and no comparison is kept or dropped on it. A target counts as readable with at least 1,000 control counts, and 22 of the 25 miRNA-making targets fell as expected.

## Supplement experiments

| Comparison | Accession | miRNA | Cells and time | Treated | Control | Treated cells | Control cells | Replicates per side | Design |
|---|---|---|---|---|---|---|---|---|---|
| `miR-18a_KO_A375` | GSE318353 | miR-18a-5p lost | A375 melanoma cells | knockout clone | wild type | 9,760 | 10,488 | 3 pseudo-replicates | `~ condition` |
| `miR-449a_miR-340-5p_RD18_60h` | GSE282611 | miR-449a and miR-340-5p raised together | RD18 rhabdomyosarcoma cells, single nuclei, 60 h | two miRNAs delivered together | lipofectamine only | 7,911 | 6,417 | 3 pseudo-replicates | `~ condition` |
| `let-7i_miR-452_OE_iPSC-CM` | GSE135595 | let-7i-5p and miR-452-5p raised together | iPSC heart muscle cells, glucose media | overexpression | empty vector | 354 | 144 | 3 pseudo-replicates | `~ condition` |

| Comparison | Responsive | No response | Unresolved | Seed effect and rank of 2,022 families |
|---|---|---|---|---|
| `miR-18a_KO_A375` | 1,935 | 9,851 | 7,440 | miR-18a-5p: 0.019, rank 882 |
| `miR-449a_miR-340-5p_RD18_60h` | 1,698 | 5,370 | 12,158 | miR-449a: -0.043, rank 1,226; miR-340-5p: 0.083, rank 117 |
| `let-7i_miR-452_OE_iPSC-CM` | 9 | 13,989 | 5,228 | let-7i-5p: 0.012, rank 877; miR-452-5p: -0.006, rank 1,490 |

Sources and cell filters:

| Comparison | Paper | Cell filter |
|---|---|---|
| `miR-18a_KO_A375` | 10.1002/advs.202515158 | 200 or more genes, 1,000 or more counts, at most 20% mitochondrial (authors) |
| `miR-449a_miR-340-5p_RD18_60h` | none recorded | deposit's Cell Ranger filter |
| `let-7i_miR-452_OE_iPSC-CM` | 10.1038/s41467-019-12482-1 | cells in the authors' metadata: at most 40% mitochondrial, 200 or more genes, 2,000 or more counts |

Caveats:

- `miR-18a_KO_A375`: one knockout clone, so clone effects cannot be separated from the loss of miR-18a; per-library matrices are not deposited (2 KO and 8 WT libraries merged), so three technical pseudo-replicates.
- `miR-449a_miR-340-5p_RD18_60h`: one library per arm: three technical pseudo-replicates, labels mean changed in this culture; two miRNAs delivered together; control is lipofectamine only, no control mimic; single nuclei.
- `let-7i_miR-452_OE_iPSC-CM`: one library per arm: three technical pseudo-replicates, labels mean changed in this culture; two miRNAs raised together; few cells.

## Comparisons not selected

| Comparisons | Count | Why not selected |
|---|---:|---|
| miR-199a GFP high against low | 3 | Fitted, left out of the install: about 280 GB of reads and alignment for sorted cells |
| miR-98 at days 8 and 14 | 2 | The seed effect is no longer visible, so this choice used a result |
| Reporter sorts (let-7 switch, miR-126 sensor) | 7 | Cells were sorted on a reporter and no miRNA was changed |
| `MIR9-2HG` enhancer deletions | 4 | The miRNA drop was never measured |
| Other pathway comparisons | 297 | 158 are of 8 other pathway genes (`DIS3L2`, `LIN28A`, `LIN28B`, `PRKRA`, `TARBP2`, `TUT4`, `TUT7`, `ZSWIM8`); 139 fail the two rules above |
| Non-coding host genes | 360 | 28 have a profile and at least 50 cells, and none of those hosts an expressed miRNA |
| Protein-coding host genes | 13,177 | Held back; about 66 have an expressed miRNA, and the loss of the host protein is mixed with any miRNA effect |
| Other fitted comparisons | 655 | 613 are from `GSE181338`, which is held out (308 knockouts, 305 host genes); the rest are other designs |
| Not fitted | 736 | Too few cells or replicate units; the reason is in the final table |

Not selected means left out of the evaluation. Every fitted comparison keeps its table, and these lists can change.

## How the data were processed

The same steps run for every comparison, and only the source files and the replicate units differ between datasets.

1. **Download.** `download.sh` fetches the 24 source files in `sources.tsv` and checks their sizes and MD5s. Main is about 2 GB, supplement about 3 GB and machinery about 228 GB. Reads are not realigned: the counts are the authors' or, for two screens, the HKU reprocessing.
2. **Cell store.** One adapter per dataset writes all cells to one file with raw counts, each cell's sample, batch, replicate unit, group and QC result. Features map to unversioned Ensembl 115 gene IDs through the GTF, and features with no unique ID are left out.
3. **Cell filter.** Each dataset uses the filter its authors stated, listed in the tables above. Where the authors deposit filtered cells, every deposited cell passes. Doublets are not removed.
4. **Replicate units.** In order of preference: real replicates such as donors, days or independent cultures; for screens, sets of about 8 consecutive lanes or the authors' batches; with one sample per group, three pseudo-replicates at random.
5. **Pseudobulk.** The passing cells are summed per batch, replicate unit and group. Units below the minimum (20 cells in main and supplement, 3 cells per lane set in machinery screens) are dropped, then units in batches that do not hold both groups, and each side needs at least 2 units. Fractional counts are rounded after summing, because DESeq2 takes integers.
6. **DESeq2 and shrinkage.** The sums go through the bulk script `run_deseq2_counts.R` with `~ batch + condition` when batches differ and otherwise `~ condition` (DESeq2 1.50.2, apeglm 1.32.0, R 4.5.3).
7. **Labels.** The bulk script `label_states.R` gives the shrunken log2FC, the s-values and the states on the 19,226 listed genes, oriented so that positive means the gene moved like a target. A gene is responsive when the oriented shrunken change is positive and the s-value against d is below 0.05, and no-response when the whole shrunken 95% interval of the oriented change is below d.

Genes a dataset does not hold, such as those the authors filtered out, have no counts and no statistics.

The commands are in `README.md`:

```bash
pipelines/single_cell/setup.sh
pipelines/single_cell/download.sh --groups main
pipelines/single_cell/run.sh --groups main
```

### Reproduction

All 70 selected gene tables were rebuilt with these scripts on 1 November 2026, and all are byte-identical to the final tables. The fresh test downloaded and built 18 comparisons from an empty folder in 33 minutes: main, supplement and the 10 machinery comparisons of two datasets. All 62 machinery comparisons were then rebuilt from copies of the source files on disk, so the downloads behind the other 52 were not repeated.

## Data checks

These checks were written down before they were run. None of them changes a table.

Control check. Two groups of cells with no perturbation between them should give about 0 responsive genes. Numbers are genes called responsive in the direction a knockout would give, then in the direction an overexpression would give.

| Two sides, nothing perturbed between them | Cells per side | Knockout direction | Overexpression direction |
|---|---:|---:|---:|
| One library split at random, miR-98 control | 1,857 | 0 | 0 |
| One library split at random, miR-874 control | 6,557 | 0 | 0 |
| One library split at random, antimiR-365 control | 326 | 0 | 1 |
| One hashed group split at random, saline | 2,936 | 3 | 1 |
| Two lane sets of one screen, K562 | about 2,100 | 0 | 0 |
| Two clones, organoid week 12 | 9,619 and 7,439 | 391 | 226 |
| Two libraries of one day, miR-98 controls day 14 | 3,713 and 3,486 | 991 | 1,720 |
| Two libraries of two days, miR-98 controls days 2 and 3 | 1,982 and 3,924 | 4,632 | 4,743 |

In screens, 300 fake comparisons of control cells in 20 screens gave at most 22 responsive genes, and the median was 0 in 16 of the 20 screens. In three screens real comparisons could not be told from these false calls at the median: X-Atlas HEK293T (12 false against 4.5 real), X-Atlas HCT116 (1 and 1) and Nourreddine KOLF, which is not selected (0.5 and 0).

miR-192 by cell type. The share of T cells is 48% with the mimic and 34% with saline. No cell type shows a seed effect, and most responsive genes are in monocytes (1,198), where a response to the nanoparticle cannot be separated because the control is saline.

## Evaluation

The evaluation scores each predictor on these tables with the bulk evaluator's measures.

- AUROC: chance that a responsive gene scores above a no-response gene, on comparisons with at least 30 genes in each state.
- Spearman rho: agreement of the score's ranking with the oriented shrunken log2FC.
- Decoy z-score: Spearman for the tested miRNA against different-seed miRNAs, in standard deviations. For machinery comparisons, the expressed miRNA families (at least 0.5% of reads) weighted by share against 20 sets of families below 0.01%.

In the single-cell plot, each box shows the 21 bulk experiments and each dot one main experiment. In the machinery plot, each dot is one screen.
