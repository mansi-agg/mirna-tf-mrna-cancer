# Cancer cell-line ChIP-seq reference sets

Rebuilds the melanoma (SKCM) and breast cancer (BRCA) cell-line ChIP-seq lists used to refine TF target genes. The old lists were made by hand; these are made in code from public metadata, and every experiment gets a decision and a reason.

```bash
Rscript cell_lines/build_cell_lines.R
```

The first run downloads two metadata files into `cell_lines/raw/` (about 480 MB, not tracked in git):

- [Cellosaurus](https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt): each cell line's disease, category (cancer line, immortalised, iPSC, ...), synonyms, and flags for misidentified, drug-selected or genome-modified lines.
- [ChIP-Atlas experiment list](https://chip-atlas.dbcls.jp/data/metadata/experimentList.tab): every ChIP-seq experiment with its cell type, title and GEO/SRA sample metadata.

## Rules

Every hg38 "TFs and others" experiment in the matching ChIP-Atlas tissue (Epidermis for SKCM, Breast for BRCA), or on a cell line with the matching disease, is checked in this order. The first rule it fails is its reason:

| # | Rule | Reason shown |
|---|---|---|
| 1 | ChIP target is a protein, not a tag or DNA mark | ChIP target is a tag or DNA mark, not a TF |
| 2 | Cell type is a specific, unambiguous cell line in Cellosaurus | Cell type is not a named cell line / matches several Cellosaurus lines |
| 3 | Cellosaurus category is "Cancer cell line" | Not a cancer cell line (...) |
| 4 | Cellosaurus disease is melanoma (not uveal, ocular or conjunctival) or breast carcinoma | No disease recorded / Wrong disease |
| 5 | ChIP-Atlas tissue agrees with the disease, to catch name clashes such as the AML line "ME-1" | Name matches a ... line but ChIP-Atlas tissue is ... |
| 6 | Not misidentified or contaminated | Misidentified or contaminated line |
| 7 | Line not recorded in Cellosaurus as drug-resistant ("Selected for resistance to"), carrying an added gene ("Genetic integration") or a knockout ("Knockout cell") | Line is drug-resistant, has an added gene or a knockout (Cellosaurus) |
| 8 | Experiment title and sample metadata don't mention knockdown, knockout, shRNA, siRNA, CRISPR, overexpression, inducible or dox, resistant lines, and so on | Modified cell line (keywords) |
| 9 | No tag or knockdown/knockout labels (FLAG, HA, GFP, KO, KD, sh..., si...) | Tagged or knockdown/knockout sample (keywords) |

Drug and hormone treatments (E2, DMSO, inhibitors, ...) are **not** excluded. They are listed in `treatment_keywords` so you can decide.

## Outputs

| File | What it is |
|---|---|
| `<CANCER>_experiments_review.csv` | Every experiment considered, with `decision`, `reason`, the Cellosaurus fields, matched keywords, title and sample metadata. **Review this.** |
| `<CANCER>_cell_lines.csv` | Included experiments. Has the `SRX ID` and `Track type` columns, so it can replace the old file in `R/config.R` (`CELL_LINES_FILE`) |
| `<CANCER>_tf_coverage.csv` | Per TF: number of included experiments and distinct cell lines |
| `manual_overrides.csv` | Your decisions after checking by hand (`SRX ID`, `decision`, `note`). These override the rules; rerun the script after editing |

## Results (Cellosaurus 22 Oct 2026, ChIP-Atlas 14 Oct 2026)

| | Old list | New list |
|---|---|---|
| **SKCM** TF experiments | 557 (35 cell types) | 94 (16 cell lines) |
| **SKCM** TFs with data in `data_files/TF_targets/` | 36 | 13 (8 in only one cell line) |
| **BRCA** TF experiments | 4,495 (66 cell types) | 3,210 (40 cell lines) |
| **BRCA** TFs with data in `data_files/TF_targets/` | 134 | 123 (83 in only one cell line) |

- **SKCM:** most of the old list was not melanoma: A-673, IMR-32, HT-1080, SK-N-MC, SU-CCS-1, WI-38, LN-229 and others. Once they are removed, very few TFs have melanoma ChIP-seq, which supports Miles's doubt about the cell-line filter for melanoma.
- **BRCA:** the old list was mostly right. The main removals are modified or tagged MCF-7, T-47D and MDA-MB-231 experiments, vague labels such as "Breast cancer cells", and MCF10A-derived lines.

## To check by hand

- **Probably genuine, but tissue mislabelled in ChIP-Atlas:** COLO 829 (melanoma; ChIP-Atlas tissue "Digestive tract") and SUM 185PE (breast; "Others"). Add them to `manual_overrides.csv` if you agree.
- **HCC1187** (breast cancer line): Cellosaurus records no disease for it.
- **Vague cell types** ("Breast cancer cells", "Breast cancer", "Melanoma", "Melanoma cell line", "Mammary epithelial cells"): excluded unless resolved on GEO, as Peinan suggested. Each experiment's title and metadata in the review file usually say which line it is.
- **Keyword false positives:** skim the "Modified cell line" and "Tagged" rows. A control sample from a knockdown study (e.g. "siControl") is excluded on purpose, because it was still transfected.
- **Fusion and drug-resistant derivatives** with their own ChIP-Atlas names (e.g. FULVR, C70R, FULR) are excluded as "not a named cell line".
