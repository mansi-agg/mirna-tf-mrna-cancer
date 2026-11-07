# Build the reference sets of cancer cell-line ChIP-seq experiments used to
# refine TF target genes, from public metadata only (no manual Excel steps).
#
# Sources (downloaded into cell_lines/raw/ if missing):
#   Cellosaurus      https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt
#   ChIP-Atlas       https://chip-atlas.dbcls.jp/data/metadata/experimentList.tab
#
# For each cancer type, every hg38 "TFs and others" ChIP-Atlas experiment whose
# cell type belongs to that tissue, or maps to a cell line with that disease, is
# given a decision and a reason:
#   1. The ChIP target must be a protein, not a tag or DNA mark ("Epitope tags",
#      "Biotin", "BrdU", ...)
#   2. Cell type must map to a Cellosaurus cell line (vague labels such as
#      "Breast cancer cells" are excluded unless resolved in manual_overrides.csv)
#   3. Cellosaurus disease must match the cancer type, and the ChIP-Atlas tissue
#      must agree (guards against name clashes, e.g. the AML line "ME-1")
#   4. Cellosaurus category must be "Cancer cell line" (drops transformed,
#      immortalised, iPSC, hybridoma ... lines)
#   5. The line must not be flagged as misidentified / contaminated
#   6. The line must not be recorded in Cellosaurus as selected for drug
#      resistance, given an added gene, or had a gene knocked out
#   7. The experiment must not use a genetically modified, tagged or
#      drug-resistant version of the line (keywords in title and sample metadata)
# Experiments resolved by hand go in manual_overrides.csv, which wins over these rules.
#
# Outputs, per cancer type, in cell_lines/:
#   <CANCER>_cell_lines.csv          included experiments (drop-in for the pipeline:
#                                    has `SRX ID` and `Track type` columns)
#   <CANCER>_experiments_review.csv  every experiment considered, with decision and reason
#   <CANCER>_tf_coverage.csv         per TF: number of experiments and distinct cell lines
#
# Run from the project root:  Rscript cell_lines/build_cell_lines.R

library(here)
library(tidyverse)

OUT_DIR <- here("cell_lines")
RAW_DIR <- here("cell_lines", "raw")
CELLOSAURUS_FILE <- file.path(RAW_DIR, "cellosaurus.txt")
CHIP_ATLAS_FILE <- file.path(RAW_DIR, "experimentList.tab")
MANUAL_OVERRIDES_FILE <- file.path(OUT_DIR, "manual_overrides.csv")

GENOME <- "hg38"

# ChIP-Atlas "Cell type class" and Cellosaurus disease terms for each cancer.
# TCGA-SKCM is cutaneous melanoma, so eye melanomas are excluded.
CANCERS <- list(
  SKCM = list(
    label = "melanoma",
    chip_atlas_classes = c("Epidermis"),
    disease_include = "melanoma",
    disease_exclude = "uveal|ocular|conjunctival|choroid"
  ),
  BRCA = list(
    label = "breast cancer",
    chip_atlas_classes = c("Breast"),
    disease_include = "breast.*(carcinoma|cancer)|(carcinoma|cancer).*breast",
    disease_exclude = "phyllodes|neuroendocrine"
  )
)

# Experiments mentioning these in their title or sample metadata are treated as
# modified versions of the cell line. Short tokens are matched case-sensitively
# as whole words to avoid false hits (e.g. "HA" inside other words).
MODIFICATION_PATTERN <- paste(
  "knock.?down", "knock.?out", "knock.?in", "shRNA", "siRNA", "sgRNA", "gRNA", "CRISPR",
  "Cas9", "dCas9", "over.?express", "transfect", "transduc", "lentivir", "retrovir",
  "inducible", "doxycycline", "tet.?on", "tet.?off", "degron", "dTAG", "auxin",
  "mutant", "empty vector", "scrambl", "non.?targeting", "stably express",
  "resistant", "resistance", "LTED", "tamoxifen.?resistant",
  sep = "|"
)
TAG_PATTERN <- "\\b(FLAG|HA|V5|GFP|eGFP|EGFP|YFP|mCherry|Myc-tag|BirA|AID|FKBP|KO|KD|OE|sh[A-Z0-9]+|si[A-Z0-9]{2,})\\b"

# ChIP-Atlas "TFs and others" tracks that are tags or DNA marks, not a TF
NON_TF_TRACKS <- c("Epitope tags", "Biotin", "GFP", "Cas9", "BrdU", "5-mC", "5-hmC",
                   "8-Hydroxydeoxyguanosine", "G-quadruplex")

# Drug or hormone treatments are flagged for review but not excluded
TREATMENT_PATTERN <- "treat|estradiol|E2\\b|DMSO|vehicle|inhibitor|dexamethasone|TNF|IFN|hypoxi|starv|stimulat"


# ---- Download ----

download_if_missing <- function(url, path) {
  if (file.exists(path)) return(invisible(path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  message("Downloading ", url)
  options(timeout = 3600)
  download.file(url, path, mode = "wb")
}

download_if_missing("https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt", CELLOSAURUS_FILE)
download_if_missing("https://chip-atlas.dbcls.jp/data/metadata/experimentList.tab", CHIP_ATLAS_FILE)


# ---- Read Cellosaurus ----

# Turn the Cellosaurus text file into a table with one row per human cell line.
#
# In the file, each cell line is a block of lines, and each line starts with a
# two-letter code saying what it holds (ID = name, AC = accession, SY = synonyms,
# DI = disease, OX = species, HI = parent line, CA = category, CC = comment).
#
# Steps:
#   1. Keep only the lines we need, including four kinds of comment (CC):
#      "Problematic cell line", plus three that mark a modified line:
#      "Selected for resistance to" (drug-resistant version), "Genetic integration"
#      (gene or construct added) and "Knockout cell" (gene knocked out).
#   2. Work out which cell line each kept line belongs to. Each cell line's
#      block starts with its "ID" line, so counting the ID lines seen so far
#      gives every line a cell line number (`entry`): 1 for the first cell
#      line's lines, 2 for the second's, and so on.
#   3. Label each line with its code (`field`), giving the two kinds of comment
#      their own labels, PR (problematic) and MO (modified), so they don't mix.
#   4. For each code, join each cell line's values into one text (a line can
#      have several synonyms or diseases), lined up by cell line number.
#   5. Build the table, tidy the disease names, and keep human lines only.
read_cellosaurus <- function(path) {
  lines <- readLines(path)
  lines <- lines[grepl("^(ID|AC|SY|DI|OX|HI|CA)   ", lines) |
                   grepl("^CC   (Problematic cell line|Selected for resistance to|Genetic integration|Knockout cell):", lines)]
  entry <- cumsum(startsWith(lines, "ID   "))
  field <- substr(lines, 1, 2)
  field[startsWith(lines, "CC   Problematic")] <- "PR" # problematic
  field[grepl("^CC   (Selected for resistance to|Genetic integration|Knockout cell):", lines)] <- "MO" # modified
  value <- substring(lines, 6)

  get <- function(f) tapply(value[field == f], entry[field == f], paste, collapse = " | ")
  ids <- unique(entry)
  pick <- function(f) { x <- get(f); as.vector(x[as.character(ids)]) }

  tibble(
    cellosaurus_name = pick("ID"),
    cellosaurus_id   = pick("AC"),
    synonyms         = pick("SY"),
    disease          = gsub("(NCIt|ORDO); [^;]+; ", "", pick("DI")),
    species          = pick("OX"),
    parent_line      = pick("HI"),
    category         = pick("CA"),
    problematic      = sub("^Problematic cell line: ", "", pick("PR")),
    line_modification = pick("MO")
  ) %>%
    filter(str_detect(species, "Homo sapiens"))
}

# Lower-case, alphanumeric only, so "SK-MEL-28", "SKMEL28" and "SK MEL 28" match
normalise_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))

# Every name and synonym of every line -> its Cellosaurus accession
build_name_lookup <- function(cellosaurus) {
  primary <- cellosaurus %>%
    transmute(key = normalise_name(cellosaurus_name), cellosaurus_id, is_primary = TRUE)
  synonyms <- cellosaurus %>%
    filter(!is.na(synonyms)) %>%
    transmute(key = synonyms, cellosaurus_id) %>%
    separate_rows(key, sep = "; ") %>%
    mutate(key = normalise_name(key), is_primary = FALSE)

  bind_rows(primary, synonyms) %>%
    filter(key != "") %>%
    distinct(key, cellosaurus_id, .keep_all = TRUE)
}

# Map each ChIP-Atlas cell type to one Cellosaurus line. A primary-name match
# wins over a synonym match; names that still match several lines are ambiguous.
match_cell_types <- function(cell_types, name_lookup) {
  tibble(cell_type = unique(cell_types), key = normalise_name(unique(cell_types))) %>%
    left_join(name_lookup, by = "key", relationship = "many-to-many") %>%
    group_by(cell_type) %>%
    filter(if (any(is_primary %in% TRUE)) is_primary %in% TRUE else TRUE) %>%
    summarise(n_matches = n_distinct(cellosaurus_id, na.rm = TRUE),
              cellosaurus_id = if (n_matches == 1) first(na.omit(cellosaurus_id)) else NA_character_,
              .groups = "drop")
}


# ---- Read ChIP-Atlas ----

# hg38 "TFs and others" experiments, with the title and sample metadata as text
read_chip_atlas_tf_experiments <- function(path, genome) {
  lines <- readLines(path)
  lines <- lines[grepl(paste0("\t", genome, "\tTFs and others\t"), lines, fixed = TRUE)]
  fields <- strsplit(lines, "\t", fixed = TRUE)

  # Column 9 is the title; columns 10 onwards are the sample metadata
  get_title <- function(x) {
    if (length(x) >= 9) x[9] else ""
  }
  get_metadata <- function(x) {
    if (length(x) > 9) paste(x[10:length(x)], collapse = " | ") else ""
  }

  tibble(
    `SRX ID`           = map_chr(fields, 1),
    Genome             = map_chr(fields, 2),
    `Track type class` = map_chr(fields, 3),
    `Track type`       = map_chr(fields, 4),
    `Cell type class`  = map_chr(fields, 5),
    `Cell type`        = map_chr(fields, 6),
    title              = map_chr(fields, get_title),
    sample_metadata    = map_chr(fields, get_metadata)
  )
}


# ---- Decide on each experiment ----

classify_experiments <- function(experiments, cellosaurus, cell_type_map, settings, overrides) {
  disease_ok <- function(d) {
    !is.na(d) & str_detect(tolower(d), settings$disease_include) &
      !str_detect(tolower(d), settings$disease_exclude)
  }
  relevant_lines <- cellosaurus$cellosaurus_id[disease_ok(cellosaurus$disease)]

  experiments %>%
    left_join(cell_type_map, by = c("Cell type" = "cell_type")) %>%
    filter(`Cell type class` %in% settings$chip_atlas_classes |
             cellosaurus_id %in% relevant_lines) %>%
    left_join(cellosaurus, by = "cellosaurus_id") %>%
    mutate(
      # Search the title and the metadata values, not the field names
      text = paste(title, gsub("(^|\\| )[^|=]+=", " ", sample_metadata)),
      modification_keywords = map_chr(str_extract_all(text, regex(MODIFICATION_PATTERN, ignore_case = TRUE)),
                                      ~ paste(unique(tolower(.x)), collapse = ", ")),
      tag_keywords = map_chr(str_extract_all(text, TAG_PATTERN), ~ paste(unique(.x), collapse = ", ")),
      treatment_keywords = map_chr(str_extract_all(text, regex(TREATMENT_PATTERN, ignore_case = TRUE)),
                                   ~ paste(unique(tolower(.x)), collapse = ", ")),
      reason = case_when(
        `Track type` %in% NON_TF_TRACKS ~ "ChIP target is a tag or DNA mark, not a TF",
        n_matches > 1 ~ "Cell type matches several Cellosaurus lines",
        is.na(cellosaurus_id) ~ "Cell type is not a named cell line",
        category != "Cancer cell line" ~ paste0("Not a cancer cell line (", category, ")"),
        is.na(disease) ~ "No disease recorded in Cellosaurus (check by hand)",
        !disease_ok(disease) ~ "Wrong disease",
        !`Cell type class` %in% settings$chip_atlas_classes ~
          paste0("Name matches a ", settings$label, " line but ChIP-Atlas tissue is ", `Cell type class`, " (check by hand)"),
        !is.na(problematic) ~ "Misidentified or contaminated line",
        !is.na(line_modification) ~ "Line is drug-resistant, has an added gene or a knockout (Cellosaurus)",
        modification_keywords != "" ~ "Modified cell line (keywords)",
        tag_keywords != "" ~ "Tagged or knockdown/knockout sample (keywords)",
        TRUE ~ "Included"
      ),
      decision = if_else(reason == "Included", "include", "exclude")
    ) %>%
    apply_manual_overrides(overrides) %>%
    select(decision, reason, `SRX ID`, Genome, `Track type class`, `Track type`,
           `Cell type class`, `Cell type`, cellosaurus_id, cellosaurus_name, disease,
           category, parent_line, problematic, line_modification, modification_keywords, tag_keywords,
           treatment_keywords, title, sample_metadata) %>%
    arrange(desc(decision), `Cell type`, `Track type`)
}

# manual_overrides.csv (optional): columns `SRX ID`, decision ("include"/"exclude"), note.
# Use it for experiments resolved by hand, e.g. vague cell types checked on GEO.
apply_manual_overrides <- function(df, overrides) {
  if (is.null(overrides)) return(df)
  df %>%
    left_join(select(overrides, `SRX ID`, override_decision = decision, override_note = note),
              by = "SRX ID") %>%
    mutate(decision = coalesce(override_decision, decision),
           reason = if_else(is.na(override_decision), reason,
                            paste0("Manual: ", override_note, " (was: ", reason, ")"))) %>%
    select(-override_decision, -override_note)
}

# Per TF: experiments and distinct cell lines among the included experiments
tf_coverage <- function(included) {
  included %>%
    group_by(TF = `Track type`) %>%
    summarise(n_experiments = n(),
              n_cell_lines = n_distinct(cellosaurus_id),
              cell_lines = paste(sort(unique(cellosaurus_name)), collapse = "; "),
              .groups = "drop") %>%
    arrange(desc(n_cell_lines), desc(n_experiments))
}


# ---- Run ----

message("Reading Cellosaurus")
cellosaurus <- read_cellosaurus(CELLOSAURUS_FILE)
name_lookup <- build_name_lookup(cellosaurus)

message("Reading ChIP-Atlas experiment list")
experiments <- read_chip_atlas_tf_experiments(CHIP_ATLAS_FILE, GENOME)
cell_type_map <- match_cell_types(experiments$`Cell type`, name_lookup)

overrides <- if (file.exists(MANUAL_OVERRIDES_FILE)) {
  read_csv(MANUAL_OVERRIDES_FILE, show_col_types = FALSE, comment = "#")
} else NULL

for (cancer in names(CANCERS)) {
  review <- classify_experiments(experiments, cellosaurus, cell_type_map, CANCERS[[cancer]], overrides)
  included <- filter(review, decision == "include")

  write_csv(review, file.path(OUT_DIR, paste0(cancer, "_experiments_review.csv")), na = "")
  write_csv(select(included, `SRX ID`, Genome, `Track type class`, `Track type`,
                   `Cell type class`, `Cell type`, cellosaurus_id, treatment_keywords),
            file.path(OUT_DIR, paste0(cancer, "_cell_lines.csv")), na = "")
  write_csv(tf_coverage(included), file.path(OUT_DIR, paste0(cancer, "_tf_coverage.csv")))

  message("\n", cancer, ": ", nrow(included), " of ", nrow(review), " experiments included (",
          n_distinct(included$cellosaurus_id), " cell lines, ",
          n_distinct(included$`Track type`), " TFs)")
  print(count(review, reason, sort = TRUE), n = Inf)
}
