library(tidyverse, quietly = TRUE)

CACHE_DIR <- "/home/ood-coffm049/hbcdData/imputed/gp_cache"
ANCESTRY_FILE <- "/shared/release/hbcd/hbcd/rawdata/phenotype/sed_basic_demographics.tsv"

# ── Read typed sample list ──────────────────────────────────────────────
si <- read_table(file.path(CACHE_DIR, "sample_order_typed.txt"),
                 col_names = c("sid", "type"), col_types = "cc")
cat("sample_order_typed.txt:", nrow(si), "rows\n")
print(head(si, 5))

# ── Read demographics ───────────────────────────────────────────────────
dm <- read_delim(ANCESTRY_FILE, show_col_types = FALSE, delim = "\t")
cat("\nDemographics columns:", paste(names(dm), collapse = ", "), "\n")
cat("session_id values:", paste(unique(dm$session_id), collapse = ", "), "\n")
cat("Demographics rows:", nrow(dm), "\n")

# ── Check for the columns we need ───────────────────────────────────────
for (col in c("participant_id", "session_id",
              "sed_basic_demographics_child_race",
              "sed_basic_demographics_screen_mother_race")) {
    cat("Column", col, "present:", col %in% names(dm), "\n")
}

# ── Child lookup ────────────────────────────────────────────────────────
cat("\nses-V02 rows:", sum(dm$session_id == "ses-V02", na.rm = TRUE), "\n")
cat("ses-V02 + non-NA child race:",
    sum(dm$session_id == "ses-V02" & !is.na(dm$sed_basic_demographics_child_race), na.rm = TRUE), "\n")

child_lookup <- dm %>%
    filter(!is.na(.data[["sed_basic_demographics_child_race"]]),
           session_id == "ses-V02") %>%
    count(participant_id, sed_basic_demographics_child_race) %>%
    group_by(participant_id) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(sid = str_remove(tolower(participant_id), "^sub-"),
              ancestry = as.character(sed_basic_demographics_child_race))

cat("\nchild_lookup:", nrow(child_lookup), "rows\n")
print(head(child_lookup, 5))

# ── Mother lookup ───────────────────────────────────────────────────────
cat("\nses-V01 rows:", sum(dm$session_id == "ses-V01", na.rm = TRUE), "\n")
cat("ses-V01 + non-NA mother race:",
    sum(dm$session_id == "ses-V01" & !is.na(dm$sed_basic_demographics_screen_mother_race), na.rm = TRUE), "\n")

mother_lookup <- dm %>%
    filter(!is.na(.data[["sed_basic_demographics_screen_mother_race"]]),
           session_id == "ses-V01") %>%
    count(participant_id, sed_basic_demographics_screen_mother_race) %>%
    group_by(participant_id) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(sid = str_remove(tolower(participant_id), "^sub-"),
              ancestry = as.character(sed_basic_demographics_screen_mother_race))

cat("\nmother_lookup:", nrow(mother_lookup), "rows\n")
print(head(mother_lookup, 5))

# ── Test matches ────────────────────────────────────────────────────────
c_sids <- si$sid[si$type == "c"]
matched_c <- match(c_sids, child_lookup$sid)
cat("\nC-type samples:", length(c_sids), "\n")
cat("  Matched in child_lookup:", sum(!is.na(matched_c)), "\n")
if (length(c_sids) > 0) {
    cat("  First 10 sids:", paste(head(c_sids, 10), collapse = ", "), "\n")
    cat("  First 10 child_lookup sids:", paste(head(child_lookup$sid, 10), collapse = ", "), "\n")
}

m_sids <- si$sid[si$type == "m"]
matched_m <- match(m_sids, mother_lookup$sid)
cat("\nM-type samples:", length(m_sids), "\n")
cat("  Matched in mother_lookup:", sum(!is.na(matched_m)), "\n")
if (length(m_sids) > 0) {
    cat("  First 10 sids:", paste(head(m_sids, 10), collapse = ", "), "\n")
    cat("  First 10 mother_lookup sids:", paste(head(mother_lookup$sid, 10), collapse = ", "), "\n")
}

# ── What about raw participant_id values? ───────────────────────────────
cat("\n--- Raw participant_id from demographics (first 10) ---\n")
print(head(unique(dm$participant_id), 10))
cat("--- Raw sid from sample_order_typed (first 10) ---\n")
print(head(si$sid, 10))
