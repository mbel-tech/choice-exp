# =============================================================================
# cld_letter_policy.R
# =============================================================================
# CANONICAL compact-letter-display policy for treatment x interval cells.
#
# Sourced by BOTH figures_step5_module.R (which draws the letters) and
# export_tukey_interval_module.R (which tabulates them). One definition, so a
# letter in the exported document cannot disagree with the same letter in the
# figure. Do not re-implement this rule anywhere else.
#
# -----------------------------------------------------------------------------
# THE RULE
# -----------------------------------------------------------------------------
# Cells are walked in a fixed reading order -- control first, then interval
# ascending -- and letters are renamed a, b, c, ... in order of FIRST
# APPEARANCE. Because control/interval 1 is first in that order, its letters are
# the first ones seen, so its group ALWAYS begins with "a".
#
# emmeans/multcomp assign letters in an order that depends on the fitted means,
# so without this step the same experimental cell carries a different letter in
# different panels and "a" lands on whichever cell happened to be lowest. That
# is legal but unreadable across a multi-panel figure: the rule exists so a
# reader can compare panels.
#
# Remapping is a RELABELLING ONLY. It permutes letter identities; it never
# changes which cells share a letter, so no statistical claim is affected.
#
# -----------------------------------------------------------------------------
# WHY LETTERS MUST NOT BE TRUNCATED BEFORE THIS RUNS  (defect fixed 2026-08-16)
# -----------------------------------------------------------------------------
# .mg_read_cld_tp() used to cut .group to its first two characters as a display
# tidy-up, BEFORE the remap. With six cells, three-letter groups are common and
# the cut silently changed the statistics being displayed. Real example, IID at
# interval level:
#
#     true groups     a | b | c | abc | abc | abc
#     after cut       a | b | c | ab  | ab  | ab      <-- "c" destroyed
#
# The three "abc" cells genuinely do NOT differ from the "c" cell; after
# truncation they no longer share a letter with it, so the panel rendered
# non-different cells as though they might differ. NND was affected the same
# way (a | ab | abc | bc | bc | c).
#
# It also broke label merging downstream: .mg_cld_label_rows() merges the two
# treatments at an interval when their letter STRINGS match, so truncating
# "abc" to "ab" made matching cells look mismatched and stacked two labels where
# one belonged. The manual cld_nudges workaround on the IID panel was written to
# hide exactly that symptom.
#
# Letters are therefore kept in full. Three characters is not a legibility
# problem; a wrong significance claim is.
# =============================================================================

.CLD_TREAT_ORDER <- c("control", "exercise choice")

# Relabel so the reading order control/1, control/2, control/3, exercise/1, ...
# introduces letters a, b, c, ... in that order.
.cld_remap_control_first <- function(cld_df,
                                     tp_col      = "tp_num",
                                     treat_col   = "treatment",
                                     group_col   = ".group",
                                     treat_order = .CLD_TREAT_ORDER) {
  if (is.null(cld_df) || !nrow(cld_df)) return(cld_df)
  if (!all(c(tp_col, treat_col, group_col) %in% names(cld_df))) return(cld_df)

  ord <- order(
    match(trimws(tolower(as.character(cld_df[[treat_col]]))), treat_order),
    suppressWarnings(as.numeric(as.character(cld_df[[tp_col]])))
  )

  seen <- character(0)
  for (i in ord) {
    lets <- strsplit(trimws(as.character(cld_df[[group_col]][i])), "")[[1]]
    for (l in lets) if (nzchar(l) && !(l %in% seen)) seen <- c(seen, l)
  }
  if (!length(seen)) return(cld_df)
  if (length(seen) > length(letters)) {
    warning("More CLD groups (", length(seen), ") than available letters; ",
            "leaving letters unchanged.")
    return(cld_df)
  }

  map <- stats::setNames(letters[seq_along(seen)], seen)
  cld_df[[group_col]] <- vapply(as.character(cld_df[[group_col]]), function(g) {
    lets <- strsplit(trimws(g), "")[[1]]
    r    <- map[lets]
    r    <- r[!is.na(r)]
    if (!length(r)) return(trimws(g))
    paste(sort(unique(r)), collapse = "")
  }, character(1L), USE.NAMES = FALSE)

  cld_df
}

# -----------------------------------------------------------------------------
# HOUSE RULE 2: no letters at all when nothing differs
# -----------------------------------------------------------------------------
# If every cell carries the same single group, no pair differs and the letters
# convey nothing -- six identical "a"s invite the reader to hunt for a contrast
# that is not there. Panels in that state are drawn WITHOUT a compact letter
# display, and their caption says so.
.cld_is_uniform <- function(cld_df, group_col = ".group") {
  if (is.null(cld_df) || !nrow(cld_df) || !group_col %in% names(cld_df)) return(FALSE)
  g <- unique(trimws(as.character(cld_df[[group_col]])))
  length(g) == 1L
}

# -----------------------------------------------------------------------------
# HOUSE RULE 3: no redundant letters  (added 2026-08-21)
# -----------------------------------------------------------------------------
# multcomp's insert-and-absorb algorithm produces a display that is CORRECT but
# not always MINIMAL: a cell can end up carrying a letter that conveys nothing,
# because every non-difference the letter would express is already expressed by
# another letter on the same cell. Real example -- bout_rate_flow:
#
#     control  5-25  ab      exercise choice  5-25  ad
#     control 45-65  c       exercise choice 45-65  abd   <-- "a" is redundant
#     control 85-105 b       exercise choice 85-105 d
#
# The "a" on exercise-choice/45-65 links it to control/5-25 and to
# exercise-choice/5-25. Both links are already carried -- by "b" and by "d"
# respectively -- so removing "a" changes no claim the panel makes. MV flagged
# it on sight ("U dont need the a in the exercise choice middle trial").
#
# WHAT THIS DOES. A compact letter display is valid iff, for every pair of
# cells, the cells share at least one letter EXACTLY WHEN the pair is not
# significantly different. Deleting a letter can only remove shared letters, so
# it can never turn a "differs" pair into a "same" pair -- only the non-different
# pairs need re-checking. The reduction therefore walks every letter occurrence
# in a fixed order, deletes it whenever the display stays valid without it, and
# repeats to a fixed point. The order is deterministic (cells in the canonical
# reading order, letters alphabetically) so the result is reproducible.
#
# SAFETY. If the ORIGINAL display is already inconsistent with the supplied
# p-value matrix -- which would mean the letters and the p-values came from
# different adjustments -- nothing is removed and a warning is issued. Reducing
# against a p-matrix that does not match the letters would be worse than leaving
# a redundant letter in place. A cell is never reduced to zero letters.
#
# Run this BEFORE .cld_remap_control_first(): removing a letter can empty out a
# letter identity entirely, and the remap is what closes the resulting gap in
# the a, b, c, ... sequence.

# Is `groups` (named character vector, cell id -> letter string) a valid display
# for the "same" relation encoded in `same` (logical matrix, cells x cells)?
.cld_display_valid <- function(groups, same) {
  ids <- names(groups)
  lets <- lapply(groups, function(g) strsplit(trimws(g), "")[[1]])
  for (i in seq_along(ids)) for (j in seq_along(ids)) {
    if (j <= i) next
    shares <- length(intersect(lets[[i]], lets[[j]])) > 0L
    if (shares != isTRUE(same[i, j])) return(FALSE)
  }
  TRUE
}

# `pmat`: symmetric numeric matrix of adjusted p-values, dimnames = cell ids,
# in the same order as `groups`. `alpha`: significance threshold.
.cld_drop_redundant_letters <- function(groups, pmat, alpha = 0.05, label = "") {
  ids <- names(groups)
  if (length(ids) < 2L) return(groups)
  if (!identical(rownames(pmat), ids) || !identical(colnames(pmat), ids))
    stop("cld reduction: p-matrix dimnames do not match the cell ids")

  same <- pmat >= alpha
  diag(same) <- TRUE

  if (!.cld_display_valid(groups, same)) {
    warning("CLD reduction skipped", if (nzchar(label)) paste0(" [", label, "]") else "",
            ": the letters supplied are not consistent with the p-value matrix ",
            "(different adjustment?). Letters left exactly as produced.")
    return(groups)
  }

  repeat {
    dropped <- FALSE
    for (i in seq_along(ids)) {
      lets <- sort(unique(strsplit(trimws(groups[[i]]), "")[[1]]))
      if (length(lets) < 2L) next          # never empty a cell
      for (l in lets) {
        trial <- groups
        keep  <- setdiff(strsplit(trimws(trial[[i]]), "")[[1]], l)
        if (!length(keep)) next
        trial[[i]] <- paste(sort(unique(keep)), collapse = "")
        if (.cld_display_valid(trial, same)) {
          groups  <- trial
          dropped <- TRUE
          break                            # re-read this cell's letters
        }
      }
    }
    if (!dropped) break
  }
  groups
}

# Convenience wrapper for a CLD data frame plus a pairwise table.
#
# `pairs_df` must carry integer columns `i` and `j` -- 1-based indices into
# `cell_ids` -- and a `p.value` column. Indices rather than parsed contrast
# labels on purpose: emmeans renders a cell as "control timepoint_f1" when a
# factor level looks numeric and as "control 1" when it does not, so any label
# parser is one factor-level rename away from silently matching nothing. The
# caller gets the indices from coef() on the contrast object, which is exact.
.cld_reduce_df <- function(cld_df, cell_ids, pairs_df, group_col = ".group",
                           alpha = 0.05, label = "") {
  n <- length(cell_ids)
  if (n < 2L) return(cld_df)
  if (!all(c("i", "j", "p.value") %in% names(pairs_df)))
    stop("cld reduction: pairs_df needs columns i, j, p.value")

  pmat <- matrix(NA_real_, n, n, dimnames = list(cell_ids, cell_ids))
  diag(pmat) <- 1
  for (k in seq_len(nrow(pairs_df))) {
    a <- pairs_df$i[k]; b <- pairs_df$j[k]
    if (is.na(a) || is.na(b) || a < 1 || b < 1 || a > n || b > n) next
    pmat[a, b] <- pmat[b, a] <- as.numeric(pairs_df$p.value[k])
  }
  if (anyNA(pmat)) {
    warning("CLD reduction skipped", if (nzchar(label)) paste0(" [", label, "]") else "",
            ": ", sum(is.na(pmat)), " pairwise p-values could not be matched to cells.")
    return(cld_df)
  }

  groups <- stats::setNames(as.list(trimws(as.character(cld_df[[group_col]]))), cell_ids)
  before <- unlist(groups, use.names = FALSE)
  groups <- .cld_drop_redundant_letters(groups, pmat, alpha = alpha, label = label)
  after  <- unlist(groups, use.names = FALSE)
  if (!identical(before, after))
    message("CLD reduction", if (nzchar(label)) paste0(" [", label, "]") else "",
            ": ", paste(before, collapse = " "), "  ->  ", paste(after, collapse = " "))
  cld_df[[group_col]] <- after
  cld_df
}

# -----------------------------------------------------------------------------
# ALTERNATIVE LETTER ORDER: lone "a" first  (added 2026-08-21)
# -----------------------------------------------------------------------------
# House rule 1 puts "a" on control/interval 1, but cannot make it a SINGLE
# letter -- the number of letters on a cell is fixed by the significance
# pattern, not by the labelling (see HOUSE RULE 1 below for the worked
# example). In a panel where control/interval 1 belongs to two non-difference
# groups, that leaves the display with no plain "a" anywhere, which reads
# oddly: a compact letter display is normally scanned from a single "a".
#
# This alternative puts "a" on the group of the FIRST cell in the canonical
# reading order that carries exactly one letter -- for bout_rate_flow that is
# control/interval 2, the cell nearest control/interval 1 both in reading order
# and in plot distance. Remaining groups then take b, c, d ... in order of
# first appearance, as usual.
#
# Like every other relabelling here this is cosmetic: it permutes letter
# identities and never changes which cells share a letter, so no statistical
# claim moves. Use it per panel, not globally -- most panels already have a
# lone "a" in the natural order and do not need it.
.cld_remap_lone_first <- function(cld_df,
                                  tp_col      = "tp_num",
                                  treat_col   = "treatment",
                                  group_col   = ".group",
                                  treat_order = .CLD_TREAT_ORDER) {
  if (is.null(cld_df) || !nrow(cld_df)) return(cld_df)
  if (!all(c(tp_col, treat_col, group_col) %in% names(cld_df))) return(cld_df)

  ord <- order(
    match(trimws(tolower(as.character(cld_df[[treat_col]]))), treat_order),
    suppressWarnings(as.numeric(as.character(cld_df[[tp_col]])))
  )
  lets_of <- function(i) strsplit(trimws(as.character(cld_df[[group_col]][i])), "")[[1]]

  # the group that will become "a": the first single-letter cell in reading order
  target <- NULL
  for (i in ord) {
    l <- lets_of(i)
    if (length(l) == 1L) { target <- l[1]; break }
  }
  if (is.null(target)) return(.cld_remap_control_first(cld_df, tp_col, treat_col,
                                                       group_col, treat_order))

  seen <- c(target)
  for (i in ord) for (l in lets_of(i)) if (nzchar(l) && !(l %in% seen)) seen <- c(seen, l)
  if (length(seen) > length(letters)) {
    warning("More CLD groups (", length(seen), ") than available letters; ",
            "leaving letters unchanged.")
    return(cld_df)
  }

  map <- stats::setNames(letters[seq_along(seen)], seen)
  cld_df[[group_col]] <- vapply(as.character(cld_df[[group_col]]), function(g) {
    r <- map[strsplit(trimws(g), "")[[1]]]
    r <- r[!is.na(r)]
    if (!length(r)) return(trimws(g))
    paste(sort(unique(r)), collapse = "")
  }, character(1L), USE.NAMES = FALSE)

  cld_df
}

# -----------------------------------------------------------------------------
# HOUSE RULE 1: control / interval 1 carries "a"
# -----------------------------------------------------------------------------
# .cld_remap_control_first() guarantees its group BEGINS with "a". It cannot
# guarantee that group is "a" ALONE, and no relabelling can:
#
#   the number of letters on a cell = the number of maximal non-difference
#   groups it belongs to, which is a property of the significance pattern, not
#   of the labelling.
#
# Worked example -- school area. The only pair that genuinely differs is
# control/interval-2 vs exercise-choice/interval-3. With exactly one differing
# pair the display is forced: those two cells take distinct single letters, and
# every other cell -- control/interval-1 included -- must carry BOTH, because it
# fails to differ from either. Printing "a" there instead of "ab" would assert
# that control/interval-1 differs from control/interval-2, which the model does
# not support.
#
# So the rule is enforced as far as it is valid (first letter is always "a") and
# .cld_check_control_first() reports, rather than silently rewrites, the cases
# where more than one letter is required.

# Post-condition guard. Verifies the invariant the whole rule exists to create,
# so a silent regression in the CLD source surfaces as a warning at build time
# rather than as a wrong letter in a published figure.
.cld_check_control_first <- function(cld_df, label = "",
                                     tp_col    = "tp_num",
                                     treat_col = "treatment",
                                     group_col = ".group",
                                     mode      = "control_first") {
  # Under the lone-first order the invariant deliberately does not hold, so
  # checking it there would report a warning for the intended state.
  if (identical(mode, "lone_first")) return(invisible(TRUE))
  if (is.null(cld_df) || !nrow(cld_df)) return(invisible(TRUE))
  if (!all(c(tp_col, treat_col, group_col) %in% names(cld_df)))
    return(invisible(TRUE))
  tp <- suppressWarnings(as.numeric(as.character(cld_df[[tp_col]])))
  i  <- which(trimws(tolower(as.character(cld_df[[treat_col]]))) == "control" &
              tp == min(tp, na.rm = TRUE))
  if (!length(i)) return(invisible(TRUE))
  g <- trimws(as.character(cld_df[[group_col]][i[1]]))
  ok <- nzchar(g) && substr(g, 1, 1) == "a"
  if (!ok)
    warning("CLD policy violated", if (nzchar(label)) paste0(" [", label, "]") else "",
            ": control/interval 1 carries '", g, "', expected a group starting with 'a'.")
  # Not a violation -- see HOUSE RULE 1 above -- but worth surfacing, because it
  # is the one case where the printed label departs from the plain "a" the house
  # style leads a reader to expect.
  if (ok && nchar(g) > 1L)
    message("CLD note", if (nzchar(label)) paste0(" [", label, "]") else "",
            ": control/interval 1 carries '", g, "' (", nchar(g), " letters). ",
            "It belongs to that many non-difference groups, so a single letter ",
            "would misstate the result; not rewritten.")
  invisible(ok)
}
