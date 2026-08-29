# Refit all fourteen reported behavioural outcomes from the DEPOSITED CSV alone
# and check each treatment F/df against the published value. This is the claim
# the deposit README makes, so it has to be tested against the deposited file
# rather than against the pipeline.
suppressPackageStartupMessages({
  library(readr); library(lme4); library(lmerTest); library(pbkrtest)
})
d <- readr::read_csv(file.path(PROJECT_ROOT, "zenodo_dataset/pref_trials.csv"),
                     show_col_types = FALSE)
d$treatment  <- factor(d$treatment, levels = c("control", "exercise choice"))
d$interval_f <- factor(d$interval)
d$phys_trial <- factor(d$phys_trial)
cat("deposited rows:", nrow(d), " cols:", ncol(d), "\n\n")

tr <- function(y, how) switch(how, none = y, log1p = log1p(y), sqrt = sqrt(y), stop("?"))

SPEC <- list(
  list("alr_flow",              "none",  NA,             17.503896,  14.000000),
  list("alr_high",              "none",  NA,              7.529668,  14.000000),
  list("alr_medium",            "none",  NA,              5.275043,  14.000000),
  list("alr_low",               "none",  NA,              2.178357,  14.000000),
  list("crossings_per_session", "log1p", NA,             45.003611,  14.000000),
  list("bouts_per_min",         "none",  NA,             59.884535,  14.000000),
  list("longest_flow_bout_s",   "none",  "fish_density", 86.675038,  13.000000),
  list("longest_calm_bout_s",   "none",  "fish_density",  3.419979,  13.000000),
  list("mean_flow_bout_s",      "none",  "fish_density",  0.585329,  13.357711),
  list("mean_calm_bout_s",      "none",  "fish_density",  0.107205,  12.902739),
  list("mean_nnd_cm",           "log1p", NA,             10.697117,  14.000000),
  list("mean_iid_cm",           "sqrt",  NA,              3.192528,  14.000000),
  list("mean_school_area_cm2",  "sqrt",  NA,              6.480284,  14.000000),
  list("mean_school_speed_cm_s","log1p", NA,              0.745289,  14.000000)
)

ok <- 0L
cat(sprintf("%-24s %5s %12s %12s %9s %9s  %s\n",
            "outcome", "n", "F_refit", "F_published", "df2", "p", "verdict"))
for (sp in SPEC) {
  col <- sp[[1]]; how <- sp[[2]]; extra <- sp[[3]]; Fpub <- sp[[4]]; df2pub <- sp[[5]]
  dd <- d[!is.na(d[[col]]), ]
  dd$y <- tr(as.numeric(dd[[col]]), how)
  rhs <- "treatment * interval_f"
  if (!is.na(extra)) rhs <- paste(rhs, "+", extra)
  fit <- lmerTest::lmer(stats::as.formula(paste("y ~", rhs, "+ (1 | phys_trial)")),
                        data = dd, REML = TRUE)
  av <- anova(fit, ddf = "Kenward-Roger")
  Fr <- av["treatment", "F value"]; df2 <- av["treatment", "DenDF"]
  pr <- av["treatment", "Pr(>F)"]
  good <- abs(Fr - Fpub) < 1e-3 && abs(df2 - df2pub) < 1e-3
  ok <- ok + good
  cat(sprintf("%-24s %5d %12.6f %12.6f %9.4f %9.2g  %s\n",
              col, nrow(dd), Fr, Fpub, df2, pr, if (good) "MATCH" else "*** MISMATCH ***"))
}
cat(sprintf("\n%d / %d reported behavioural outcomes reproduced from pref_trials.csv alone\n",
            ok, length(SPEC)))
if (ok != length(SPEC)) quit(status = 1)
