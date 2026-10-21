# ============================================================================
# YALE ED: PRIOR UTILIZATION, TRIAGE ACUITY AND INSURANCE STATUS
# Single-script, encounter-level analytical epidemiology study
# Data: 5v_cleandf.RData from Hong, Haimovich & Taylor (PLOS ONE, 2018)
# Research aims:
#   1. Association of prior-year ED utilization with current ED -> admission.
#   2. Effect modification by Emergency Severity Index (ESI).
#   3. Effect modification by insurance status.
# Analyses: descriptive tables, nonlinear adjusted logistic regression,
# interaction tests, marginal standardized predictions and contrasts,
# and sensitivity adjustment for prior-year admissions.
# IMPORTANT: De-identified public release has NO patient ID; some encounters
# come from the same patient. Model-based SEs, CIs and P values do NOT account
# for within-patient clustering and can be anti-conservative.
# The observational results are ASSOCIATIONS, not causal effects.
# ============================================================================

set.seed(20261021)
DATA_FILE <- "5v_cleandf.RData"       # Place this beside this script.
OUT_DIR <- "yale_ed_results"
VISITS_FOR_PREDICTION <- c(0, 1, 3, 6)
MAX_STANDARDIZATION_SAMPLE <- 20000L # Per ESI/insurance stratum, for speed.

if (!file.exists(DATA_FILE)) {
  stop("Cannot find ", DATA_FILE, ". Set working directory to the folder containing the RData file.")
}
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load only needed columns -----------------------------------------
# The full RData contains 560,486 x 972 variables and can use several GB RAM.
e <- new.env(parent = emptyenv())
objects <- load(DATA_FILE, envir = e)
if (!"df" %in% objects || !is.data.frame(e$df)) {
  stop("Expected the original Yale dataset to be a data.frame named 'df'. Found: ",
       paste(objects, collapse = ", "))
}
required <- c("disposition", "n_edvisits", "esi", "insurance_status", "age")
missing_required <- setdiff(required, names(e$df))
if (length(missing_required)) stop("Required columns missing: ", paste(missing_required, collapse = ", "))
optional <- c("n_admissions", "gender", "race", "ethnicity", "arrivalmode",
              "dep_name", "arrivalmonth", "chfnonhp", "diabmelnoc",
              "chrkidneydisease", "copd", "maligneopls")
d <- e$df[, intersect(c(required, optional), names(e$df)), drop = FALSE]
rm(e)
invisible(gc())

# Record missingness before recoding.
missingness <- data.frame(
  variable = names(d),
  missing_n = vapply(d, function(x) sum(is.na(x) | trimws(as.character(x)) == ""), numeric(1)),
  total_n = nrow(d)
)
missingness$missing_pct <- 100 * missingness$missing_n / missingness$total_n
write.csv(missingness, file.path(OUT_DIR, "00_missingness.csv"), row.names = FALSE)

# ---- 2. Construct the analysis dataset ----------------------------------
d$admitted <- ifelse(as.character(d$disposition) == "Admit", 1L,
                     ifelse(as.character(d$disposition) == "Discharge", 0L, NA_integer_))
d$prior_ed <- suppressWarnings(as.numeric(as.character(d$n_edvisits)))
d$age <- suppressWarnings(as.numeric(as.character(d$age)))
d$esi <- trimws(as.character(d$esi))

before_n <- nrow(d)
d <- d[!is.na(d$admitted) & is.finite(d$prior_ed) & d$prior_ed >= 0 &
       d$esi %in% as.character(1:5) & (is.na(d$age) | (d$age >= 18 & d$age <= 110)), , drop = FALSE]
if (nrow(d) < 1000) stop("Unexpectedly few valid encounters; check the variable coding.")
d$esi <- factor(d$esi, levels = as.character(1:5))
d$esi <- stats::relevel(d$esi, ref = "3")
d$log_ed <- log1p(d$prior_ed)

# Use a pragmatic missingness indicator for age rather than excluding every
# encounter with an unknown age; categorical missingness is made explicit.
d$age_missing <- as.integer(is.na(d$age))
d$age[is.na(d$age)] <- median(d$age, na.rm = TRUE)

# Reduce sparse categorical levels for stable main-effects/interaction models.
# Returns common observed levels plus "Other" (and possibly "Unknown").
group_cat <- function(x, keep = 5L, min_n = 1000L) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- "Unknown"
  freq <- sort(table(x), decreasing = TRUE)
  selected <- head(names(freq)[freq >= min_n], keep)
  if (!length(selected)) selected <- head(names(freq), keep)
  y <- ifelse(x %in% selected, x, "Other")
  factor(y, levels = unique(c(selected, if (any(y == "Other")) "Other")))
}

# Mandatory insurance exposure; keep the five biggest categories and pool rest.
d$insurance_grp <- group_cat(d$insurance_status, keep = 5L, min_n = 1000L)
if (nlevels(d$insurance_grp) < 2L) stop("Not enough distinct insurance groups.")

covars <- character(0)
for (nm in c("gender", "race", "ethnicity", "arrivalmode", "dep_name", "arrivalmonth")) {
  if (nm %in% names(d)) {
    newname <- paste0(nm, "_grp")
    d[[newname]] <- group_cat(d[[nm]], keep = if (nm == "arrivalmonth") 12L else 6L,
                              min_n = if (nm == "arrivalmonth") 1L else 500L)
    if (nlevels(d[[newname]]) > 1L) covars <- c(covars, newname)
  }
}

# A small set of recorded comorbidity-history indicators, as possible
# clinical case-mix adjustments. These are EHR history flags, not diagnoses
# newly made during the current ED encounter.
for (nm in c("chfnonhp", "diabmelnoc", "chrkidneydisease", "copd", "maligneopls")) {
  if (nm %in% names(d)) {
    v <- suppressWarnings(as.numeric(as.character(d[[nm]])))
    v[!is.finite(v)] <- 0
    newname <- paste0("hx_", nm)
    d[[newname]] <- as.integer(v > 0)
    if (length(unique(d[[newname]])) > 1L) covars <- c(covars, newname)
  }
}
if (length(unique(d$age_missing)) > 1L) covars <- c(covars, "age_missing")

# Keep only model variables to save memory; all output tables are derived here.
keep <- unique(c("admitted", "prior_ed", "log_ed", "esi", "insurance_grp", "age",
                 "n_admissions", covars))
d <- d[, intersect(keep, names(d)), drop = FALSE]
invisible(gc())

meta <- data.frame(
  measure = c("encounters_in_file", "encounters_analyzed", "encounters_excluded",
              "admission_fraction", "prior_ED_median", "prior_ED_95th_percentile"),
  value = c(before_n, nrow(d), before_n - nrow(d),
            mean(d$admitted), median(d$prior_ed),
            as.numeric(quantile(d$prior_ed, 0.95)))
)
write.csv(meta, file.path(OUT_DIR, "01_sample_summary.csv"), row.names = FALSE)
write.csv(data.frame(insurance_group = levels(d$insurance_grp),
                     n = as.integer(table(d$insurance_grp))),
          file.path(OUT_DIR, "02_insurance_groups.csv"), row.names = FALSE)

# ---- 3. Descriptive prevalence -------------------------------------------
d$prior_visit_bin <- cut(d$prior_ed, breaks = c(-Inf, 0, 1, 2, 4, 9, Inf),
                          labels = c("0", "1", "2", "3-4", "5-9", "10+"))
crude_summary <- function(group_vars) {
  g <- interaction(d[, group_vars, drop = FALSE], drop = TRUE, sep = " | ")
  ns <- as.integer(table(g))
  events <- as.numeric(tapply(d$admitted, g, sum))
  frac <- events / ns
  # Wilson binomial CI: descriptive ONLY, does not correct patient clustering.
  z <- stats::qnorm(.975)
  den <- 1 + z^2 / ns
  ctr <- (frac + z^2 / (2 * ns)) / den
  half <- z * sqrt(frac * (1 - frac) / ns + z^2 / (4 * ns^2)) / den
  data.frame(group = names(table(g)), n = ns, admitted = events,
             rate = frac, lower = pmax(0, ctr - half), upper = pmin(1, ctr + half))
}
write.csv(crude_summary(c("esi", "prior_visit_bin")),
          file.path(OUT_DIR, "03_crude_ESI_by_prior_visits.csv"), row.names = FALSE)
write.csv(crude_summary(c("insurance_grp", "prior_visit_bin")),
          file.path(OUT_DIR, "04_crude_insurance_by_prior_visits.csv"), row.names = FALSE)

# ---- 4. Logistic models: nonlinear visit history and interactions --------
# No prior admissions in PRIMARY models: admissions are a subset of previous
# ED encounters. Adding them while hypothetically setting previous visits to
# zero would create impossible predictor combinations. Prior admissions are
# included in a separately labeled sensitivity analysis below.
base_terms <- c("splines::ns(age, df = 4)", "esi", "insurance_grp", covars)
base_formula <- stats::as.formula(paste("admitted ~", paste(base_terms, collapse = " + ")))
main_formula <- update(base_formula, . ~ . + splines::ns(log_ed, df = 3))
esi_formula <- update(main_formula, . ~ . + splines::ns(log_ed, df = 3):esi)
ins_formula <- update(main_formula, . ~ . + splines::ns(log_ed, df = 3):insurance_grp)
full_formula <- update(main_formula, . ~ . +
  splines::ns(log_ed, df = 3):esi + splines::ns(log_ed, df = 3):insurance_grp)

cat("Fitting logistic models on", nrow(d), "encounters...\n")
m0 <- glm(base_formula, data = d, family = binomial())
m1 <- glm(main_formula, data = d, family = binomial())
m2 <- glm(esi_formula, data = d, family = binomial())
m3 <- glm(ins_formula, data = d, family = binomial())
m4 <- glm(full_formula, data = d, family = binomial())

compare_lrt <- function(smaller, bigger, hypothesis) {
  ll0 <- as.numeric(logLik(smaller)); ll1 <- as.numeric(logLik(bigger))
  df0 <- attr(logLik(smaller), "df"); df1 <- attr(logLik(bigger), "df")
  stat <- max(0, 2 * (ll1 - ll0)); df_diff <- df1 - df0
  data.frame(hypothesis = hypothesis, chi_sq = stat, df = df_diff,
             p_value = if (df_diff > 0) pchisq(stat, df_diff, lower.tail = FALSE) else NA_real_)
}
lrt <- do.call(rbind, list(
  compare_lrt(m0, m1, "Prior ED visits add information"),
  compare_lrt(m1, m2, "Prior ED visits x ESI (without insurance interaction)"),
  compare_lrt(m1, m3, "Prior ED visits x insurance (without ESI interaction)"),
  compare_lrt(m3, m4, "Prior ED visits x ESI conditional on insurance interaction"),
  compare_lrt(m2, m4, "Prior ED visits x insurance conditional on ESI interaction"),
  compare_lrt(m1, m4, "Both interaction terms together")
))
write.csv(lrt, file.path(OUT_DIR, "05_interaction_tests_EXPLORATORY.csv"), row.names = FALSE)
write.csv(data.frame(model = c("Baseline", "ED history", "ED history x ESI",
                                "ED history x insurance", "Both interactions"),
                     AIC = c(AIC(m0), AIC(m1), AIC(m2), AIC(m3), AIC(m4))),
          file.path(OUT_DIR, "06_model_AIC.csv"), row.names = FALSE)
capture.output(summary(m4), file = file.path(OUT_DIR, "07_full_logistic_model.txt"))
if (anyNA(coef(m4))) {
  stop("Full model has aliased coefficients. Inspect insurance/ESI cell sizes or collapse sparse groups.")
}
rm(m0, m1, m2, m3)
invisible(gc())

# ---- 5. Adjusted predicted admission probabilities -----------------------
# We average predictions within each observed ESI/insurance stratum, while
# holding its observed case mix fixed. This is observational standardization,
# NOT an intervention or a causal effect. Subsample up to 20,000 encounters
# per stratum for computational feasibility (same sampled patients per contrast).
# CI uses delta method with GLM model-based vcov, NO patient clustering.
V <- vcov(m4)
b <- coef(m4)
T4 <- stats::delete.response(stats::terms(m4))

one_prediction <- function(subset_df, number_prior_visits) {
  nd <- subset_df
  nd$log_ed <- log1p(number_prior_visits)
  X <- model.matrix(T4, data = nd, contrasts.arg = m4$contrasts, xlev = m4$xlevels)
  X <- X[, names(b), drop = FALSE]
  p <- plogis(drop(X %*% b))
  gradient <- colMeans(X * as.vector(p * (1 - p)))
  list(probability = mean(p), gradient = gradient)
}

std_rows <- list(); contrast_rows <- list()
for (stratum_name in c("esi", "insurance_grp")) {
  for (level in levels(d[[stratum_name]])) {
    idx <- which(d[[stratum_name]] == level)
    if (length(idx) < 500L) next
    if (length(idx) > MAX_STANDARDIZATION_SAMPLE) idx <- sample(idx, MAX_STANDARDIZATION_SAMPLE)
    sub <- d[idx, , drop = FALSE]
    predictions <- lapply(VISITS_FOR_PREDICTION, function(v) one_prediction(sub, v))
    for (j in seq_along(VISITS_FOR_PREDICTION)) {
      p <- predictions[[j]]$probability
      grad <- predictions[[j]]$gradient
      se <- sqrt(max(0, as.numeric(crossprod(grad, V %*% grad))))
      std_rows[[length(std_rows) + 1L]] <- data.frame(
        group_type = stratum_name, group = as.character(level),
        sample_n = nrow(sub), prior_ED_visits = VISITS_FOR_PREDICTION[j],
        adjusted_probability = p,
        lower95 = pmax(0, p - 1.96 * se), upper95 = pmin(1, p + 1.96 * se))
    }
    # Risk difference: predicted admission probability with 3 vs 0 past visits.
    a <- predictions[[which(VISITS_FOR_PREDICTION == 3)]]
    z <- predictions[[which(VISITS_FOR_PREDICTION == 0)]]
    dif <- a$probability - z$probability
    grad_dif <- a$gradient - z$gradient
    se_dif <- sqrt(max(0, as.numeric(crossprod(grad_dif, V %*% grad_dif))))
    contrast_rows[[length(contrast_rows) + 1L]] <- data.frame(
      group_type = stratum_name, group = as.character(level),
      risk_difference_3_vs_0 = dif,
      lower95 = dif - 1.96 * se_dif, upper95 = dif + 1.96 * se_dif)
  }
}
std <- do.call(rbind, std_rows)
contrasts <- do.call(rbind, contrast_rows)
write.csv(std, file.path(OUT_DIR, "08_standardized_probabilities.csv"), row.names = FALSE)
write.csv(contrasts, file.path(OUT_DIR, "09_standardized_risk_differences.csv"), row.names = FALSE)

# ---- 6. Sensitivity: add previous hospital admissions --------------------
if ("n_admissions" %in% names(d)) {
  d$prior_admissions <- suppressWarnings(as.numeric(as.character(d$n_admissions)))
  usable <- is.finite(d$prior_admissions) & d$prior_admissions >= 0 &
            d$prior_admissions <= d$prior_ed
  write.csv(data.frame(n_total = nrow(d), n_consistent_prior_admissions = sum(usable)),
            file.path(OUT_DIR, "10_prior_admissions_consistency.csv"), row.names = FALSE)
  if (sum(usable) > 1000 && length(unique(d$prior_admissions[usable])) >= 2L) {
    ds <- d[usable, , drop = FALSE]
    ds$log_prev_admissions <- log1p(ds$prior_admissions)
    ms <- glm(update(full_formula, . ~ . + log_prev_admissions),
              data = ds, family = binomial())
    capture.output(summary(ms),
                   file = file.path(OUT_DIR, "11_prior_admissions_adjusted_sensitivity.txt"))
    rm(ms, ds)
  }
}

# ---- 7. Simple presentation plots (base R; no extra packages) -----------
plot_standardized <- function(which_group, title, filename) {
  zz <- std[std$group_type == which_group, , drop = FALSE]
  grps <- unique(zz$group)
  png(filename, width = 1250, height = 820, res = 140)
  par(mar = c(5.0, 5.5, 3.2, 1.5))
  plot(NA, xlim = range(VISITS_FOR_PREDICTION), ylim = range(c(zz$lower95, zz$upper95)),
       xlab = "Number of prior-year ED visits", ylab = "Adjusted probability of admission",
       main = title, xaxt = "n")
  axis(1, at = VISITS_FOR_PREDICTION)
  grid(col = "gray90")
  colors <- grDevices::hcl.colors(length(grps), palette = "Dark 3")
  for (i in seq_along(grps)) {
    q <- zz[zz$group == grps[i], , drop = FALSE]
    q <- q[order(q$prior_ED_visits), , drop = FALSE]
    lines(q$prior_ED_visits, q$adjusted_probability, type = "b", lwd = 2,
          pch = 19, col = colors[i])
  }
  legend("topleft", legend = grps, col = colors, lty = 1, pch = 19,
         cex = .8, bty = "n")
  dev.off()
}
plot_standardized("esi", "ED utilization and admission by ESI",
                  file.path(OUT_DIR, "12_ESI_adjusted_probabilities.png"))
plot_standardized("insurance_grp", "ED utilization and admission by insurance",
                  file.path(OUT_DIR, "13_insurance_adjusted_probabilities.png"))

# ---- 8. Record interpretation/limitations --------------------------------
writeLines(c(
  "Study: Prior emergency department utilization and hospital admission:",
  "       effect modification by ESI triage severity and insurance status.",
  "", "READ FIRST:",
  "1. Unit of analysis = ED encounter, NOT unique patient.",
  "2. Patient IDs were removed in the public dataset. Repeat encounters are",
  "   not independent; logistic model CIs/P-values can be too narrow/small.",
  "3. This is an observational ASSOCIATION study, not a causal analysis.",
  "4. Insurance categories are as recorded and collapsed for adequate counts.",
  "5. ESI 1 = most urgent, ESI 5 = least urgent.",
  "6. Prior ED visits counts visits in the preceding year; model applies log1p.",
  "7. For adjusted probabilities, up to 20,000 observed encounters per stratum",
  "   were used for case-mix standardization; uncertainty is model-based only.",
  "8. '3 versus 0 visits' contrasts are standardized associations, not effects",
  "   of changing a person's ED use. Significant interactions alone do not",
  "   establish clinical importance or healthcare inequity.",
  "9. The prior-admissions-adjusted model is a sensitivity analysis, because",
  "   previous admissions are a subset of previous ED visits.",
  "10. No within-person longitudinal outcomes can be reconstructed here.",
  "", "Core results:",
  "03/04: observed admission proportions; 05: interaction tests;",
  "08: adjusted admission probabilities; 09: adjusted risk differences.",
  "", "Original data: Hong WS, Haimovich AD, Taylor RA (2018), PLOS ONE",
  "https://doi.org/10.1371/journal.pone.0201016"
), con = file.path(OUT_DIR, "README_interpretation.txt"))

cat("\nFinished. All tables/plots saved in:", normalizePath(OUT_DIR), "\n")
cat("Main outputs: 05_interaction_tests_EXPLORATORY.csv, ",
    "08_standardized_probabilities.csv, 09_standardized_risk_differences.csv\n", sep = "")
