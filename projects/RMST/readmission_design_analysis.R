# Historical pilot: frailtypack::readmission. Times use the published t.stop scale.
tau <- 365                         # days
budget <- 300 * tau                # expected patient-days for a future cohort
n_boot <- 1000L
workers <- 4L                      # 1 for serial execution
seed <- 20261117L
out <- file.path(getwd(), "readmission_design_results")

stopifnot(tau > 0, is.finite(tau), budget >= 2 * tau, is.finite(budget),
          n_boot >= 2L, workers >= 1L)
data("readmission", package = "frailtypack", envir = environment())
d <- readmission[order(readmission$id, readmission$t.stop), ]
last <- !duplicated(d$id, fromLast = TRUE)
stopifnot(!anyNA(d[c("id", "t.start", "t.stop", "event", "death")]),
          all(d$event == as.integer(!last)), all(d$death %in% 0:1),
          all(d$death[!last] == 0), all(d$t.stop > d$t.start))
p <- d[last, c("id", "t.stop", "death")]
admissions <- d[d$event == 1L, ]
second <- admissions[ave(admissions$t.stop, admissions$id, FUN = seq_along) == 2L, ]
t2 <- second$t.stop[match(p$id, second$id)]
# U = time to second readmission or death. Death remains observable after U.
p$u <- pmin(t2, p$t.stop, na.rm = TRUE)
p$u_event <- !is.na(t2) | p$death == 1L

# Exact integration weights for KM steps on [left,right), in days.
x <- seq(log(10), 0, length.out = 31)
cuts <- sort(unique(c(0, tau, pmin(p$u, tau), pmin(p$t.stop, tau))))
left <- head(cuts, -1L)
width <- diff(cuts)
cost_weight <- sapply(x / tau, function(a) {
  if (a == 0) width else exp(-a * left) * (-expm1(-a * width)) / a
})
# E(Z^2) = 2 integral S_U(t) integral_0^t 1/G_L(s) ds dt.
moment_weight <- sapply(x / tau, function(a) {
  if (a == 0) return(width * (2 * left + width))
  2 * (exp(a * left) * expm1(a * width) / a - width) / a
})
L <- colSums(cost_weight)
setup <- list(L = L, width = width, cost = cost_weight, moment = moment_weight,
              u = match(pmin(p$u, tau), cuts), d = match(pmin(p$t.stop, tau), cuts),
              ue = p$u_event & p$u <= tau, de = p$death == 1 & p$t.stop <= tau,
              budget = budget)

# KM assumes independent historical censoring. G_L is set by the future design.
fit <- function(i, s) {
  km <- function(time_bin, event) {
    k <- length(s$width) + 1L
    exits <- tabulate(time_bin[i], k)
    deaths <- tabulate(time_bin[i[event[i]]], k)
    risk <- rev(cumsum(rev(exits)))
    survival <- head(cumprod(1 - deaths / pmax(1, risk)), -1L)
    if (any(head(risk - exits, -1L) == 0 & survival > 0))
      stop("Insufficient follow-up through tau for this KM estimate.")
    survival
  }
  su <- km(s$u, s$ue)
  sd <- km(s$d, s$de)
  mu <- sum(su * s$width)
  v <- drop(crossprod(su, s$moment)) - mu^2
  if (any(v < -1e-7)) stop("Negative variance. Check the integration weights.")
  v <- pmax(0, v)
  cost <- cbind(planned = s$L, actual = drop(crossprod(sd, s$cost)))
  n <- floor(s$budget / cost)
  H <- cost * v
  risk <- v / n
  j <- apply(H, 2L, which.min)
  list(mu = mu, v = v, cost = cost, n = n, H = H, risk = risk, j = j)
}

point <- fit(seq_len(nrow(p)), setup)
set.seed(seed)
indices <- replicate(n_boot, sample.int(nrow(p), replace = TRUE))
bootstrap_chunk <- function(index, s, fit_one) {
  t(vapply(seq_len(ncol(index)), function(b) {
    est <- fit_one(index[, b], s)
    c(mu = est$mu, est$j)
  }, numeric(3L)))
}
workers <- min(workers, n_boot)
if (workers == 1L) {
  boot <- bootstrap_chunk(indices, setup, fit)
} else {
  # Independent R processes; no socket connections or worker-side random seeds.
  run_workers <- function() {
    jobs <- list()
    on.exit(for (job in jobs) if (job$is_alive()) job$kill(), add = TRUE)
    for (chunk in parallel::splitIndices(n_boot, workers)) {
      jobs[[length(jobs) + 1L]] <- callr::r_bg(bootstrap_chunk,
        args = list(index = indices[, chunk, drop = FALSE], s = setup, fit_one = fit))
    }
    do.call(rbind, lapply(jobs, function(job) { job$wait(); job$get_result() }))
  }
  boot <- run_workers()
}
colnames(boot) <- c("mu", "planned", "actual")
ci <- quantile(boot[, "mu"], c(.025, .975), names = FALSE)
cohort <- data.frame(patients = nrow(p), horizon_days = tau,
  readmissions_by_tau = sum(admissions$t.stop <= tau),
  patients_with_second_readmission_by_tau = sum(t2 <= tau, na.rm = TRUE),
  deaths_by_tau = sum(p$death == 1 & p$t.stop <= tau),
  censored_alive_before_tau = sum(p$death == 0 & p$t.stop < tau),
  restricted_mean_days = point$mu, bootstrap_lower_95 = ci[1], bootstrap_upper_95 = ci[2])

summarize_design <- function(basis) {
  risk <- point$risk[, basis]
  rules <- c(short = 1L, H_optimal = point$j[[basis]],
             finite_budget_optimal = which.min(risk), full = length(L))
  j <- sort(unique(unname(rules)))
  data.frame(cost_basis = basis,
    rule = vapply(j, function(k) paste(names(rules)[rules == k], collapse = "+"), ""),
    planned_days = L[j], expected_cost_days = point$cost[j, basis],
    budget_patient_days = budget, n = point$n[j, basis], predicted_se_days = sqrt(risk[j]),
    excess_mse_pct = 100 * (risk[j] / min(risk) - 1),
    bootstrap_H_selection_pct = vapply(j, function(k) 100 * mean(boot[, basis] == k), 0.0))
}
designs <- do.call(rbind, lapply(c("planned", "actual"), summarize_design))
dir.create(out, showWarnings = FALSE, recursive = TRUE)
write.csv(cohort, file.path(out, "cohort_estimate.csv"), row.names = FALSE)
write.csv(designs, file.path(out, "design_comparison.csv"), row.names = FALSE)
pdf(file.path(out, "design_criterion.pdf"), width = 7, height = 4.8)
matplot(L, point$H / tau^3, type = "l", lty = 1:2, lwd = 2,
        col = c("#204870", "#AD4D2E"), xlab = "Planned mean follow-up (days)",
        ylab = expression(hat(H)(L) / tau^3))
points(L[point$j], point$H[cbind(point$j, 1:2)] / tau^3, pch = 19,
       col = c("#204870", "#AD4D2E"))
legend("topright", c("Planned-duration cost", "Cost through death"), bty = "n",
       lty = 1:2, lwd = 2, col = c("#204870", "#AD4D2E"))
invisible(dev.off())
