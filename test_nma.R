# test_nma.R
#
# Stat-core smoke test for the FREQUENTIST network meta-analysis path of app.R.
#
# The app's primary engine is Bayesian (rjags + the external JAGS binary),
# which cannot be exercised in this environment. However, the app ALSO builds a
# frequentist netmeta object on the Network branch via netmeta::netmeta(), using
# the same (TE, seTE, treat1, treat2, studlab) contract. This test runs a tiny,
# known 3-treatment network through netmeta() and asserts structural properties
# with stopifnot(). It does NOT touch JAGS/rjags.
#
# Run: & "C:\Program Files\R\R-4.6.0\bin\Rscript.exe" test_nma.R

suppressMessages(library(netmeta))

# Tiny known network mirroring the app's network_example_smd sample:
# treatments A, B, C; contrasts A-B, B-C, A-C (fully connected triangle).
net_df <- data.frame(
  studlab = c("S1", "S2", "S3", "S4", "S5", "S6"),
  treat1  = c("A", "A", "B", "B", "A", "A"),
  treat2  = c("B", "B", "C", "C", "C", "C"),
  TE      = c(0.15, 0.22, 0.18, 0.30, 0.10, 0.20),
  seTE    = c(0.11, 0.10, 0.15, 0.14, 0.13, 0.14),
  stringsAsFactors = FALSE
)

treatments <- unique(c(net_df$treat1, net_df$treat2))

nm <- netmeta(
  TE = TE, seTE = seTE,
  treat1 = treat1, treat2 = treat2,
  studlab = studlab, data = net_df,
  sm = "SMD",
  common = TRUE, random = FALSE,
  reference.group = treatments[1]
)

# --- Property checks -------------------------------------------------------
# 1) Object class
stopifnot(inherits(nm, "netmeta"))

# 2) Treatment / study counts match the input network
stopifnot(nm$n == 3)              # 3 treatments (A, B, C)
stopifnot(nm$m == nrow(net_df))   # 6 pairwise comparisons

# 3) Fixed-effect treatment-effect matrix is square over the treatments,
#    skew-symmetric (TE.A.B == -TE.B.A), and zero on the diagonal.
TEmat <- nm$TE.common
stopifnot(nrow(TEmat) == 3, ncol(TEmat) == 3)
stopifnot(all(abs(diag(TEmat)) < 1e-8))
stopifnot(max(abs(TEmat + t(TEmat))) < 1e-8)

# 4) All pooled effects are finite (network is connected -> estimable).
stopifnot(all(is.finite(TEmat)))

# 5) netrank (used by the app's "Netrank Table" plot) returns a ranking
#    covering all 3 treatments with P-scores in [0, 1].
nr <- netrank(nm)
pscores <- nr$ranking.common
stopifnot(length(pscores) == 3)
stopifnot(all(pscores >= 0 - 1e-9 & pscores <= 1 + 1e-9))

# 6) netsplit (app's "Netsplit Forest Plot") runs without error on this network.
ns <- netsplit(nm)
stopifnot(inherits(ns, "netsplit"))

cat("test_nma.R: all netmeta smoke-test checks PASSED\n")
cat(sprintf("  treatments=%d comparisons=%d\n", nm$n, nm$m))
cat("  P-scores (common):", paste(round(pscores, 3), collapse = ", "), "\n")
