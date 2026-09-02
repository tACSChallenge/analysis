#  Bayesian analyses for H1 and H2, tACS Challenge
#
#  The method used here is not preregistered. See the bayes_readme.txt 
#
#  INPUT
#    bayes_H1H2_delta.csv, written by tACS_exportBayesData.m

library(BayesFactor)

# Edit to point at the CSV written by tACS_exportBayesData.m
csv_path <- "bayes_H1H2_delta.csv"

if (!file.exists(csv_path)) {
  stop(sprintf("File not found: %s\nRun tACS_exportBayesData.m first.", csv_path))
}

data <- read.csv(csv_path)
data$delta <- data$Amplitude - data$PermMean

conditions  <- c("A", "B", "C")
cond_labels <- c(A = "Occipital", B = "Retinal", C = "Cutaneous")


#  H1 

cat("\n - H1: One-sample Bayes Factor (JZS, r = sqrt(2)/2) -\n\n")

for (cond in conditions) {
  d <- data$delta[data$Condition == cond]
  d <- d[!is.na(d)]

  bf   <- ttestBF(x = d, mu = 0)
  bf10 <- as.numeric(extractBF(bf)$bf)
  bf01 <- 1 / bf10
  tt   <- t.test(d, mu = 0)

  cat(sprintf("%-10s | N = %3d | mean delta = %+.4f | t = %+.2f | p = %.3f\n",
              cond_labels[[cond]], length(d), mean(d),
              unname(tt$statistic), tt$p.value))
  cat(sprintf("           | BF10 = %.4g | BF01 = %.4g\n\n", bf10, bf01))
}


#  H2
cat("\n - H2: paired bayes factor, occipital vs control -\n\n")

wide <- reshape(data[, c("SubjectID", "Condition", "delta")],
                idvar = "SubjectID", timevar = "Condition", direction = "wide")

for (ctrl in c("B", "C")) {
  x <- wide[["delta.A"]]
  y <- wide[[paste0("delta.", ctrl)]]
  keep <- !is.na(x) & !is.na(y)

  bf   <- ttestBF(x = x[keep], y = y[keep], paired = TRUE)
  bf10 <- as.numeric(extractBF(bf)$bf)
  bf01 <- 1 / bf10
  tt   <- t.test(x[keep], y[keep], paired = TRUE)

  cat(sprintf("Occipital vs %-10s | N = %3d | mean difference = %+.4f | t = %+.2f\n",
              cond_labels[[ctrl]], sum(keep), mean(x[keep] - y[keep]),
              unname(tt$statistic)))
  cat(sprintf("                        | BF10 = %.4g | BF01 = %.4g\n\n", bf10, bf01))
}


#  Prior robustness

r_values <- c(0.1, 0.2, 0.3, 0.4, 0.5, sqrt(2)/2, 1.0, 1.2, 1.5)

cat("\n Prior Robustness, H1\n")

robustness <- data.frame()

for (cond in conditions) {
  d <- data$delta[data$Condition == cond]
  d <- d[!is.na(d)]

  cat(sprintf("\n--- %s ---\n", cond_labels[[cond]]))
  cat(sprintf("%-10s %-10s\n", "r", "BF01"))

  for (r in r_values) {
    bf   <- ttestBF(x = d, mu = 0, rscale = r)
    bf01 <- 1 / as.numeric(extractBF(bf)$bf)
    tag  <- if (abs(r - sqrt(2)/2) < 1e-6) "  <- default" else ""
    cat(sprintf("%-10.3f %-10.3f%s\n", r, bf01, tag))
    robustness <- rbind(robustness,
                        data.frame(Condition = cond_labels[[cond]], r = r, BF01 = bf01))
  }
}

write.csv(robustness, "bayes_H1_prior_robustness.csv", row.names = FALSE)

