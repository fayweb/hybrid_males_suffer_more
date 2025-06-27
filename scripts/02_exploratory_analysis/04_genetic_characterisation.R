# Load your data (assuming field_mice is already loaded)

# 1. Basic sample size and HI distribution
n_mice <- nrow(field_mice)
hi_summary <- summary(field_mice$HI)
hi_range <- range(field_mice$HI, na.rm = TRUE)

# 2. Analyze genotyping success (number of loci amplified per mouse)
# Count successfully amplified loci
loci_success <- table(field_mice$HI_NLoci)
loci_summary <- summary(field_mice$HI_NLoci)

# Calculate percentages for different amplification thresholds
n_10plus_loci <- sum(field_mice$HI_NLoci >= 10, na.rm = TRUE)
pct_10plus_loci <- round(n_10plus_loci / n_mice * 100, 1)

n_4plus_loci <- sum(field_mice$HI_NLoci >= 4, na.rm = TRUE)
pct_4plus_loci <- round(n_4plus_loci / n_mice * 100, 1)

# Get minimum loci amplified
min_loci <- min(field_mice$HI_NLoci, na.rm = TRUE)

# 3. Check for genotyping bias across HI range
# Create HI bins to check if amplification success varies with HI
field_mice$HI_bin <- cut(field_mice$HI,
                         breaks = seq(0, 1, by = 0.1),
                         include.lowest = TRUE)

# Average loci success by HI bin
bias_check <- aggregate(HI_NLoci ~ HI_bin,
                        data = field_mice,
                        FUN = function(x) c(mean = mean(x), n = length(x)))

# 4. Calculate expected heterozygosity
field_mice$hHe <- 2 * field_mice$HI * (1 - field_mice$HI)
hHe_summary <- summary(field_mice$hHe)

# 5. Count mice by subspecies category
n_domesticus <- sum(field_mice$HI < 0.1, na.rm = TRUE)  # mostly domesticus
n_musculus <- sum(field_mice$HI > 0.9, na.rm = TRUE)    # mostly musculus
n_hybrids <- sum(field_mice$HI >= 0.1 & field_mice$HI <= 0.9, na.rm = TRUE)

# Print results for methods section
cat("GENETIC CHARACTERISATION RESULTS\n")
cat("================================\n")
cat("Total mice genotyped:", n_mice, "\n")
cat("HI range:", hi_range[1], "-", hi_range[2], "\n")
cat("Median HI:", median(field_mice$HI, na.rm = TRUE), "\n\n")

cat("GENOTYPING SUCCESS:\n")
cat("Loci amplified per mouse:\n")
print(loci_success)
cat("\nMinimum loci amplified:", min_loci, "\n")
cat("Mice with ≥10 loci:", n_10plus_loci, "(", pct_10plus_loci, "%)\n")
cat("Mice with ≥4 loci:", n_4plus_loci, "(", pct_4plus_loci, "%)\n\n")

cat("SUBSPECIES DISTRIBUTION:\n")
cat("Predominantly domesticus (HI < 0.1):", n_domesticus, "\n")
cat("Predominantly musculus (HI > 0.9):", n_musculus, "\n")
cat("Hybrids (0.1 ≤ HI ≤ 0.9):", n_hybrids, "\n\n")

# 6. Create visualization to check for bias
library(ggplot2)
p <- ggplot(field_mice, aes(x = HI, y = HI_NLoci)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(x = "Hybrid Index",
       y = "Number of Loci Amplified",
       title = "Genotyping Success Across Hybrid Index") +
  theme_minimal()
print(p)

# 7. Statistical test for bias
# Test if amplification success correlates with HI
cor_test <- cor.test(field_mice$HI, field_mice$HI_NLoci)
cat("Correlation between HI and amplification success:\n")
cat("r =", round(cor_test$estimate, 3), ", p =", round(cor_test$p.value, 4), "\n")

